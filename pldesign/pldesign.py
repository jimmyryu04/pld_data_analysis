#!/usr/bin/env python3
#
# Copyright (c) 2025 Seoul National University
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#

from .codon_dfa import generate_codon_dfa_file, mask_codon_in_dfa
from .codon_usage import get_codon_usage
from .adaptive_masking import find_masking_regions
from .resite import RestrictionSiteExclusion
from . import PLDESIGN_BINARY, PLDESIGN_DIR, CODON_USAGE_FILE
from fcntl import fcntl, F_GETFL, F_SETFL
from itertools import count
from concurrent import futures
from select import select
from io import StringIO
import tempfile
import random
import pandas as pd
import numpy as np
import subprocess as sp
from tqdm import tqdm
from Bio import Seq
from ViennaRNA import fold
import sys
import os

UTR3_FULL_TOPOLOGY = False # Set to True to use full CDS DFA topology in 3' UTR
PARALINEARDESIGN_CAPABILITY = 2

DEBUG_SCRIPT_TEMPLATE = """\
#!/bin/bash
export _PLD_REQUEST=1
DFA=<<-END
{dfa}
END

cd "{cwd}"
echo "$DFA" | {command}

exit $?
"""

class SequenceValidationError(Exception):
    pass

class RestrictedPatternFoundError(SequenceValidationError):
    @property
    def matches(self):
        return self.args[1]

def validate_optimization_result(protein, utr3, rnaseq, rnastr, exclude_patterns=[]):

    cdslen = len(protein) * 3
    utrlen = len(utr3)

    translated = Seq.Seq(rnaseq[:cdslen]).translate()
    if translated != protein:
        raise SequenceValidationError(
                'The optimized sequence does not match the protein sequence.')

    if utr3 and rnaseq[-utrlen:].upper() != utr3:
        raise SequenceValidationError(
                'The optimized sequence does not match the 3\' UTR sequence.')

    if cdslen + utrlen - len(rnaseq) > 3: # allow 3 nt difference for stop codon
        raise SequenceValidationError(
                'The optimized sequence is shorter than expected.')

    if rnastr.count('(') != rnastr.count(')'):
        raise SequenceValidationError(
                'The optimized sequence has unmatched base pairings.')

    # Check the validity of basepairings
    fivep_positions = []
    for i, c in enumerate(rnastr):
        if c == '(':
            fivep_positions.append(i)
            continue
        elif c != ')':
            continue

        p5 = fivep_positions.pop()
        pairnucs = ''.join(sorted((rnaseq[p5] + rnaseq[i]).upper()))
        if pairnucs not in ('AU', 'CG', 'GU'):
            raise SequenceValidationError(
                    f'Invalid basepairing at positions {p5+1} and {i+1}.')

    # Check the existence of sequence patterns to avoid
    if exclude_patterns:
        resite_matcher = RestrictionSiteExclusion(exclude_patterns)
        matches = list(resite_matcher.find_sites(rnaseq))
        if matches:
            raise RestrictedPatternFoundError(
                    'Found a sequence pattern in the optimized sequence.',
                    matches)

def trim_utr3_into_frame(seq):
    taillen = len(seq) % 3
    return (seq, '') if taillen == 0 else (seq[:-taillen], seq[-taillen:])

class WorkloadAdjustedProgressBar:

    def __init__(self, total, **kwargs):
        self.loadmax_pct = self.workload_func(total) / 100
        self.total = total
        self.pbar = None
        self.kwargs = kwargs

        self.prev_n = 0
        self.prev_pct = 0.0

    def workload_func(self, x):
        return x ** 2

    def close(self):
        self.pbar.close()

    def update(self, i):
        if self.pbar is None:
            self.pbar = tqdm(total=100, unit='%', **self.kwargs)

        progress = self.prev_n + i
        progress_pct = self.workload_func(progress) / self.loadmax_pct
        self.prev_n = progress

        self.pbar.n = round(progress_pct, 2)
        self.pbar.refresh()


class ParalinearDesign:

    MAX_RETRIES_RESITE_REMOVAL = 10
    RESITE_REMOVAL_MASKING_WEIGHT = -200000000

    def __init__(self, protein, utr3, ldlambda, caigamma, codonusage_filename,
                 exclude_patterns=[], minimum_cai_weight=0, minimum_rscu=0,
                 allow_no_stop=False, skip_vienna_reevaluation=False,
                 preamble='', quiet=False, debug_opts={}):
        self.protein = protein
        self.utr3 = utr3
        self.ldlambda = ldlambda
        self.caigamma = caigamma
        self.codonusage_filename = codonusage_filename
        self.exclude_patterns = exclude_patterns
        self.minimum_cai_weight = minimum_cai_weight
        self.minimum_rscu = minimum_rscu
        self.allow_no_stop = allow_no_stop
        self.skip_vienna_reevaluation = skip_vienna_reevaluation
        self.quiet = quiet
        self.codonusage_tmpfile = None
        self.debug_opts = debug_opts
        self.preamble = preamble

        if self.codonusage_filename is None:
            self.codonusage_filename = CODON_USAGE_FILE

    def __enter__(self):
        self.prepare_env()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if self.codonusage_tmpfile is not None:
            if os.path.exists(self.codonusage_tmpfile.name):
                os.unlink(self.codonusage_tmpfile.name)
            self.codonusage_tmpfile.close()

    def prepare_env(self):
        utr3 = self.utr3.replace('T', 'U').upper().strip()
        utr3_trimmed, utr3_tail = trim_utr3_into_frame(utr3)

        requires_modified_codonusage = (
            abs(self.caigamma - 1) > 1e-6 or
            self.minimum_cai_weight > 0 or self.minimum_rscu > 0)

        if requires_modified_codonusage:
            codonusage_filename = self.prepare_gamma_transformed_codonusage()
        else:
            codonusage_filename = self.codonusage_filename

        envvars = os.environ.copy()
        envvars['_PLD_REQUEST'] = '1'

        self.run_params = {
            'envvars': envvars,
            'utr3': utr3,
            'utr3_trimmed': utr3_trimmed,
            'utr3_tail': utr3_tail,
            'codonusage_filename': codonusage_filename,
        }

        self.generate_dfa(self.run_params, [])

        # For the CAI calculation, we need to use the original codon usage table
        self.codon_usage = get_codon_usage(self.codonusage_filename)

    def generate_dfa(self, params, remove_paths, rng=random):
        modifiers = []
        if remove_paths:
            def remove_paths_from_dfa(dfa, dfaseq, utr3seq,
                                      weight=self.RESITE_REMOVAL_MASKING_WEIGHT):
                for path in remove_paths:
                    if path['full_codons']:
                        selected = [rng.choice(path['full_codons'])]
                    else:
                        selected = [rng.choice(path['partial_codons'])]

                    for aapos, codon in selected:
                        mask_codon_in_dfa(dfa, aapos, codon, weight)

            modifiers.append(remove_paths_from_dfa)

        dfafile = StringIO()

        generate_codon_dfa_file(
            open(params['codonusage_filename']), self.ldlambda, self.protein,
            params['utr3_trimmed'], dfafile, self.allow_no_stop, modifiers,
            utr3_topology_compat=UTR3_FULL_TOPOLOGY)

        params['dfa'] = dfafile.getvalue().encode()

    def prepare_gamma_transformed_codonusage(self):
        freqtbl = pd.read_csv(self.codonusage_filename,
                              names=['codon', 'aa', 'freq'], comment='#')

        # Mask rare codons upon request
        maxfreq = freqtbl.groupby('aa')['freq'].max()
        meanfreq = freqtbl.groupby('aa')['freq'].mean()
        rscu = freqtbl['freq'] / freqtbl['aa'].map(meanfreq)
        cai_w = freqtbl['freq'] / freqtbl['aa'].map(maxfreq)
        tomask = (rscu < self.minimum_rscu) | (cai_w < self.minimum_cai_weight)
        if tomask.any():
            freqtbl.loc[tomask, 'freq'] = 1e-100

        # Apply gamma transformation and re-normalize
        maxfreq = freqtbl.groupby('aa')['freq'].max()
        freqtbl['freq_g'] = (freqtbl['freq'] / freqtbl['aa'].map(maxfreq)) ** self.caigamma
        aacaisum = freqtbl.groupby('aa')['freq_g'].sum()
        freqtbl['freq_g'] = (freqtbl['freq_g'] / freqtbl['aa'].map(aacaisum)).round(6)

        # Write the transformed codon usage table to a temporary file
        tmpfile = tempfile.NamedTemporaryFile(mode='w', delete=False)
        print('#,,', file=tmpfile)
        freqtbl[['codon', 'aa', 'freq_g']].to_csv(tmpfile, header=False, index=False)
        tmpfile.flush()
        tmpfile.seek(0)

        self.codonusage_tmpfile = tmpfile
        return tmpfile.name

    def run(self, masks, random_rejection_rate=0.0, random_seed=0, validate=True):
        params = self.run_params
        if self.exclude_patterns:
            validate=True

        rng = random.Random(random_seed)

        paths_to_remove = []
        for i in range(self.MAX_RETRIES_RESITE_REMOVAL):
            try:
                return self._run(masks, params, random_rejection_rate, random_seed, validate)
            except RestrictedPatternFoundError as exc:
                paths_to_remove.extend(exc.matches)
                if i == 0:
                    params = params.copy() # avoid modifying the class variable

                self.generate_dfa(params, paths_to_remove, rng=rng)

        return {
            'error': 'Failed to generate sequence without restricted patterns '
                     f'after {self.MAX_RETRIES_RESITE_REMOVAL} attempts',
            'restricted_sites': paths_to_remove,
            'seed': random_seed,
            'masks': masks,
        }

    def _run(self, masks, params, random_rejection_rate, random_seed, validate):
        verbose = '0'
        masks_formatted = ' '.join(f'{s5},{e5},{s3},{e3}' for (s5, e5), (s3, e3) in masks)

        ld_args = [PLDESIGN_BINARY, str(self.ldlambda), verbose,
                   self.codonusage_filename, masks_formatted,
                   str(random_seed), str(random_rejection_rate)]

        if self.debug_opts.get('prepare_ld_debug'):
            self.prepare_ld_debug(ld_args, params['dfa'])
            raise SystemExit

        pbar = WorkloadAdjustedProgressBar(
                    total=len(self.protein) * 3 + len(params['utr3_trimmed']),
                    file=sys.stderr, desc='ParalinearDesign', disable=self.quiet)
        j = 0
        ret = []

        RESULT_MARKER = 'mRNA sequence:  '
        if not self.quiet:
            print('=> Preparing to run ParalinearDesign optimization...', file=sys.stderr,
                  end='\r', flush=True)

        with sp.Popen(ld_args, stdin=sp.PIPE, stdout=sp.PIPE, env=params['envvars'],
                      cwd=PLDESIGN_DIR) as proc:
            # Send the sequence and DFA
            proc.stdin.write(params['dfa'])
            proc.stdin.close()

            # Check compability
            header = proc.stdout.readline()
            if header.startswith(b'PLD-capability'):
                capability = int(header.decode().split(':')[-1])
                if capability < PARALINEARDESIGN_CAPABILITY:
                    raise ValueError('Required ParalinearDesign capability is not '
                                     'available.')
            else:
                raise ValueError('Customized version of ParalinearDesign is required.')

            # Read the updates and pass to the progress bar
            for chunk in self.read_live_updates(proc.stdout):
                if chunk.startswith('j='):
                    newj = int(chunk[2:])
                    if newj - j >= 5:
                        pbar.update(newj - j)
                        j = newj
                elif chunk.startswith(RESULT_MARKER):
                    ret.append(chunk)

        pbar.update(pbar.total)
        pbar.close()
        ret = ''.join(ret)

        if RESULT_MARKER not in ret:
            raise ValueError('Unexpected output from ParalinearDesign.')

        lines = ret.split(RESULT_MARKER)[1].splitlines()
        rnaseq = lines[0].strip() + params['utr3_tail'].lower()
        rnaseq_cds = rnaseq[:len(self.protein) * 3]
        rnastr_ld = lines[1].split(': ')[1].strip() + '.' * len(params['utr3_tail'])

        if validate:
            validate_optimization_result(self.protein, params['utr3'], rnaseq, rnastr_ld,
                                         self.exclude_patterns)

        # CAI from LD can be inaccurate when caigamma != 1
        # cai = float(lines[2].split(': ')[2].strip())
        cai = round(self.codon_usage.calc_cai(rnaseq_cds), 3)

        if self.skip_vienna_reevaluation:
            rnastr = '.' * len(self.preamble) + rnastr_ld
            mfe = float(lines[2].split(': ')[1].split()[0])
        else:
            rnastr, mfe = fold(self.preamble + rnaseq)
            mfe = round(mfe, 3)

        return {
            'sequence': self.preamble + rnaseq,
            'structure': rnastr,
            'mfe': mfe,
            'cai': cai,
            'masks': masks,
        }

    def run_in_thread(self, masks, random_rejection_rate, seed, validate=True):
        try:
            return self.run(masks, random_rejection_rate, seed, validate=validate)
        except Exception as e:
            return {'error': e.args[0], 'seed': seed, 'masks': masks}

    def run_parallel_sampling(self, num_random_samples, threads, masks,
                              random_rejection_rate=0.0, random_seed_start=0,
                              adaptive_sampling_params=None,
                              validate=True, dont_retry_on_failure=False,
                              quiet=False):
        res = []
        jobs = []

        with futures.ThreadPoolExecutor(max_workers=threads) as executor:
            if adaptive_sampling_params is not None:
                return AdaptiveSamplingSession(self, executor, adaptive_sampling_params,
                            num_random_samples, random_rejection_rate, random_seed_start,
                            validate, dont_retry_on_failure, quiet)()

            def submit_job_for_next_seed():
                seed = next(gen_seed)
                job = executor.submit(
                    self.run_in_thread, masks, random_rejection_rate, seed, validate)
                jobs.append(job)

            def iter_results(jobs):
                while jobs:
                    done, _ = futures.wait(jobs, timeout=1,
                                           return_when=futures.FIRST_COMPLETED)
                    for job in done:
                        jobs.remove(job)
                        yield job.result()

            gen_seed = count(random_seed_start)

            for _ in range(num_random_samples):
                submit_job_for_next_seed()

            msg_will_try = '' if dont_retry_on_failure else ' Will retry with a new seed.'
            num_failures = 0
            with tqdm(total=num_random_samples, unit='sample', file=sys.stderr,
                      desc='Sampling', disable=quiet) as pbar:
                for r in iter_results(jobs):
                    pbar.update()
                    if 'error' in r:
                        num_failures += 1
                        print(f'(ok={len(res)}/fail={num_failures}) '
                              f'Failed for seed {r["seed"]}.' + msg_will_try,
                              f'[{r["error"]}]', file=sys.stderr)
                        if not dont_retry_on_failure:
                            pbar.total += 1
                            submit_job_for_next_seed()
                    else:
                        res.append(r)

        return res

    @staticmethod
    def parse_masks(args):
        if not args:
            return []

        masks = []
        for arg in args:
            start5, end5, start3, end3 = map(int, arg.split(','))
            masks.append(((start5, end5), (start3, end3)))

        return masks

    @staticmethod
    def read_live_updates(stdout, bufsize=8192):
        fd = stdout.fileno()
        flags = fcntl(fd, F_GETFL)
        fcntl(fd, F_SETFL, flags | os.O_NONBLOCK)

        buf = []

        while True:
            r, _, _ = select([fd], [], [], 0.1)
            if not r:
                continue

            data = os.read(fd, bufsize)
            if not data:
                break

            if b'\r' in data:
                chunks = data.split(b'\r')
                yield (b''.join(buf) + chunks[0]).decode()

                for chunk in chunks[1:-1]:
                    yield chunk.decode()

                buf[:] = [chunks[-1]]
            else:
                buf.append(data)

        if buf:
            yield b''.join(buf).decode()

    def prepare_ld_debug(self, ld_args, dfa):
        import shlex

        output_filename = 'pldesign_debug.sh'
        with open(output_filename, 'w') as f:
            command = ' '.join(map(shlex.quote, ld_args))
            dfa = dfa.decode()
            f.write(DEBUG_SCRIPT_TEMPLATE.format(command=command, dfa=dfa, cwd=PLDESIGN_DIR))

        print(f'Wrote debug script to {output_filename}.')


class AdaptiveSamplingSession:

    SAMPLER_OPTIONS = [
        'iterations', 'max_gap_size', 'min_group_size', 'max_group_size',
        'min_mask_regions', 'max_mask_regions', 'turnover_rate',
    ]

    def __init__(self, ldenv, executor, adaptive_sampling_params, num_random_samples,
                 random_rejection_rate, random_seed_start, validate,
                 dont_retry_on_failure, quiet):
        self.ldenv = ldenv
        self.num_random_samples = num_random_samples
        self.random_rejection_rate = random_rejection_rate
        self.random_seed_start = random_seed_start
        self.executor = executor
        self.validate = validate
        self.dont_retry_on_failure = dont_retry_on_failure
        self.quiet = quiet

        for key in self.SAMPLER_OPTIONS:
            setattr(self, key, adaptive_sampling_params[key])

    def initialize(self):
        self.jobs = []
        self.res = []
        self.gen_seed = count(self.random_seed_start)

        self.msg_will_try = (
            '' if self.dont_retry_on_failure else ' Will retry with a new seed.')
        self.num_failures = 0

    def submit_job_for_next_seed(self, masks):
        seed = next(self.gen_seed)
        job = self.executor.submit(
            self.ldenv.run_in_thread, masks, self.random_rejection_rate, seed, self.validate)
        self.jobs.append(job)

    def iter_results(self):
        while self.jobs:
            done, _ = futures.wait(self.jobs, timeout=1, return_when=futures.FIRST_COMPLETED)
            for job in done:
                self.jobs.remove(job)
                yield job.result()

    def msg(self, *args):
        if not self.quiet:
            print(*args, file=sys.stderr)

    def get_masking_regions(self, structure):
        return find_masking_regions(structure, self.max_gap_size,
                                    self.min_group_size, self.max_group_size)

    def __call__(self):
        self.initialize()

        self.msg('=> Evaluating the sequence without masks...')
        self.submit_job_for_next_seed([])
        res = self.wait_for_results()
        assert len(res) == 1

        for iter_no in range(1, self.iterations + 1):
            masking_regions = [(self.get_masking_regions(r['structure']), r['masks']) for r in res]

            self.msg(f'=> Adaptive samping iteration {iter_no}')
            self.res.clear()
            self.jobs.clear()

            for _ in range(self.num_random_samples):
                new_source_masks, old_source_masks = random.choice(masking_regions)
                maskindices = list(range(len(new_source_masks)))

                num_masks = int(np.random.uniform(self.min_mask_regions, self.max_mask_regions) *
                                len(new_source_masks))
                keeping_count = int(len(old_source_masks) * self.turnover_rate)
                keeping_count = max(0, min(keeping_count, num_masks - 1))
                turnover_count = num_masks - keeping_count

                selected_new_masks = np.random.choice(maskindices, turnover_count, replace=False)
                selected_new_masks = sorted(selected_new_masks)
                selected_new_masks = [new_source_masks[i] for i in selected_new_masks]

                if keeping_count > 0:
                    maskindices = list(range(len(old_source_masks)))
                    selected_old_masks = np.random.choice(maskindices, keeping_count,
                                                        replace=False)
                    selected_old_masks = sorted(selected_old_masks)
                    selected_old_masks = [old_source_masks[i] for i in selected_old_masks]
                else:
                    selected_old_masks = []

                selected_masks = sorted(selected_new_masks + selected_old_masks)

                self.submit_job_for_next_seed(selected_masks)

            res = self.wait_for_results()
            assert len(res) == self.num_random_samples

        return res

    def wait_for_results(self):
        with tqdm(total=len(self.jobs), unit='sample', file=sys.stderr,
                  desc='Sampling', disable=self.quiet) as pbar:
            for r in self.iter_results():
                pbar.update()
                if 'error' in r:
                    self.num_failures += 1
                    print(f'(ok={len(self.res)}/fail={self.num_failures}) '
                          f'Failed for seed {r["seed"]}.' + self.msg_will_try,
                          f'[{r["error"]}]', file=sys.stderr)
                    if not self.dont_retry_on_failure:
                        pbar.total += 1
                        self.submit_job_for_next_seed(r['masks'])
                else:
                    self.res.append(r)

        return self.res
