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

from .codon_dfa import sanitize_sequences
from .codon_usage import reverse_translate_most_common
from .pldesign import ParalinearDesign
from . import CODON_USAGE_FILE as DEFAULT_CODONUSAGE
import pandas as pd
import argparse
import sys

def parse_arguments():
    parser = argparse.ArgumentParser(description='ParalinearDesign sequence optimization tool')

    # Argument groups
    input_output = parser.add_argument_group('Input/Output Options')
    optimization = parser.add_argument_group('Optimization Parameters')
    sampling = parser.add_argument_group('Random Sampling Options')
    adaptive = parser.add_argument_group('Adaptive Sampling Options')
    performance = parser.add_argument_group('Performance Options')

    # Input/output options
    input_output.add_argument('--protein', '-p', type=str, required=True,
                        help='Protein sequence')
    input_output.add_argument('--utr3', '-u', type=str, default='',
                        help='3\' UTR sequence')
    input_output.add_argument('--output', '-o', type=argparse.FileType('w'), default=sys.stdout,
                        help='Output file (default: stdout)')
    input_output.add_argument('--fasta', action='store_true',
                        help='Output in FASTA format')
    input_output.add_argument('--fasta-cds-only', action='store_true',
                        help='Output only CDS in FASTA format')
    input_output.add_argument('--allow-no-stop', action='store_true',
                        help='Allow sequences with no stop codon')
    input_output.add_argument('--skip-vienna-reevaluation', action='store_true',
                        help='Don\'t run ViennaRNA for folding reevaluation')
    input_output.add_argument('--quiet', '-q', action='store_true',
                        help='Suppress progress output')

    # Optimization parameters
    optimization.add_argument('--lambda', '-l', dest='ldlambda', type=float, default=1.0,
                        help='Lambda parameter for LinearDesign (default: 1.0)')
    optimization.add_argument('--gamma', '-g', dest='caigamma', type=float, default=1.0,
                        help='Gamma for CAI weight transformation (default: 1.0)')
    optimization.add_argument('--codonusagefile', '-c', type=argparse.FileType('r'), default=DEFAULT_CODONUSAGE,
                        help='LinearDesign codon usage file')
    optimization.add_argument('--minimum-cai-weight', type=float, default=0.0,
                        help='Don\'t use codons under this CAI weight (default: 0.0)')
    optimization.add_argument('--minimum-rscu', type=float, default=0.0,
                        help='Don\'t use codons under this RSCU (default: 0.0)')
    optimization.add_argument('--mask', '-m', action='append', default=[],
                        help='Mask regions in format "start5,end5,start3,end3"')
    optimization.add_argument('--exclude', action='append', default=[],
                        help='Exclude sequence patterns in optimized sequences '
                             '(e.g., "GAATTC")')
    optimization.add_argument('--exclude-start', type=int, default=0,
                        help='Exclude N residues from the beginning from optimization')

    # Random sampling options
    sampling.add_argument('--num-random-samples', '-r', type=int, default=1,
                        help='Number of random samples to generate')
    sampling.add_argument('--random-rejection-rate', type=float, default=0.01,
                        help='Rejection rate for random sampling (default: 0.01)')
    sampling.add_argument('--random-seed-start', type=int, default=0,
                        help='Starting random seed (default: 0)')
    sampling.add_argument('--dont-retry-on-failure', action='store_true',
                        help='Don\'t retry on LinearDesign error')

    # Adaptive sampling options
    adaptive.add_argument('--adaptive-sampling', type=int, default=0,
                        help='Number of adaptive sampling iterations (default: 0)')
    adaptive.add_argument('--adaptive-sampling-max-gap-size', type=int, default=2,
                        help='Maximum gap size for finding stems in adaptive masking (default: 2)')
    adaptive.add_argument('--adaptive-sampling-min-group-size', type=int, default=5,
                        help='Minimum size for groups in adaptive masking (default: 5)')
    adaptive.add_argument('--adaptive-sampling-max-group-size', type=int, default=8,
                        help='Maximum size for groups in adaptive masking (default: 8)')
    adaptive.add_argument('--adaptive-sampling-min-mask-regions', type=float, default=0.1,
                        help='Minimum fraction of mask regions (default: 0.1)')
    adaptive.add_argument('--adaptive-sampling-max-mask-regions', type=float, default=0.3,
                        help='Maximum fraction of mask regions (default: 0.3)')
    adaptive.add_argument('--adaptive-sampling-turnover-rate', type=float, default=0.3,
                        help='Turnover rate of masks in subsequent iterations (default: 0.3)')

    # Performance options
    performance.add_argument('--threads', '-t', type=int, default=1,
                        help='Number of threads to use (default: 1)')

    # Debugging options
    parser.add_argument('--prepare-ld-debug', action='store_true',
                        help=argparse.SUPPRESS)

    args = parser.parse_args()
    masks = ParalinearDesign.parse_masks(args.mask)

    if masks and args.adaptive_sampling > 0:
        raise ValueError('Adaptive sampling and manual masking cannot be used '
                         'simultaneously.')

    if (args.num_random_samples > 1 and args.adaptive_sampling > 0 and
            args.random_rejection_rate < 0):
        raise ValueError('Random rejection rate must be set with a non-negative '
                         'probability when using multiple random samples.')

    return args, masks


def main():
    args, masks = parse_arguments()

    codonusage_filename = args.codonusagefile.name

    protein = args.protein.upper()
    utr3 = args.utr3.upper().replace('T', 'U')
    protein, utr3 = sanitize_sequences(protein, utr3, args.allow_no_stop)

    exclude_start = args.exclude_start
    removed_nterm = removed_nterm_cds = ''
    if exclude_start > 0:
        removed_nterm = protein[:exclude_start]
        removed_nterm_cds = reverse_translate_most_common(removed_nterm, codonusage_filename)
        protein = protein[exclude_start:]

    adaptive_sampling_params = None
    if args.adaptive_sampling > 0:
        adaptive_sampling_params = {
            'iterations': args.adaptive_sampling,
            'max_gap_size': args.adaptive_sampling_max_gap_size,
            'min_group_size': args.adaptive_sampling_min_group_size,
            'max_group_size': args.adaptive_sampling_max_group_size,
            'min_mask_regions': args.adaptive_sampling_min_mask_regions,
            'max_mask_regions': args.adaptive_sampling_max_mask_regions,
            'turnover_rate': args.adaptive_sampling_turnover_rate,
        }

    debug_opts = {
        'prepare_ld_debug': args.prepare_ld_debug,
    }

    with ParalinearDesign(protein, utr3, args.ldlambda, args.caigamma, codonusage_filename,
                          args.exclude, args.minimum_cai_weight, args.minimum_rscu,
                          args.allow_no_stop, args.skip_vienna_reevaluation,
                          removed_nterm_cds, args.num_random_samples > 1 or args.quiet,
                          debug_opts) as ldenv:

        if args.num_random_samples == 1:
            r = ldenv.run(masks, args.random_rejection_rate if args.random_seed_start != 0 else 0,
                          args.random_seed_start, validate=True)
            res = [r]
        else:
            res = ldenv.run_parallel_sampling(args.num_random_samples, args.threads,
                                              masks, args.random_rejection_rate,
                                              args.random_seed_start,
                                              adaptive_sampling_params,
                                              validate=True,
                                              dont_retry_on_failure=args.dont_retry_on_failure,
                                              quiet=args.quiet)

    if args.fasta or args.fasta_cds_only:
        if args.fasta_cds_only:
            seqoutput_slice = slice(0, len(protein) * 3 + len(removed_nterm_cds))
        else:
            seqoutput_slice = slice(None)

        for i, r in enumerate(res):
            print(f'>opt_{i+1} mfe={r["mfe"]} cai={r["cai"]}', file=args.output)
            print(r['sequence'][seqoutput_slice], file=args.output)
    else:
        output_columns = ['sequence', 'structure', 'mfe', 'cai']
        pd.DataFrame(res)[output_columns].to_csv(args.output, index=False)

    print('Optimization finished.', file=sys.stderr)

if __name__ == '__main__':
    main()
