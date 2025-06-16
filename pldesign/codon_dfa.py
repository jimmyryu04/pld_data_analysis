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

from .codon_usage import get_codon_usage
import pandas as pd
import numpy as np
import sys
from collections import defaultdict

PROTEIN_ALPHABET = 'ACDEFGHIKLMNPQRSTVWY'
STOP = '*'
STOP_CODONS = ['UAA', 'UAG', 'UGA']
SCORE_MASKING = -1e100


class CodonDFABuilder:

    # The orders in coding wheel of LinearDesign
    DFA_ORDER1 = 'CUAG'
    DFA_ORDER2 = 'CUAG'
    DFA_ORDER3 = 'GCUA'

    def __init__(self, codonusage, weight, utr3_topology_compat=False):
        self.codon_usage = codonusage
        self.weight = weight
        self.utr3_topology_compat = utr3_topology_compat

        self.codon_graphs = self.build_codon_graphs()

    def build_codon_graphs(self):
        tbl = self.codon_usage.codon_table.copy()
        tbl['codon1'] = tbl['codon'].str[0]
        tbl['codon2'] = tbl['codon'].str[:2]

        # Sort the table by the order of the coding wheel
        tbl['c1key'] = tbl['codon1'].map(lambda x: self.DFA_ORDER1.index(x[-1]))
        tbl['c2key'] = tbl['codon2'].map(lambda x: self.DFA_ORDER2.index(x[-1]))
        tbl['c3key'] = tbl['codon'].map(lambda x: self.DFA_ORDER3.index(x[-1]))

        tbl = tbl.sort_values(by=['c1key', 'c2key', 'c3key'])
        tbl['freq'] = tbl['freq'].clip(lower=1e-30) # Avoid division by zero

        return {
            aa: self.build_codon_graph_aa(codons)
            for aa, codons in tbl.groupby('aa')}

    def build_codon_graph_aa(self, codons):
        # Calculate the weights for each codon and intermediate state
        c1freqs = (
            codons.drop_duplicates(subset='codon1', keep='first')[['codon1', 'freq']]
                    .rename(columns={'codon1': 'codon', 'freq': 'freq'})
                    .reset_index(drop=True))
        c2freqs = (
            codons.drop_duplicates(subset='codon2', keep='first')[['codon2', 'freq']]
                    .rename(columns={'codon2': 'codon', 'freq': 'freq'})
                    .reset_index(drop=True))
        c3freqs = codons[['codon', 'freq']].reset_index(drop=True)
        startpoint = pd.DataFrame([{'codon': '', 'freq': c1freqs['freq'].max()}])

        cfreqs = pd.concat([startpoint, c1freqs, c2freqs, c3freqs]).reset_index()
        cfreqs['level'] = cfreqs['codon'].str.len()
        cfreqs['parent'] = cfreqs['codon'].str[:-1]
        cfreqs.loc[0, 'parent'] = '-'
        freqmax = cfreqs.groupby('parent')['freq'].max()

        cfreqs['weight'] = (np.log(cfreqs['freq'] / cfreqs['parent'].map(freqmax))
                            * self.weight * 100).round(3)

        # Build the node index
        cfreqs['nodeidx'] = cfreqs.apply(lambda x: x['index'], axis=1)
        nodeidx = cfreqs.set_index('codon')['nodeidx'].to_dict()

        # Build the edges
        edges = defaultdict(list)
        for _, crow in cfreqs.iterrows():
            if crow['level'] == 3:
                continue

            state = crow['codon']
            children = cfreqs[cfreqs['parent'] == crow['codon']]
            for _, child in children.iterrows():
                emission = child['codon'][-1]
                weight = child['weight']
                edges[state].append((emission, weight,
                                     child['index'] if child['level'] < 3 else 0))

        return nodeidx, dict(edges)

    def generate_single_codon(self, aa, start_idx, output, fixed_triplet=None):
        if output is None:
            output = defaultdict(list)

        nodeidx, edges = self.codon_graphs[aa]
        nodes = ['']

        # Generate the DFA
        for i in range(3):
            nextnodes = []

            for node in nodes:
                nidx = (start_idx + i, nodeidx[node])

                for emission, weight, cidx in edges[node]:
                    if fixed_triplet is not None:
                        if fixed_triplet[i] != emission:
                            weight = SCORE_MASKING
                        else:
                            weight = 0
                    output[nidx].append((emission, weight, cidx))
                    nextnodes.append(node + emission)

            nodes = nextnodes

        return output

    def generate_untranslated_dfa(self, triplet, start_idx, output):
        if output is None:
            output = defaultdict(list)

        if self.utr3_topology_compat:
            aa = self.codon_usage.codon_to_aa[triplet]
            self.generate_single_codon(aa, start_idx, output, triplet)
        else:
            # Generate the DFA
            for i in range(3):
                output[start_idx + i, 0].append((triplet[i], 0.0, 0))

        return output

    def generate_dfa(self, protseq, utr3seq):
        # Remove the trailing UTR sequence that does not fit into a codon
        utr3seq = utr3seq[:len(utr3seq) // 3 * 3]

        # Generate the DFA for the protein sequence
        output = None
        for i, aa in enumerate(protseq):
            output = self.generate_single_codon(aa, i * 3, output)

        # Generate the DFA for the UTR sequence
        utrprotseq = []
        for i in range(0, len(utr3seq), 3):
            triplet = utr3seq[i:i+3]
            concated_pos = i + len(protseq) * 3
            utrprotseq.append(self.codon_usage.codon_to_aa[triplet])
            output = self.generate_untranslated_dfa(triplet, concated_pos, output)

        utrprotseq = ''.join(utrprotseq)

        return output, protseq + utrprotseq, utr3seq

    def format_dfa_lineardesign(self, dfa):
        # Generate the output
        output = []

        for nodeid in sorted(dfa.keys()):
            output.append('-')
            output.append('N {0} {1}'.format(*nodeid))
            for emission, weight, ptr in dfa[nodeid]:
                output.append(f'R {emission} {weight} {ptr}')

        return '\n'.join(output)

    builder_cache = {}
    @classmethod
    def get_builder(kls, codonusagefile, ldlambda, utr3_topology_compat=False):
        cache_key = (codonusagefile, ldlambda, utr3_topology_compat)
        if cache_key in CodonDFABuilder.builder_cache:
            return kls.builder_cache[cache_key]

        usage = get_codon_usage(codonusagefile)
        builder = CodonDFABuilder(usage, ldlambda, utr3_topology_compat)
        kls.builder_cache[cache_key] = builder
        return builder

def sanitize_sequences(protein, utr3, allow_no_stop=False):
    # Check if the UTR sequence has other letters than A, C, G, U
    if not all(c in 'ACGU' for c in utr3):
        raise ValueError('UTR sequence must only contain A, C, G, U')

    # Check if the protein sequence contains only valid amino acids
    if not all(c in PROTEIN_ALPHABET or c == '*' for c in protein):
        raise ValueError('Protein sequence contains invalid amino acids')

    if not allow_no_stop:
        # Check if the sequences contain stop codons
        if not protein.endswith(STOP) and utr3[:3] not in STOP_CODONS:
            print('>> Protein sequence does not contain a stop codon. Adding '
                  'a stop codon to the end of the protein sequence.',
                  file=sys.stderr)

            protein += STOP

    return protein, utr3

def generate_codon_dfa_file(codonusagefile, ldlambda, protein, utr3, output,
                            allow_no_stop, dfa_modifiers=[],
                            utr3_topology_compat=False):
    protein, utr3 = sanitize_sequences(protein, utr3, allow_no_stop)

    dfabuilder = CodonDFABuilder.get_builder(codonusagefile, ldlambda,
                                             utr3_topology_compat)

    dfa, dfaseq, utr3seq_fixed = dfabuilder.generate_dfa(protein, utr3)
    for modifier in dfa_modifiers:
        modifier(dfa, dfaseq, utr3seq_fixed)
    formatted_dfa = dfabuilder.format_dfa_lineardesign(dfa)

    print(dfaseq, file=output)
    print(utr3seq_fixed, file=output)
    print(formatted_dfa, file=output)

def mask_codon_in_dfa(dfa, aapos, codon, weight=-1e100):
    assert type(aapos) is int

    # Follow the graph to find the last edge for the codon
    c1_edges = dfa[aapos * 3, 0]

    c2_nodeidx = next(nidx for base, _, nidx in c1_edges if base == codon[0])
    c2_edges = dfa[aapos * 3 + 1, c2_nodeidx]

    c3_nodeidx = next(nidx for base, _, nidx in c2_edges if base == codon[1])
    c3_edges = dfa[aapos * 3 + 2, c3_nodeidx]

    # Mask the codon by changing the weight of the last edge
    for i in range(len(c3_edges)):
        base, orig_weight, nidx = c3_edges[i]
        if base == codon[2]:
            c3_edges[i] = (base, weight, nidx)
            break

    # Check if the DFA is still valid
    for j in (0, 1):
        edges = dfa.get((aapos * 3 + 2, j), ())
        if any(w > weight for base, w, nidx in edges):
            break
    else:
        # Restore the original weight
        c3_edges[i] = (base, orig_weight, nidx)
        raise ValueError('Cannot mask codon - no feasible coding sequence remaining.')

if __name__ == '__main__':
    import click

    @click.command()
    @click.option('--codonusagefile', '-c', type=click.File('r'),
                default='codon_usage_freq_table_human.csv')
    @click.option('--lambda', '-l', 'ldlambda', type=float, default=1.0)
    @click.option('--protein', '-p', type=str, required=True)
    @click.option('--utr3', '-u', type=str, default='')
    @click.option('--output', '-o', type=click.File('w'), default='-')
    @click.option('--allow-no-stop', is_flag=True, default=False)
    @click.option('--utr3-topology-compat', is_flag=True, default=False)
    def main(codonusagefile, ldlambda, protein, utr3, output, allow_no_stop,
            utr3_topology_compat):
        generate_codon_dfa_file(codonusagefile, ldlambda, protein, utr3, output,
                                allow_no_stop, [], utr3_topology_compat)

    main()
