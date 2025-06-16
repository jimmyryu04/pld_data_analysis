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

from .codon_dfa import mask_codon_in_dfa
import re

IUPAC_PATTERNS = {
    'R': r'[AG]',
    'Y': r'[CU]',
    'S': r'[GC]',
    'W': r'[AU]',
    'K': r'[GU]',
    'M': r'[AC]',
    'B': r'[CGU]',
    'D': r'[AGU]',
    'H': r'[ACU]',
    'V': r'[ACG]',
    'N': r'[ACGU]'
}


class RestrictionSiteExclusion:

    def __init__(self, patterns):
        self.patterns = self.compile_patterns(patterns)

    def compile_patterns(self, patterns):
        substituted = []

        for pat in patterns:
            pat = pat.upper().replace('T', 'U')
            pat = ''.join([IUPAC_PATTERNS.get(base, base) for base in pat])
            substituted.append(pat)
        
        return re.compile('|'.join(substituted))

    def find_sites(self, seq):
        # seq = seq.upper().replace('T', 'U')

        for match in self.patterns.finditer(seq):
            start = match.start()
            end = match.end()

            # Interval including the partial matches
            overlap_start = start // 3
            overlap_end = (end + 2) // 3

            # Interval only including the full codons
            full_start = (start + 2) // 3
            full_end = end // 3

            full_codons = list(range(full_start, full_end))
            partial_codons = (
                list(range(overlap_start, full_start)) +
                list(range(full_end, overlap_end))
            )

            yield {
                'start': start,
                'end': end,
                'match': seq[start:end],
                'full_codons': [(pos, seq[pos * 3:pos * 3 + 3])
                                for pos in full_codons],
                'partial_codons': [(pos, seq[pos * 3:pos * 3 + 3])
                                   for pos in partial_codons],
            }


if __name__ == '__main__':
    from pldesign import codon_usage, codon_dfa
    from Bio import Seq

    seq = (
        'AUGGUGUCCAAGGGCGAGGAGCUCUUCACAGGAGUCGUGCCCAUUCUGGUGGAGCUGGAUGGCGACGUGAAC'
        'GGGCACAAGUUCUCCGUGAGCGGGGAGGGCGAGGGUGACGCCACGUACGGGAAGCUUACCCUGAAGUUCAUC'
        'UGUACUACAGGGAAGCUUCCCGUGCCGUGGCCCACCCUCGUCACUACCCUGACUUACGGAGUGCAGUGCUUC'
        'AGUCGCUAUCCAGACCACAUGAAGCAGCACGACUUCUUCAAGAGCGCCAUGCCCGAGGGAUACGUCCAGGAG'
        'CGUACCAUCUUCUUUAAGGACGAUGGCAAUUACAAGACCCGGGCCGAGGUGAAGUUCGAGGGCGACACCUUG'
        'GUGAACCGGAUCGAGCUGAAGGGCAUCGACUUUAAGGAGGAUGGUAACAUCCUGGGCCACAAGUUGGAGUAC'
        'AACUACAACUCCCAUAACGUCUAUAUUAUGGCGGAUAAGCAGAAGAAUGGCAUCAAGGUCAACUUCAAGAUC'
        'CGCCAUAAUAUAGAGGAUGGGAGUGUAUGAcacacugacguaugucgggccuacgcaggucauuguuacacc'
    )
    patterns = ['UUUUUU', 'GCTTAC']

    cdsend = sum(c.isupper() for c in seq)
    prot = Seq.Seq(seq[:cdsend]).translate()

    usage = codon_usage.CodonUsage('codon_usage_freq_table_human.csv')
    builder = codon_dfa.CodonDFABuilder(usage, 1)
    dfa, fakeprot, u3norm = builder.generate_dfa(
        prot, seq[cdsend:].upper().replace('T', 'U'))

    reex = RestrictionSiteExclusion(patterns)
    matches = reex.find_sites(seq)
    for match in matches:
        for aapos, codon in match['codons']:
            mask_codon_in_dfa(dfa, aapos, codon, -9999999999)
