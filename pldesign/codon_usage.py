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

import pandas as pd
import numpy as np
from functools import lru_cache


class CodonUsage:

    MINIMUM_FREQ = 1e-6

    def __init__(self, codonusagefile):
        self.codonusagefile = codonusagefile

        self.load_codon_usage_table()
        self.calc_cai_weights()

    def load_codon_usage_table(self):
        self.codon_table = pd.read_csv(self.codonusagefile, names=['codon', 'aa', 'freq'],
                                       header=None, skiprows=1)
        self.codon_to_aa = dict(zip(self.codon_table['codon'], self.codon_table['aa']))
        self.aa_to_codons = self.codon_table.groupby('aa')['codon'].apply(list).to_dict()

        assert len(self.aa_to_codons) == 21, 'Invalid number of amino acids'

    def calc_cai_weights(self):
        tbl = self.codon_table.set_index('codon')
        self.log_cai_w = np.log(
            tbl['freq'].clip(self.MINIMUM_FREQ) /
            tbl.groupby('aa')['freq'].transform('max')).to_dict()

    def calc_cai(self, seq):
        assert len(seq) % 3 == 0

        _ = self.log_cai_w
        sumlogw = sum(_[seq[i:i+3]] for i in range(0, len(seq), 3))
        return np.exp(sumlogw / (len(seq) // 3))

@lru_cache()
def get_codon_usage(codonusagefile):
    return CodonUsage(codonusagefile)

def reverse_translate_most_common(seq, codonusagefile):
    cu = get_codon_usage(codonusagefile)
    bestcodons = cu.codon_table.sort_values('freq').groupby('aa')['freq'].idxmax()
    aa_to_bestcodons = cu.codon_table.loc[bestcodons].set_index('aa')['codon'].to_dict()
    return ''.join(aa_to_bestcodons[aa] for aa in seq)