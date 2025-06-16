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

from importlib.resources import files, as_file

with as_file(files('pldesign')) as pldesign_prefix:
    PLDESIGN_DIR = f'{pldesign_prefix}'
    CODON_USAGE_FILE = f'{pldesign_prefix}/codon_usage_freq_table_human.csv'
    PLDESIGN_BINARY = f'{pldesign_prefix}/paralinear_design'

del as_file, files, pldesign_prefix

def run_design(protein, utr3='', ldlambda=0.1, masks=[], caigamma=1,
               codon_usage_path=CODON_USAGE_FILE, exclude_patterns=[],
               minimum_cai_weight=0, minimum_rscu=0, allow_no_stop=False,
               skip_vienna_reevaluation=False, random_rejection_rate=0.0,
               random_seed=0, validate=True, quiet=False):

    from .pldesign import ParalinearDesign

    with ParalinearDesign(
            protein=protein, utr3=utr3, ldlambda=ldlambda, caigamma=caigamma,
            codonusage_filename=codon_usage_path, exclude_patterns=exclude_patterns,
            allow_no_stop=allow_no_stop, minimum_rscu=minimum_rscu,
            minimum_cai_weight=minimum_cai_weight,
            skip_vienna_reevaluation=skip_vienna_reevaluation,
            quiet=quiet) as design:
        return design.run(masks, random_rejection_rate=random_rejection_rate,
                          random_seed=random_seed, validate=validate)