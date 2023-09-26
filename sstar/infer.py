# Apache License Version 2.0
# Copyright 2023 Xin Huang
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


import os
from sstar.preprocess import process_data


def infer(vcf_file, ref_ind_file, tgt_ind_file, anc_allele_file, win_len, win_step, thread, match_bonus, max_mismatch, mismatch_penalty, model_file, output_dir, output_prefix, algorithm=None):
    """
    """
    os.makedirs(output_dir, exist_ok=True)
    feature_file = output_dir + '/' + output_prefix + '.features'

    process_data(vcf_file=vcf_file, ref_ind_file=ref_ind_file, tgt_ind_file=tgt_ind_file,
                 anc_allele_file=anc_allele_file, output=feature_file, thread=thread,
                 win_len=win_len, win_step=win_step, match_bonus=match_bonus, max_mismatch=max_mismatch, 
                 mismatch_penalty=mismatch_penalty)

    if algorithm == 'logistic_regression':
        _infer_logistic_regression()
    elif algorithm == 'extra_trees':
        _infer_extra_trees()
    elif (algorithm == 'sstar') or (algorithm is None):
        _infer_sstar()
    else:
        raise Exception(f'The {algorithm} algorithm is NOT available!')


def _infer_logistic_regression():
    pass


def _infer_extra_trees():
    pass


def _infer_sstar():
    pass


if __name__ == '__main__':
    infer(vcf_file="./examples/data/real_data/sstar.example.biallelic.snps.vcf.gz", ref_ind_file="./examples/data/ind_list/ref.ind.list", tgt_ind_file="./examples/data/ind_list/tgt.ind.list", 
          anc_allele_file=None, win_len=50000, win_step=50000, thread=8, match_bonus=5000, max_mismatch=5, mismatch_penalty=-10000, model_file="./examples/pre-trained/test.logistic.regression.model", 
          output_dir="./sstar/test", output_prefix="test", algorithm="logistic_regression")
