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


from sstar.preprocess import process_data


def infer(vcf_file, ref_ind_file, tgt_ind_file, anc_allele_file, win_len, win_step, thread, match_bonus, max_mismatch, mismatch_penalty, model_file, output_file, algorithm=None):
    """
    """
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
    infer()
