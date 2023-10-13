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
import pandas as pd
from sstar.models import LogisticRegression, ExtraTrees, Sstar


def infer(feature_file, model_file, prediction_dir, prediction_prefix, algorithm=None):
    """
    """
    feature_df = pd.read_csv(feature_file, sep="\t")
    feature_df = feature_df.drop(columns=['chrom', 'start', 'end', 'sample', 'hap'])

    if algorithm == 'logistic_regression':
        prediction_file = prediction_dir + '/' + prediction_prefix + '.logistic.regression.predicted.bed'
        model = LogisticRegression()
    elif algorithm == 'extra_trees':
        prediction_file = prediction_dir + '/' + prediction_prefix + '.extra.trees.predicted.bed'
        model = ExtraTrees()
    elif (algorithm == 'sstar') or (algorithm is None):
        prediction_file = prediction_dir + '/' + prediction_prefix + '.sstar.predicted.bed'
        model = Sstar()
    else:
        raise Exception(f'The {algorithm} algorithm is NOT available!')

    model.infer(model_file, feature_df, prediction_file)


if __name__ == '__main__':
    infer(vcf_file="./sstar/test/0/test.0.vcf", 
          ref_ind_file="./sstar/test/0/test.0.ref.ind.list", tgt_ind_file="sstar/test/0/test.0.tgt.ind.list", 
          anc_allele_file=None, win_len=50000, win_step=50000, thread=8, 
          match_bonus=5000, max_mismatch=5, mismatch_penalty=-10000, model_file="./sstar/test/test.logistic.regression.model", 
          prediction_dir="./sstar/test", prediction_prefix="test", algorithm="logistic_regression")
