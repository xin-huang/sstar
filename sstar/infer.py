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
import pickle
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sstar.models import LogisticRegression, ExtraTrees, Sstar


def infer(feature_file, model_file, prediction_dir, prediction_prefix, cutoff, algorithm=None):
    """
    """
    with open(model_file, 'rb') as f:
        trained_model = pickle.load(f)

    feature_df = pd.read_csv(feature_file, sep="\t")

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

    print(trained_model.classes_)
    labels = model.infer(trained_model, feature_df)
    feature_df['label_0_prob'] = labels[:,0]
    feature_df['label_1_prob'] = labels[:,1]
    #feature_df['label'] = model.infer(trained_model, feature_df)
    #feature_df = feature_df[feature_df['label']>cutoff][['chrom', 'start', 'end', 'sample']].sort_values(by=['sample', 'chrom', 'start', 'end'])
    
    feature_df.to_csv(prediction_file, sep="\t", index=False, header=None)


if __name__ == '__main__':
    infer(feature_file="/scratch/admixlab/xinhuang/projects/sstar2-analysis-dev/results/test_data/ArchIE_3D19/nref_50/ntgt_50/956714/0/sim.test.0.archie.features", model_file="./sstar/test/test.archie.imbalanced.scaled.model", prediction_dir="./sstar/test", prediction_prefix="test.imbalanced.scaled", cutoff=0.5, algorithm="logistic_regression")
