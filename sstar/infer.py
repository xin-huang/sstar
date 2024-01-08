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

pd.options.mode.chained_assignment = None


def infer(feature_file, model_file, prediction_dir, prediction_prefix, algorithm=None):
    """
    """
    with open(model_file, 'rb') as f:
        trained_model = pickle.load(f)

    feature_df = pd.read_csv(feature_file, sep="\t")

    if algorithm == 'logistic_regression':
        prediction_file = prediction_dir + '/' + prediction_prefix + '.logistic.regression.predictions'
        model = LogisticRegression()
    elif algorithm == 'extra_trees':
        prediction_file = prediction_dir + '/' + prediction_prefix + '.extra.trees.predictions'
        model = ExtraTrees()
    elif (algorithm == 'sstar') or (algorithm is None):
        prediction_file = prediction_dir + '/' + prediction_prefix + '.sstar.predictions'
        model = Sstar()
    else:
        raise Exception(f'The {algorithm} algorithm is NOT available!')

    classes = trained_model.classes_
    predictions = model.infer(trained_model, feature_df)
    prediction_df = feature_df[['chrom', 'start', 'end', 'sample']]

    # Rename class 0 as non-introgressed?
    # Rename class 1 as introgressed?

    for i in range(len(classes)):
        prediction_df[f'class_{classes[i]}_prob'] = predictions[:,i]
    
    prediction_df.sort_values(by=['sample', 'chrom', 'start', 'end']).to_csv(prediction_file, sep="\t", index=False)


if __name__ == '__main__':
    infer(feature_file="/scratch/admixlab/xinhuang/projects/sstar2-analysis-dev/results/test_data/ArchIE_3D19/nref_50/ntgt_50/956714/0/sim.test.0.archie.features", model_file="/scratch/admixlab/xinhuang/projects/sstar2-analysis-dev/tmp/archie.imbalanced.logistic_regression.model", prediction_dir="./sstar/test", prediction_prefix="test.imbalanced.test", algorithm="logistic_regression")
