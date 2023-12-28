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


import pickle
import pandas as pd
from sstar.models import LogisticRegression, ExtraTrees, Sstar


def train(training_data, model_file, algorithm=None):
    """
    """
    feature_df = pd.read_csv(training_data, sep="\t")

    if algorithm == 'logistic_regression':
        model = LogisticRegression()
    elif algorithm == 'extra_trees':
        model = ExtraTrees()
    elif (algorithm == 'sstar') or (algorithm is None):
        model = Sstar()
    else:
        raise Exception(f'The {algorithm} algorithm is NOT available!')

    trained_model = model.train(feature_df)

    pickle.dump(trained_model, open(model_file, "wb"))


if __name__ == '__main__':
    train(training_data="/scratch/admixlab/xinhuang/projects/sstar2-analysis-dev/tmp/sim.training.archie.imbalanced.all.labeled.features", 
          model_file="./sstar/test/test.archie.imbalanced.model", algorithm='logistic_regression')
