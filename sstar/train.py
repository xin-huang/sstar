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
from multiprocessing import Process, Queue
from sstar.preprocess import preprocess
from sstar.models import LogisticRegression, ExtraTrees, Sstar
from sstar.utils import multiprocessing_manager


def train(training_data, model_file, algorithm=None):
    """
    """
    feature_df = pd.DataFrame()
    with open(training_data, 'r') as f:
        for feature_file in f:
            df = pd.read_csv(feature_file.rstrip(), sep="\t")
            feature_df = pd.concat([feature_df, df])

    labels = feature_df['label']
    data = feature_df.drop(columns=['chrom', 'start', 'end', 'sample', 'label']).values

    if algorithm == 'logistic_regression':
        model = LogisticRegression()
    elif algorithm == 'extra_trees':
        model = ExtraTrees()
    elif (algorithm == 'sstar') or (algorithm is None):
        model = Sstar()
    else:
        raise Exception(f'The {algorithm} algorithm is NOT available!')

    trained_model = model.train(data, labels)

    pickle.dump(trained_model, open(model_file, "wb"))


if __name__ == '__main__':
    train(nrep=1000, thread=2, training_data_prefix='test', training_data_dir='./sstar/test', 
          seq_len=50000, feature_config='examples/pre-trained/test.features.yaml', model_file="./sstar/test/test.lr.model",
          is_phased=True, ploidy=2, 
          archaic_prop=0.7, not_archaic_prop=0.3, algorithm='logistic_regression')
