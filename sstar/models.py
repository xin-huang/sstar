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
# limitations under the License


import pickle
import statsmodels.api as sm
import statsmodels.formula.api as smf
import numpy as np
import pandas as pd
from abc import ABC, abstractmethod
from sklearn.linear_model import LogisticRegression as LR
from sklearn.ensemble import ExtraTreesClassifier

class Model(ABC):
    """
    """
    @abstractmethod
    def train(self, train_df, model_file):
        pass


    @abstractmethod
    def infer(self, model_file, test_df, prediction_file):
        pass


class LogisticRegression(Model):
    """
    """
    def train(self, train_df, model_file):
        """
        Description:
            Function for training of the sklearn logistic classification.

        Arguments:
            train_df pandas.DataFrame: Training data
            save_filename str: filename for output model
        """
        data = train_df.copy()
        data = data.drop(columns=['label'])
        labels = train_df['label']

        model = LR(solver="newton-cg", penalty=None, max_iter=1000)
        model.fit(data, labels.astype(int))

        pickle.dump(model, open(model_file, "wb"))


    def infer(self, model_file, test_df, prediction_file):
        """
        """
        with open(model_file, 'rb') as f:
            model = pickle.load(f)

        data = test_df.copy()
        labels = model.predict(data)
        test_df['label'] = labels

        test_df.to_csv(prediction_file, sep="\t", index=False)


class ExtraTrees(Model):
    """
    """
    
    def train(self, train_df, model_file):
        """
        Description:
            Function for training of the sklearn extraTrees classifier.

        Arguments:
            train_df pandas.DataFrame: Training data
            model_file str: filename for output model
        """
        data = train_df.copy()
        data = data.drop(columns=['label'])
        labels = train_df['label']

        model = ExtraTreesClassifier()
        model.fit(data, labels.astype(int))

        pickle.dump(model, open(model_file, "wb"))


    def infer(self, model_file, test_df, prediction_file):
        """
        Description:
            Function for predicting using the sklearn extraTrees classifier.

        Arguments:
            model_file str: filename for input model
            test_df pandas.DataFrame: test data
            prediction_file str: filename for output file containing predictions   
        """
        with open(model_file, 'rb') as f:
            model = pickle.load(f)

        data = test_df.copy()
        labels = model.predict(data)
        test_df['label'] = labels

        test_df.to_csv(prediction_file, sep="\t", index=False)


class Sstar(Model):
    """
    """
    pass
