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


import statsmodels.api as sm
import statsmodels.formula.api as smf
import numpy as np


def train_logistic_regression(train_df, model_file):
    """
    Description:
        Function for training of the statsmodels logistic classification.

    Arguments:
        train_df pandas.DataFrame: Training data
        save_filename str: filename for output model
    """
    sm_data_exog = train_df.copy()
    sm_data_exog.drop(['label'], axis=1, inplace=True)
    sm_data_exog = sm.add_constant(sm_data_exog, prepend=False)

    sm_data_endog = train_df['label']

    glm_binom = sm.GLM(sm_data_endog.astype(int), sm_data_exog.astype(float),family=sm.families.Binomial())
    result = glm_binom.fit()

    result.save(model_file)


def infer_logistic_regression():
    pass


def train_sstar():
    pass


def train_extra_trees():
    pass
