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


import pybedtools
import pandas as pd
from sstar.stats import cal_pr


def evaluate(true_tract_file, inferred_tract_file, output):
    """
    """
    true_tracts = pd.read_csv(true_tract_file, sep="\t", header=None)
    inferred_tracts = pd.read_csv(inferred_tract_file, sep="\t", header=None)

    true_tracts.columns = ['chrom', 'start', 'end', 'sample']
    inferred_tracts.columns = ['chrom', 'start', 'end', 'sample']
