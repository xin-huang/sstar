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


import pandas as pd


def label(feature_file, truth_tract_file, seq_len, archaic_prop, not_archaic_prop, output):
    """
    """
    feature_df = pd.read_csv(feature_file, sep="\t")

    try:
        truth_tract_df = pd.read_csv(truth_tract_file, sep="\t", header=None)
    except pd.errors.EmptyDataError:
        feature_df['label'] = 0.0
    else:
        truth_tract_df.columns = ['chrom', 'start', 'end', 'sample']
        truth_tract_df['len'] = truth_tract_df['end'] - truth_tract_df['start']
        truth_tract_df = truth_tract_df.groupby(by=['sample'])['len'].sum().reset_index()
        truth_tract_df['prop'] = truth_tract_df['len'] / seq_len
        truth_tract_df['label'] = truth_tract_df.apply(lambda row: _add_label(row, archaic_prop, not_archaic_prop), axis=1)
        feature_df = feature_df.merge(truth_tract_df.drop(columns=['len', 'prop']),
                                      left_on=['sample'], right_on=['sample'], how='left').fillna(0)
    finally:
        feature_df.to_csv(output, sep="\t", index=False)


def _add_label(row, archaic_prop, not_archaic_prop):
    """
    """
    if row['prop'] > archaic_prop: return 1.0
    elif row['prop'] < not_archaic_prop: return 0.0
    else: return -1.0


if __name__ == '__main__':
    label(feature_file="", truth_tract_file="", seq_len=50000, archaic_prop=0.7, not_archaic_prop=0.3, output="")
