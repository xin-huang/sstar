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
import numpy as np
from sstar.stats import cal_pr


def evaluate(truth_tract_file, inferred_tract_file, output):
    """
    """
    truth_tracts = pd.read_csv(truth_tract_file, sep="\t", header=None)
    inferred_tracts = pd.read_csv(inferred_tract_file, sep="\t", header=None)

    truth_tracts.columns = ['chrom', 'start', 'end', 'sample']
    inferred_tracts.columns = ['chrom', 'start', 'end', 'sample']

    truth_tracts_samples = truth_tracts['sample'].unique()
    inferred_tracts_samples = inferred_tracts['sample'].unique()

    print(truth_tracts_samples)
    print(inferred_tracts_samples)

    res = pd.DataFrame(columns=['sample', 'precision', 'recall'])

    for s in np.intersect1d(truth_tracts_samples, inferred_tracts_samples):
        ind_truth_tracts = truth_tracts[truth_tracts['sample'] == s][['chrom', 'start', 'end']]
        ind_inferred_tracts = inferred_tracts[inferred_tracts['sample'] == s][['chrom', 'start', 'end']]

        ind_truth_tracts = pybedtools.BedTool.from_dataframe(ind_truth_tracts).sort().merge()
        ind_inferred_tracts = pybedtools.BedTool.from_dataframe(ind_inferred_tracts).sort().merge()
        if s == 'tsk_50_1': 
            print(s)
            print(ind_truth_tracts)
            print(ind_inferred_tracts)

        ntruth_tracts = sum([x.stop - x.start for x in (ind_truth_tracts)])
        ninferred_tracts = sum([x.stop - x.start for x in (ind_inferred_tracts)])
        ntrue_positives = sum([x.stop - x.start for x in ind_inferred_tracts.intersect(ind_truth_tracts)])

        precision, recall = cal_pr(ntruth_tracts, ninferred_tracts, ntrue_positives)

        res.loc[len(res.index)] = [s, precision, recall]

    for s in np.setdiff1d(truth_tracts_samples, inferred_tracts_samples):
        # ninferred_tracts = 0
        res.loc[len(res.index)] = [s, np.nan, 0]

    for s in np.setdiff1d(inferred_tracts_samples, truth_tracts_samples):
        # ntruth_tracts = 0
        res.loc[len(res.index)] = [s, 0, np.nan]

    res.sort_values(by=['sample']).to_csv(output, sep="\t", index=False)


if __name__ == '__main__':
    evaluate(truth_tract_file="sim.test.0.truth.tracts.bed", inferred_tract_file="tsk.cutoff.0.5.predicted.bed", output="./sstar/test/test.performance")
