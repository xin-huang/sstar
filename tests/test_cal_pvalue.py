# Apache License Version 2.0
# Copyright 2022 Xin Huang
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

import allel
import pytest
import numpy as np
from sstar.cal_pvalue import cal_pvalue, _cal_pvalue, _cal_pvalue_ind, _read_score_file, _read_ref_match_pct_file, _query_ref_match_pct_naive
from sstar.utils import read_data

@pytest.fixture
def data():
    pytest.ref_ind_file = "./examples/data/ind_list/ref.ind.list"
    pytest.tgt_ind_file = "./examples/data/ind_list/tgt.ind.list"
    pytest.src_ind_file = "./examples/data/ind_list/nean.ind.list"
    pytest.vcf = "./tests/data/test.pvalue.data.vcf"
    pytest.score_file = "./tests/results/test.pvalue.score.exp.results"
    pytest.output = "./tests/results/test.pvalue.results"
    pytest.exp_output = "./tests/results/test.pvalue.exp.results"
    pytest.ref_match_pct = "./tests/data/test.pvalue.ref.match.pct"

def test_cal_pvalues(data):
    cal_pvalue(pytest.vcf, pytest.ref_ind_file, pytest.tgt_ind_file, pytest.src_ind_file, None, pytest.output, 1, pytest.score_file, pytest.ref_match_pct, None, low_memory=False, mapped_len_esp=1000, len_esp=1000, var_esp=1, sfs_esp=0.02)
    with open(pytest.output, 'r') as f:
        res = [l for l in f]
    with open(pytest.exp_output, 'r') as f:
        exp_res = [l for l in f]

    assert res == exp_res

    cal_pvalue(pytest.vcf, pytest.ref_ind_file, pytest.tgt_ind_file, pytest.src_ind_file, None, pytest.output, 1, pytest.score_file, pytest.ref_match_pct, None, low_memory=True, mapped_len_esp=1000, len_esp=1000, var_esp=1, sfs_esp=0.02)
    with open(pytest.output, 'r') as f:
        res = [l for l in f]
    with open(pytest.exp_output, 'r') as f:
        exp_res = [l for l in f]

    assert res == exp_res

def test_cal_pvalue_ind(data):
    ref_data, ref_samples, tgt_data, tgt_samples, src_data, src_samples = read_data(pytest.vcf, pytest.ref_ind_file, pytest.tgt_ind_file, pytest.src_ind_file, None)
    chr_names = ref_data.keys()

    data, windows, samples = _read_score_file(pytest.score_file, chr_names, tgt_samples)
    ref_match_pct = _read_ref_match_pct_file(pytest.ref_match_pct)
    for s in samples[0:1]:
        i = samples.index(s)
        res = _cal_pvalue_ind(data[s], i, None, tgt_data, src_data, src_samples, ref_match_pct, len(samples), _query_ref_match_pct_naive, mapped_len_esp=1000, len_esp=1000, var_esp=1, sfs_esp=0.02)

    assert res == ['21\t9400000\t9450000\tNA06986\tNA\tNean\t1\tNA\tNA\t0\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t0', '21\t9400000\t9450000\tNA06986\tNA\tNean\t2\tNA\tNA\t0\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t0', '21\t9410000\t9460000\tNA06986\tNA\tNean\t1\tNA\tNA\t0\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t0', '21\t9410000\t9460000\tNA06986\tNA\tNean\t2\tNA\tNA\t0\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t0', '21\t9420000\t9470000\tNA06986\tNA\tNean\t1\tNA\tNA\t0\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t0', '21\t9420000\t9470000\tNA06986\tNA\tNean\t2\tNA\tNA\t0\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t0']

def test_cal_pvalue():
    ref_match_pct = np.array([0.045454545454545456, 0.07142857142857142, 0.037037037037037035, 0.06666666666666667, 0.06896551724137931, 0.043478260869565216, 0.1, 0.058823529411764705, 0.05263157894736842, 0.07142857142857142])

    assert 'NA' == _cal_pvalue(ref_match_pct, 'NA')
    assert 0.1 == _cal_pvalue(ref_match_pct, 0.09090909090909091)
    assert 1.0 == _cal_pvalue(ref_match_pct, 0.0)
