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

import pytest
from sstar.cal_ref_match_pct import cal_ref_match_pct

@pytest.fixture
def data():
    pytest.ref_ind_list = "./examples/data/ind_list/ref.ind.list"
    pytest.src_ind_list = "./examples/data/ind_list/nean.ind.list"
    pytest.vcf = "./tests/data/test.pvalue.data.vcf"
    pytest.output = "./tests/results/test.ref.match.pct.results"
    pytest.exp_output = "./tests/results/test.ref.match.pct.exp.results"

def test_cal_ref_match_pct(data):
    cal_ref_match_pct(pytest.vcf, pytest.ref_ind_list, pytest.src_ind_list, None, pytest.output, 50000, 10000, 1, None)
    with open(pytest.output, 'r') as f:
        res = [l for l in f]
    with open(pytest.exp_output, 'r') as f:
        exp_res = [l for l in f]

    assert res == exp_res
