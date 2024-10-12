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
from sstar.calc_match_rate import calc_match_pct
from sstar.utils import read_data


@pytest.fixture
def data():
    pytest.ref_ind_file = "./examples/data/ind_list/ref.ind.list"
    pytest.tgt_ind_file = "./examples/data/ind_list/tgt.ind.list"
    pytest.src_ind_file = "./examples/data/ind_list/nean.ind.list"
    pytest.vcf = "./tests/data/test.match.rate.data.vcf"
    pytest.score_file = "./tests/results/test.match.rate.score.exp.results"
    pytest.output = "./tests/results/test.match.rate.results"
    pytest.exp_output = "./tests/results/test.match.rate.exp.results"


def test_calc_match_pct(data):
    calc_match_pct(
        pytest.vcf,
        pytest.ref_ind_file,
        pytest.tgt_ind_file,
        pytest.src_ind_file,
        None,
        pytest.output,
        1,
        pytest.score_file,
        None,
    )
    with open(pytest.output, "r") as f:
        res = [l for l in f]
    with open(pytest.exp_output, "r") as f:
        exp_res = [l for l in f]

    assert res == exp_res
