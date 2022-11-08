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
import numpy as np
from sstar.preprocess import process_data


@pytest.fixture
def data():
    pytest.ref_ind_list = "./tests/data/test.ref.ind.list"
    pytest.tgt_ind_list = "./tests/data/test.tgt.ind.list"
    pytest.vcf = "./tests/data/test.score.data.vcf"
    pytest.sstar_output = "./tests/results/test.sstar.out"
    pytest.archie_output = "./tests/results/test.archie.out"
    pytest.sstar_exp_output = "./tests/results/test.exp.sstar.out"
    pytest.archie_exp_output = "./tests/results/test.exp.archie.out"


def test_process_data(data):
    process_data(pytest.vcf, pytest.ref_ind_list, pytest.tgt_ind_list, None, pytest.sstar_output, 50000, 10000, 2, 5000, 5, -10000, process_archie=False)
    process_data(pytest.vcf, pytest.ref_ind_list, pytest.tgt_ind_list, None, pytest.archie_output, 50000, 10000, 2, 5000, 5, -10000, process_archie=True)

    res = [line for line in open(pytest.sstar_output, 'r').readlines()]
    exp_res = [line for line in open(pytest.sstar_exp_output, 'r').readlines()]

    assert res == exp_res

    res = [line for line in open(pytest.archie_output, 'r').readlines()]
    exp_res = [line for line in open(pytest.archie_exp_output, 'r').readlines()]

    assert res == exp_res
