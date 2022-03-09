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
from sstar.cal_s_star import cal_s_star

@pytest.fixture
def data():
    pytest.ref_ind_list = "./tests/data/test.ref.ind.list"
    pytest.tgt_ind_list = "./tests/data/test.tgt.ind.list"
    pytest.vcf = "./tests/data/test.data.vcf"
    pytest.anc_allele = "./tests/data/test.anc.allele"
    pytest.output = "./tests/results/test.score.results"
    pytest.exp_output = "./tests/results/test.score.exp.results"

def test_cal_s_star(data):
    cal_s_star(vcf=pytest.vcf, ref_ind_file=pytest.ref_ind_list, tgt_ind_file=pytest.tgt_ind_list, anc_allele_file=None, output=pytest.output, win_len=50000, win_step=10000, thread=1, match_bonus=5000, max_mismatch=5, mismatch_penalty=-10000)
    f1 = open(pytest.output, 'r')
    res = f1.read()
    f1.close()
    f2 = open(pytest.exp_output, 'r')
    exp_res = f2.read()
    f2.close()

    assert res == exp_res
