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
from sstar.get_quantile import get_quantile

@pytest.fixture
def data():
    pytest.model = "./examples/models/BonoboGhost_4K19_no_introgression.yaml"
    pytest.exp_quantile = "./tests/results/test.quantile.exp.summary"

def test_get_quantile(data):
    get_quantile(model=pytest.model, ms_dir='./ext/msdir', seeds=[1,2,3], N0=1000, nsamp=22, nreps=20000, ref_index=4, ref_size=20, tgt_index=3, tgt_size=2, mut_rate=1.2e-8, rec_rate=0.7e-8, seq_len=40000, snp_num_range=[25,30,5], output_dir='./tests/results/simulation', thread=2, all_ind_geno_dist=False)
    f1 = open('./tests/results/simulation/quantile.summary.txt', 'r')
    res = f1.read()
    f1.close()

    f2 = open(pytest.exp_quantile, 'r')
    exp_res = f2.read()
    f2.close()

    assert res == exp_res
