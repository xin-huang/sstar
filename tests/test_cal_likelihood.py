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
from sstar.cal_likelihood import cal_likelihood

@pytest.fixture
def data():
    pytest.sim_null = ["./examples/data/simulated_data/sstar.example.simulated.null.src1.pval.txt", 
                       "./examples/data/simulated_data/sstar.example.simulated.null.src2.pval.txt"]
    pytest.sim_alt1 = ["./examples/data/simulated_data/sstar.example.simulated.alt1.src1.pval.txt", 
                       "./examples/data/simulated_data/sstar.example.simulated.alt1.src2.pval.txt"]
    pytest.sim_alt2 = ["./examples/data/simulated_data/sstar.example.simulated.alt2.src1.pval.txt", 
                       "./examples/data/simulated_data/sstar.example.simulated.alt2.src2.pval.txt"]
    pytest.real_data = ["./examples/results/sstar.example.match.nean.filtered.pval.txt", 
                        "./examples/results/sstar.example.match.den.filtered.pval.txt"]
    pytest.threshold = "./examples/results/sstar.example.threshold.txt"
    pytest.output = "./tests/results/test.posterior.prob.results"
    pytest.exp_output = "./tests/results/test.posterior.prob.exp.results"

def test_cal_likelihood(data):
    cal_likelihood(sim_null=pytest.sim_null, sim_alt1=pytest.sim_alt1, sim_alt2=pytest.sim_alt2, real_data=pytest.real_data, threshold=pytest.threshold, output=pytest.output, grid_size=10, bdw_null=1.5, bdw_alt1=6, bdw_alt2=6, step_size=0.005)

    with open(pytest.output, 'r') as f:
        res = [l for l in f]
    with open(pytest.exp_output, 'r') as f:
        exp_res = [l for l in f]

    assert res == exp_res
