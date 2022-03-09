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
from sstar.cal_threshold import cal_threshold

@pytest.fixture
def data():
    pytest.simulated_data = "./examples/data/simulated_data/gravel_asn_scale_60k.simulated.data"
    pytest.score_file = "./tests/results/test.score.exp.results"
    pytest.recomb_map = "./examples/data/real_data/hum.windows.50k.10k.recomb.map"
    pytest.output = "./tests/results/test.threshold.results"
    pytest.exp_output = "./tests/results/test.threshold.exp.results"

def test_cal_threshold(data):
    cal_threshold(simulated_data=pytest.simulated_data, score_file=pytest.score_file, recomb_rate=0, recomb_map=pytest.recomb_map, quantile=0.99, output=pytest.output, k=8)

    f1 = open(pytest.output, 'r')
    res = f1.read()
    f1.close()
    f2 = open(pytest.exp_output, 'r')
    exp_res = f2.read()
    f2.close()

    assert res == exp_res
