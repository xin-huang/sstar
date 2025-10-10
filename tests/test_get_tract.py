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
from sstar.get_tract import get_tract


@pytest.fixture
def data():
    pytest.threshold = "./tests/data/test.tract.threshold"
    pytest.src1_match_pct = "./tests/data/test.tract.src1.match.rate"
    pytest.src2_match_pct = "./tests/data/test.tract.src2.match.rate"
    pytest.exp_bed = "./tests/results/test.tract.exp.bed"
    pytest.exp_bed_with_src = "./tests/results/test.tract.with.src.match.rate.exp.bed"
    pytest.exp_src1_bed = "./tests/results/test.tract.exp.src1.bed"
    pytest.exp_src2_bed = "./tests/results/test.tract.exp.src2.bed"


def test_get_tract(data):
    get_tract(
        threshold_file=pytest.threshold,
        match_pct_files=None,
        output_prefix="./tests/results/test.tract",
        diff=0,
    )
    f1 = open("./tests/results/test.tract.bed", "r")
    res = f1.read()
    f1.close()
    f2 = open(pytest.exp_bed, "r")
    exp_res = f2.read()
    f2.close()

    assert res == exp_res

    get_tract(
        threshold_file=pytest.threshold,
        match_pct_files=[pytest.src1_match_pct],
        output_prefix="./tests/results/test.tract.with.src.match.rate",
        diff=0,
    )
    f1 = open("./tests/results/test.tract.with.src.match.rate.bed", "r")
    res = f1.read()
    f1.close()
    f2 = open(pytest.exp_bed_with_src, "r")
    exp_res = f2.read()
    f2.close()

    assert res == exp_res

    get_tract(
        threshold_file=pytest.threshold,
        match_pct_files=[pytest.src1_match_pct, pytest.src2_match_pct],
        output_prefix="./tests/results/test.tract",
        diff=0,
    )
    f1 = open("./tests/results/test.tract.src1.bed", "r")
    res1 = f1.read()
    f1.close()
    f2 = open("./tests/results/test.tract.src2.bed", "r")
    res2 = f2.read()
    f2.close()
    f3 = open(pytest.exp_src1_bed, "r")
    exp_res1 = f3.read()
    f3.close()
    f4 = open(pytest.exp_src2_bed, "r")
    exp_res2 = f4.read()
    f4.close()

    assert res1 == exp_res1
    assert res2 == exp_res2
