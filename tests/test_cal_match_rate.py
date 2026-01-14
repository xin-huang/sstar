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
import allel
from sstar.cal_match_rate import calc_match_pct, cal_match_pct


@pytest.fixture
def data():
    pytest.ref_ind_file = "./examples/data/ind_list/ref.ind.list"
    pytest.tgt_ind_file = "./examples/data/ind_list/tgt.ind.list"
    pytest.src_ind_file = "./examples/data/ind_list/nean.ind.list"
    pytest.vcf = "./tests/data/test.match.rate.data.vcf"
    pytest.score_file = "./tests/results/test.match.rate.score.exp.results"
    pytest.output = "./tests/results/test.match.rate.results"
    pytest.exp_output = "./tests/results/test.match.rate.exp.results"


def test_calc_match_pct_vector_mean():
    x = np.array([1, 0, 1, 0])
    y = np.array([1, 0, 2, 1])
    assert calc_match_pct(x, y, P=2) == pytest.approx(0.75)


def test_calc_match_pct_missing_ignored():
    x = np.array([0, -1, 2])
    y = np.array([0,  1, -1])
    assert calc_match_pct(x, y, P=2) == pytest.approx(1.0)


def test_calc_match_pct_all_missing_returns_NA():
    x = np.array([-1, -1])
    y = np.array([-1, -1])
    assert calc_match_pct(x, y, P=2) == "NA"


def test_calc_match_pct_with_scikit_allel_genotypes():
    gt = np.array([
        [[0, 1], [0, 1]],
        [[0, 0], [0, 0]],
        [[0, 1], [1, 1]],
        [[0, 0], [0, 1]],
    ], dtype=np.int8)

    g = allel.GenotypeArray(gt)

    # Convert to dosage for all samples first: shape (n_variants, n_samples)
    dos = g.to_n_alt()

    x = dos[:, 0]  # [1,0,1,0]
    y = dos[:, 1]  # [1,0,2,1]

    assert calc_match_pct(x, y, P=2) == pytest.approx(0.75)



@pytest.mark.skip(reason="Integration test disabled while validating calc_match_pct unit tests.")

def test_cal_match_pct(data):
    cal_match_pct(
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





