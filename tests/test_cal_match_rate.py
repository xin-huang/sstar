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
import pandas as pd
import allel

from sstar.cal_match_rate import (
    calc_match_pct,
    cal_match_pct,
    _read_score_file,
    _cal_match_pct_ind,
)


@pytest.fixture
def data():
    """
    Centralized test inputs.

    IMPORTANT:
    - Do not store any output path here.
    - Each test uses tmp_path so outputs are never shared.
    """
    return {
        "ref_ind_file": "./examples/data/ind_list/ref.ind.list",
        "tgt_ind_file": "./examples/data/ind_list/tgt.ind.list",
        "src_ind_file": "./examples/data/ind_list/nean.ind.list",
        "vcf": "./tests/data/test.match.rate.data.vcf",
        "score_file": "./tests/results/test.match.rate.score.exp.results",
        "exp_output_phased": "./tests/results/test.match.rate.phased.exp.results",
        "exp_output_unphased": "./tests/results/test.match.rate.unphased.exp.results",
    }


# ----------------------------
# Unit tests: calc_match_pct
# ----------------------------


def test_calc_match_pct_vector_mean():
    x = np.array([1, 0, 1, 0])
    y = np.array([1, 0, 2, 1])
    # per-site: [1.0, 1.0, 0.5, 0.5] -> mean 0.75
    assert calc_match_pct(x, y, P=2) == pytest.approx(0.75)


def test_calc_match_pct_missing_ignored():
    x = np.array([0, -1, 2])
    y = np.array([0, 1, -1])
    # callable sites: only index0 -> match=1.0
    assert calc_match_pct(x, y, P=2) == pytest.approx(1.0)


def test_calc_match_pct_all_missing_returns_NA():
    x = np.array([-1, -1])
    y = np.array([-1, -1])
    assert calc_match_pct(x, y, P=2) == "NA"


def test_calc_match_pct_with_scikit_allel_genotypes():
    """
    Sanity-check: dosage extracted via scikit-allel matches the vector example.
    """
    gt = np.array(
        [
            [[0, 1], [0, 1]],  # dos: [1, 1]
            [[0, 0], [0, 0]],  # dos: [0, 0]
            [[0, 1], [1, 1]],  # dos: [1, 2]
            [[0, 0], [0, 1]],  # dos: [0, 1]
        ],
        dtype=np.int8,
    )
    g = allel.GenotypeArray(gt)
    dos = g.to_n_alt()  # (n_variants, n_samples)

    x = dos[:, 0]  # [1,0,1,0]
    y = dos[:, 1]  # [1,0,2,1]

    assert calc_match_pct(x, y, P=2) == pytest.approx(0.75)


def test_cal_match_pct_phased_matches_golden(data, tmp_path):
    """
    phased=True must match the existing golden file exactly.
    """
    out = tmp_path / "match_rate.phased.tsv"

    cal_match_pct(
        data["vcf"],
        data["ref_ind_file"],
        data["tgt_ind_file"],
        data["src_ind_file"],
        None,
        str(out),
        1,
        data["score_file"],
        None,
        phased=True,
    )

    df = pd.read_csv(out, sep="\t")
    df_expected = pd.read_csv(data["exp_output_phased"], sep="\t")

    pd.testing.assert_frame_equal(
        df,
        df_expected,
        check_dtype=False,
        check_exact=False,
        rtol=1e-6,
        atol=1e-8,
    )


def test_cal_match_pct_unphased_matches_golden(data, tmp_path):
    """
    phased=False must match the unphased golden file exactly.

    NOTE:
    - You must generate and commit ./tests/results/test.match.rate.exp.unphased.results
      once you consider the unphased logic stable/approved.
    """
    out = tmp_path / "match_rate.unphased.tsv"

    cal_match_pct(
        data["vcf"],
        data["ref_ind_file"],
        data["tgt_ind_file"],
        data["src_ind_file"],
        None,
        str(out),
        1,
        data["score_file"],
        None,
        phased=False,
    )

    df = pd.read_csv(out, sep="\t")
    df_expected = pd.read_csv(data["exp_output_unphased"], sep="\t")

    pd.testing.assert_frame_equal(
        df,
        df_expected,
        check_dtype=False,
        check_exact=False,
        rtol=1e-6,
        atol=1e-8,
    )
