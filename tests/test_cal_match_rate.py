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

"""
tests/test_cal_match_rate.py

This test suite validates:
1) Unit behavior of calc_match_pct (dosage-based match rate).
2) End-to-end cal_match_pct pipeline in both modes:
   - phased=True  -> matches phased golden file (existing)
   - phased=False -> matches unphased golden file (new)
3) No shared outputs: all integration tests write to tmp_path.
"""

import pytest
import numpy as np
import allel

from sstar.cal_match_rate import calc_match_pct, cal_match_pct


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
        # Existing repo golden (this corresponds to phased/haplotype logic).
        "exp_output_phased": "./tests/results/test.match.rate.exp.results",
        # New repo golden (you must create/commit this file after generating it once).
        "exp_output_unphased": "./tests/results/test.match.rate.exp.unphased.results",
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


# ---------------------------------------
# Helpers: integration output validation
# ---------------------------------------


def _read_lines(path: str):
    with open(path, "r") as f:
        return f.read().splitlines()


def _assert_output_sane(lines):
    """
    Validate structure and numeric invariants of cal_match_pct output.

    Guarantees:
    - correct header
    - each row has 6 columns
    - match_rate is "NA" or a float within [0,1]
    """
    assert len(lines) > 1, "Expected header + at least one row"

    header = lines[0].split("\t")
    assert header == ["chrom", "start", "end", "sample", "match_rate", "src_sample"]

    for row in lines[1:]:
        fields = row.split("\t")
        assert len(fields) == 6, f"Expected 6 columns, got {len(fields)}: {row}"
        mr = fields[4]
        if mr != "NA":
            v = float(mr)
            assert 0.0 <= v <= 1.0, f"match_rate out of bounds: {v}"


# ---------------------------------------
# Integration: phased & unphased goldens
# ---------------------------------------


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

    got = out.read_text().splitlines()
    _assert_output_sane(got)

    exp = _read_lines(data["exp_output_phased"])
    assert got == exp


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

    got = out.read_text().splitlines()
    _assert_output_sane(got)

    exp = _read_lines(data["exp_output_unphased"])
    assert got == exp
