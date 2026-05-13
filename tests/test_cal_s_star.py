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
import pandas as pd
from sstar.cal_s_star import cal_s_star, _cal_score_ind


@pytest.fixture
def data():
    pytest.ref_ind_list = "./tests/data/test.ref.ind.list"
    pytest.tgt_ind_list = "./tests/data/test.tgt.ind.list"
    pytest.vcf = "./tests/data/test.score.data.vcf"
    pytest.anc_allele = "./tests/data/test.anc.allele"
    pytest.exp_unphased_output = "./tests/results/test.score.unphased.exp.results.tsv"
    pytest.exp_phased_output = "./tests/results/test.score.phased.exp.results.tsv"


def test_cal_s_star_with_unphased_data(data, tmp_path):
    output = tmp_path / "test.score.unphased.results"

    cal_s_star(
        vcf=pytest.vcf,
        ref_ind_file=pytest.ref_ind_list,
        tgt_ind_file=pytest.tgt_ind_list,
        anc_allele_file=None,
        output=output,
        win_len=50000,
        win_step=10000,
        thread=1,
        match_bonus=5000,
        max_mismatch=5,
        mismatch_penalty=-10000,
        is_phased=False,
    )

    df = pd.read_csv(output, sep="\t")
    df_expected = pd.read_csv(pytest.exp_unphased_output, sep="\t")

    pd.testing.assert_frame_equal(
        df,
        df_expected,
        check_dtype=False,
        check_exact=False,
        rtol=1e-6,
        atol=1e-8,
    )


def test_cal_s_star_with_phased_data(data, tmp_path):
    """
    Run the SAME VCF as the unphased test, but treat GT as PHASED (haplotypes).

    We assert:
      - header is correct
      - output has at least one data row
      - both haplotype suffixes ('_1', '_2') are present in sample labels
      - phased output has at least as many rows as the unphased expected output
    """
    output = tmp_path / "test.score.phased.results"

    cal_s_star(
        vcf=pytest.vcf,
        ref_ind_file=pytest.ref_ind_list,
        tgt_ind_file=pytest.tgt_ind_list,
        anc_allele_file=None,
        output=output,
        win_len=50000,
        win_step=10000,
        thread=1,
        match_bonus=5000,
        max_mismatch=5,
        mismatch_penalty=-10000,
        is_phased=True,
    )

    df = pd.read_csv(output, sep="\t")
    df_expected = pd.read_csv(pytest.exp_phased_output, sep="\t")

    pd.testing.assert_frame_equal(
        df,
        df_expected,
        check_dtype=False,
        check_exact=False,
        rtol=1e-6,
        atol=1e-8,
    )


def test_cal_score_ind_with_too_few_snps():
    """
    len(tgt_gt) <= 2 -> early-exit branch.
    While-loop should never execute and we expect no output lines.
    """
    ref_pos = np.array([100, 200])
    tgt_pos = np.array([100, 200])
    tgt_gt = np.array([1, 1])  # non-ref at both sites

    res = _cal_score_ind(
        chr_name="chr1",
        sample_name="ind1",
        ref_pos=ref_pos,
        tgt_pos=tgt_pos,
        tgt_gt=tgt_gt,
        win_step=100,
        win_len=500,
        match_bonus=5000,
        max_mismatch=5,
        mismatch_penalty=-10000,
        chr_last_pos=tgt_pos[-1],
    )

    assert res == []


def test_cal_score_ind_with_perfect_ld_chain():
    """
    Three SNPs with identical genotypes in perfect LD.
    Should produce at least one window with numeric S* score and SNP count.
    """
    ref_pos = np.array([100, 300, 600])
    tgt_pos = np.array([100, 300, 600])
    tgt_gt = np.array([1, 1, 1])  # all non-ref, same "genotype"

    res = _cal_score_ind(
        chr_name="chr1",
        sample_name="ind1",
        ref_pos=ref_pos,
        tgt_pos=tgt_pos,
        tgt_gt=tgt_gt,
        win_step=100,
        win_len=1000,
        match_bonus=5000,
        max_mismatch=5,
        mismatch_penalty=-10000,
        chr_last_pos=tgt_pos[-1],
    )

    assert len(res) >= 1
    fields = res[0].split("\t")
    assert len(fields) == 8

    s_star_score = fields[5]
    s_star_snp_num = fields[7]

    assert s_star_score != "NA"
    assert s_star_snp_num != "NA"


def test_cal_score_ind_with_mismatch_penalty():
    """
    Force genotype distances of 1 to trigger the mismatch_penalty branch.
    """
    ref_pos = np.array([100, 300, 600])
    tgt_pos = np.array([100, 300, 600])
    # Differences of 1 between some SNPs: triggers mismatch_penalty
    tgt_gt = np.array([0, 1, 2])

    res = _cal_score_ind(
        chr_name="chr1",
        sample_name="ind1",
        ref_pos=ref_pos,
        tgt_pos=tgt_pos,
        tgt_gt=tgt_gt,
        win_step=100,
        win_len=1000,
        match_bonus=5000,
        max_mismatch=1,  # diff=1 allowed, diff=2 treated as > max_mismatch
        mismatch_penalty=-10000,
        chr_last_pos=tgt_pos[-1],
    )

    assert len(res) >= 1
    fields = res[0].split("\t")
    assert len(fields) == 8
    # S*_score may be "NA" or numeric depending on the optimal chain,
    # but the function must run and produce well-formed output.


def test_cal_score_ind_with_short_phy_distance():
    """
    SNPs very close in bp trigger the phy_dist < 10 branch (score = -inf).
    """
    ref_pos = np.array([100, 105, 300])
    tgt_pos = np.array([100, 105, 300])
    tgt_gt = np.array([1, 1, 1])

    res = _cal_score_ind(
        chr_name="chr1",
        sample_name="ind1",
        ref_pos=ref_pos,
        tgt_pos=tgt_pos,
        tgt_gt=tgt_gt,
        win_step=50,
        win_len=500,
        match_bonus=5000,
        max_mismatch=5,
        mismatch_penalty=-10000,
        chr_last_pos=tgt_pos[-1],
    )

    assert len(res) >= 1
    fields = res[0].split("\t")
    assert len(fields) == 8
