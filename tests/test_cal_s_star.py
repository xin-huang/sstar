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
from sstar.cal_s_star import cal_s_star, _cal_score_ind


@pytest.fixture
def data():
    pytest.ref_ind_list = "./tests/data/test.ref.ind.list"
    pytest.tgt_ind_list = "./tests/data/test.tgt.ind.list"
    pytest.vcf = "./tests/data/test.score.data.vcf"
    pytest.anc_allele = "./tests/data/test.anc.allele"
    pytest.output = "./tests/results/test.score.results"
    pytest.exp_output = "./tests/results/test.score.exp.results"


def test_cal_s_star_with_unphased_data(data):
    cal_s_star(
        vcf=pytest.vcf,
        ref_ind_file=pytest.ref_ind_list,
        tgt_ind_file=pytest.tgt_ind_list,
        anc_allele_file=None,
        output=pytest.output,
        win_len=50000,
        win_step=10000,
        thread=1,
        match_bonus=5000,
        max_mismatch=5,
        mismatch_penalty=-10000,
        is_phased=False,
    )
    f1 = open(pytest.output, "r")
    res = f1.read()
    f1.close()
    f2 = open(pytest.exp_output, "r")
    exp_res = f2.read()
    f2.close()

    assert res == exp_res


def test_cal_s_star_with_phased_data_same_input(data):
    """
    Run the SAME VCF as the unphased test, but treat GT as PHASED (haplotypes).

    We assert:
      - header is correct
      - output has at least one data row
      - both haplotype suffixes ('_hap1', '_hap2') are present in sample labels
      - phased output has at least as many rows as the unphased expected output
    """
    output_phased = "./tests/results/test.score.phased_same_input.results"

    cal_s_star(
        vcf=pytest.vcf,
        ref_ind_file=pytest.ref_ind_list,
        tgt_ind_file=pytest.tgt_ind_list,
        anc_allele_file=None,
        output=output_phased,
        win_len=50000,
        win_step=10000,
        thread=1,
        match_bonus=5000,
        max_mismatch=5,
        mismatch_penalty=-10000,
        is_phased=True,
    )

    # Read phased result
    with open(output_phased) as f:
        phased_lines = [l.rstrip("\n") for l in f]

    # Header + at least one data row
    assert len(phased_lines) > 1

    header = phased_lines[0].split("\t")
    assert header == [
        "chrom",
        "start",
        "end",
        "sample",
        "S*_score",
        "region_ind_SNP_number",
        "S*_SNP_number",
        "S*_SNPs",
    ]

    # Basic row integrity
    first_data = phased_lines[1].split("\t")
    assert len(first_data) == 8

    # Check that both haplotype suffixes occur across sample labels
    samples = {line.split("\t")[3] for line in phased_lines[1:]}
    assert any(s.endswith("_hap1") for s in samples)
    assert any(s.endswith("_hap2") for s in samples)

    # Compare number of rows to *expected* unphased output (static file)
    with open(pytest.exp_output) as f:
        unphased_lines = [l.rstrip("\n") for l in f]

    # Ignore headers (first line)
    phased_n = len(phased_lines) - 1
    unphased_n = len(unphased_lines) - 1

    # Phased run has at least as many windows as the unphased one
    assert phased_n >= unphased_n

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
    )

    assert len(res) >= 1
    fields = res[0].split("\t")
    assert len(fields) == 8

    s_star_score = fields[4]
    s_star_snp_num = fields[6]

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
    )

    assert len(res) >= 1
    fields = res[0].split("\t")
    assert len(fields) == 8


