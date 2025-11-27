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

import sstar.cal_match_rate as cm
from sstar.cal_match_rate import cal_match_pct
from sstar.utils import read_data  # kept from original, even if unused now


@pytest.fixture
def data():
    pytest.ref_ind_file = "./examples/data/ind_list/ref.ind.list"
    pytest.tgt_ind_file = "./examples/data/ind_list/tgt.ind.list"
    pytest.src_ind_file = "./examples/data/ind_list/nean.ind.list"
    pytest.vcf = "./tests/data/test.match.rate.data.vcf"
    pytest.score_file = "./tests/results/test.match.rate.score.exp.results"
    pytest.output = "./tests/results/test.match.rate.results"
    pytest.exp_output = "./tests/results/test.match.rate.exp.results"


def test_cal_match_pct_integration(data):
    """Original integration test: make sure we don't break existing behavior."""
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

    with open(pytest.output) as f:
        res = f.readlines()
    with open(pytest.exp_output) as f:
        exp_res = f.readlines()

    assert res == exp_res


@pytest.mark.parametrize(
    "hap1,hap2,expected_str,src_samples",
    [
        # numeric branch: (40 + 60) / 2 = 50
        (40.0, 60.0, "50.0", ["src1", "src2"]),
        # NA branch: both haplotypes NA -> hap_match_pct stays 'NA'
        ("NA", "NA", "NA", ["srcA"]),
    ],
)
def test__cal_match_pct_ind_branches(monkeypatch, hap1, hap2, expected_str, src_samples):
    """
    Directly test the inner function _cal_match_pct_ind:
    - numeric case: average of two haplotypes
    - NA case: keeps 'NA'
    """

    def fake_cal_matchpct(
        chr_name,
        mapped_intervals,
        tgt_data,
        src_data,
        tgt_ind_index,
        src_ind_index,
        hap_index,
        win_start,
        win_end,
        sample_size,
    ):
        return ["dummy", hap1 if hap_index == 0 else hap2]

    # Patch sstar.cal_match_rate.cal_matchpct
    monkeypatch.setattr(cm, "cal_matchpct", fake_cal_matchpct)

    # Minimal score line: only chr, start, end, sample, and last column are used
    data = ["chr1\t100\t200\ttgtSample\t.\t.\t0,10"]

    res = cm._cal_match_pct_ind(
        data=data,
        tgt_ind_index=0,
        mapped_intervals={},
        tgt_data={},
        src_data={},
        src_samples=src_samples,
        sample_size=len(src_samples),
    )

    # One result per source sample
    assert len(res) == len(src_samples)

    fields = res[0].split("\t")
    assert fields[0] == "chr1"
    assert fields[1] == "100"
    assert fields[2] == "200"
    assert fields[3] == "tgtSample"
    assert fields[4] == expected_str
    assert fields[5] == src_samples[0]

