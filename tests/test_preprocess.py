# Copyright 2026 Xin Huang
#
# GNU General Public License v3.0
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, please see
#
#    https://www.gnu.org/licenses/gpl-3.0.en.html

import pytest
import pandas as pd
from sstar.preprocess import preprocess


def test_preprocess(tmp_path):
    output_file = tmp_path / "check.sstar.tsv"

    preprocess(
        vcf_file="tests/data/test.vcf",
        chr_name="21",
        ref_ind_file="tests/data/test.ref.ind.list",
        tgt_ind_file="tests/data/test.tgt.ind.list",
        win_len=50000,
        win_step=10000,
        feature_config_file="tests/data/sstar.yaml",
        output_file=output_file,
        nprocess=1,
        is_phased=False,
    )

    df = pd.read_csv(output_file, sep="\t")
    df_expected = pd.read_csv("tests/data/test.sstar.expected.tsv", sep="\t")

    pd.testing.assert_frame_equal(
        df,
        df_expected,
        check_dtype=False,
        check_exact=False,
        rtol=1e-6,
        atol=1e-8, 
    )
