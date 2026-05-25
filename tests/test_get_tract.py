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
import pandas as pd
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


def test_get_tract(data, tmp_path):
    output_prefix = str(tmp_path) + "/test.tract"
    output = str(tmp_path) + "/test.tract.bed"

    get_tract(
        threshold_file=pytest.threshold,
        match_pct_files=None,
        output_prefix=output_prefix,
        diff=0,
    )

    df = pd.read_csv(output, sep="\t")
    df_expected = pd.read_csv(pytest.exp_bed, sep="\t")

    pd.testing.assert_frame_equal(
        df,
        df_expected,
        check_dtype=False,
        check_exact=False,
        rtol=1e-6,
        atol=1e-8,
    )


def test_get_tract_with_src(data, tmp_path):
    output_prefix = str(tmp_path) + "/test.tract.with.src.match.rate"
    output = str(tmp_path) + "/test.tract.with.src.match.rate.bed"

    get_tract(
        threshold_file=pytest.threshold,
        match_pct_files=[pytest.src1_match_pct],
        output_prefix=output_prefix,
        diff=0,
    )

    df = pd.read_csv(output, sep="\t")
    df_expected = pd.read_csv(pytest.exp_bed_with_src, sep="\t")

    pd.testing.assert_frame_equal(
        df,
        df_expected,
        check_dtype=False,
        check_exact=False,
        rtol=1e-6,
        atol=1e-8,
    )


def test_get_tract_with_multi_src(data, tmp_path):
    output_prefix = str(tmp_path) + "/test.tract"
    output1 = str(tmp_path) + "/test.tract.src1.bed"
    output2 = str(tmp_path) + "/test.tract.src2.bed"

    get_tract(
        threshold_file=pytest.threshold,
        match_pct_files=[pytest.src1_match_pct, pytest.src2_match_pct],
        output_prefix=output_prefix,
        diff=0,
    )

    for f1, f2 in zip([output1, output2], [pytest.exp_src1_bed, pytest.exp_src2_bed]):
        try:
            df = pd.read_csv(f1, sep="\t")
        except pd.errors.EmptyDataError:
            df = pd.DataFrame()

        try:
            df_expected = pd.read_csv(f2, sep="\t")
        except pd.errors.EmptyDataError:
            df_expected = pd.DataFrame()

        pd.testing.assert_frame_equal(
            df,
            df_expected,
            check_dtype=False,
            check_exact=False,
            rtol=1e-6,
            atol=1e-8,
        )
