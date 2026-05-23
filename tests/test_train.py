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


import pandas as pd
from sstar.train import train


def test_train(tmp_path):
    output_file = tmp_path / "check.qr.joblib"

    train(
            demes="tests/data/ArchIE_3D19_wo_intro.yaml",
            config="tests/data/test.config.yaml",
            output=output_file,
        )

    train_features_file = tmp_path / "check.qr.training.features.tsv"
    df = pd.read_csv(train_features_file, sep="\t")
    df_expected = pd.read_csv(
            "tests/exp_results/test.train.expected.features.tsv", sep="\t"
        )

    pd.testing.assert_frame_equal(
            df,
            df_expected,
            check_dtype=False,
            check_exact=False,
            rtol=1e-6,
            atol=1e-8,
        )

    assert output_file.exists(), f"{output_file} was not created"
    assert output_file.is_file(), f"{output_file} is not a file"

