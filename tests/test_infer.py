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
from sstar.infer import infer
from sstar.quantile_regression import infer as qr_infer


def test_infer(tmp_path):
    feat_file = tmp_path / "check.features.tsv"
    pred_file = tmp_path / "check.pred.tsv"
    tract_file = tmp_path / "check.inferred.tracts.bed"

    infer(
        model="tests/data/test.qr.model.onnx",
        config="tests/data/test.config.yaml",
        feat_file=feat_file,
        pred_file=pred_file,
        tract_file=tract_file,
        match_bonus=5000,
        max_mismatch=5,
        mismatch_penalty=-10000,
    )

    res_files = [feat_file, pred_file, tract_file]
    exp_res_files = [
        "tests/exp_results/test.infer.expected.features.tsv",
        "tests/exp_results/test.infer.expected.pred.tsv",
        "tests/exp_results/test.infer.expected.inferred.tracts.bed",
    ]

    for f1, f2 in zip(res_files, exp_res_files):
        df = pd.read_csv(f1, sep="\t")
        df_expected = pd.read_csv(f2, sep="\t")

        pd.testing.assert_frame_equal(
            df,
            df_expected,
            check_dtype=False,
            check_exact=False,
            rtol=1e-6,
            atol=1e-8,
        )


def test_quantile_regression_infer_skips_missing_sstar_score(tmp_path):
    features_file = tmp_path / "features.tsv"
    pred_file = tmp_path / "pred.tsv"
    tract_file = tmp_path / "tracts.bed"

    features = pd.DataFrame(
        {
            "Chromosome": ["21", "21"],
            "Start": [1, 50001],
            "End": [50000, 100000],
            "Sample": ["ind1_1", "ind1_1"],
            "S*_score": [pd.NA, 67770.0],
            "Region_ind_SNP_number": [10, 10],
        }
    )
    features.to_csv(features_file, sep="\t", index=False, na_rep="NA")

    qr_infer(
        data=features_file,
        model="tests/data/test.qr.model.onnx",
        output=pred_file,
        bed_file=tract_file,
    )

    predictions = pd.read_csv(pred_file, sep="\t")

    assert pd.isna(predictions.loc[0, "Predicted_S*_score"])
    assert pd.notna(predictions.loc[1, "Predicted_S*_score"])
