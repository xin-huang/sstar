# Copyright 2026 Xin Huang and Andrea Koça
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


import inspect, joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingRegressor


@staticmethod
def train(
    data: str,
    output: str,
    **model_params,
) -> None:
    """
    Train a gradient boosting regressor for predicting `S*_score`.

    The model is trained using `Region_ind_SNP_number` as the input feature and
    `S*_score` as the target variable. Rows with missing `S*_score` values are
    removed before training. Unsupported parameters for
    `GradientBoostingRegressor` are ignored.

    Parameters
    ----------
    data : str
        Path to the input tab-separated feature table.
    output : str
        Path to the output joblib file for the trained model.
    **model_params
        Additional keyword arguments passed to `GradientBoostingRegressor`.
    """
    is_allowed = inspect.signature(GradientBoostingRegressor).parameters
    clean_params = {k: v for k, v in model_params.items() if k in is_allowed}

    df = pd.read_csv(data, sep="\t").dropna(subset=["S*_score"])
    x = df[["Region_ind_SNP_number"]].copy()
    y = df["S*_score"].copy()

    model = GradientBoostingRegressor(**clean_params)
    model.fit(x, y)

    joblib.dump(model, output)


@staticmethod
def infer(
    data: str,
    model: str,
    output: str,
    bed_file: str,
) -> None:
    """
    Predict `S*_score` and export regions exceeding the predicted score.

    Predictions are generated for rows with non-missing
    `Region_ind_SNP_number`. Rows with missing feature values keep
    `Predicted_S*_score` as `NA`. The full table is written to `output`.
    Regions where observed `S*_score` is greater than `Predicted_S*_score`
    are written to `bed_file` in BED format. The BED start coordinate is
    converted from 1-based to 0-based by subtracting 1.

    Parameters
    ----------
    data : str
        Path to the input tab-separated feature table.
    model : str
        Path to the trained joblib model.
    output : str
        Path to the output tab-separated table with predictions.
    bed_file : str
        Path to the output BED file containing regions with observed
        `S*_score` greater than `Predicted_S*_score`.
    """
    df = pd.read_csv(data, sep="\t")
    mask = df["Region_ind_SNP_number"].notna()
    df["Predicted_S*_score"] = np.nan
    model = joblib.load(model)
    df.loc[mask, "Predicted_S*_score"] = model.predict(
        df.loc[mask, ["Region_ind_SNP_number"]]
    )

    df.sort_values(by=["Sample", "Chromosome", "Start", "End"]).to_csv(
        output, sep="\t", index=False, na_rep="NA"
    )

    bed = df.loc[
        df["S*_score"] > df["Predicted_S*_score"],
        ["Chromosome", "Start", "End", "Sample"],
    ].copy()

    bed["Start"] = bed["Start"] - 1

    bed.to_csv(bed_file, sep="\t", index=False, header=False)
