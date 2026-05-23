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


import inspect, joblib, os
import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingRegressor


@staticmethod
def train(
    data: str,
    output: str,
    **model_params,
) -> None:
    """ """
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
    """ """
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
