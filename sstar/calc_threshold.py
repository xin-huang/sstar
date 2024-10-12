# Copyright 2025 Xin Huang
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


import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.gam.api import GLMGam, BSplines
from typing import Any, Dict, List, Tuple

from pygam import LinearGAM, s, te


# @profile
def calc_threshold(
    simulated_data: str,
    score_file: str,
    recomb_rate: float,
    recomb_map: str,
    quantile: float,
    output: str,
    k: int
) -> None:
    """
    Calculate significant S* scores from simulated data.

    Parameters
    ----------
    simulated_data : str
        Name of the file containing simulated data.
    score_file : str
        Name of the file containing S* scores.
    recomb_rate : float
        Uniform recombination rate assumed across the genome.
    recomb_map : str
        Name of the file containing recombination maps.
    quantile : float
        Quantile for estimating significant S* scores.
    output : str
        Name of the output file.
    k : int
        Dimension(s) of the bases used to represent the smooth term.
    """
    if (recomb_rate != None) and (recomb_map != None):
        fit_lr = True
    else:
        fit_lr = False

    # Get model, smoother, and data ranges
    gam_tuple = _build_gam_model(simulated_data, k, fit_lr)
    res = _predict_res(gam_tuple, score_file, recomb_rate, recomb_map, quantile, fit_lr)

    header = "chrom\tstart\tend\tsample\tS*_score\texpected_S*_score\tlocal_recomb_rate\tquantile\tsignificant"
    with open(output, "w") as o:
        o.write(header + "\n")
        o.write("\n".join(res) + "\n")


def _build_gam_model(simulated_data: str, k: int, fit_lr: bool) -> Tuple[LinearGAM, None, Dict[str, Tuple[float, float]]]:
    """
    Build a GAM model using pygam (replicating mgcv behavior).

    Parameters
    ----------
    simulated_data : str
        Path to the file containing simulated data.
    k : int
        Number of basis functions for smooth terms.
    fit_lr : bool
        Whether to include local recombination rate in the model.

    Returns
    -------
    model : LinearGAM
        Fitted GAM model.
    smoother : None
        Unused; for compatibility.
    data_ranges : dict
        Min/max values for each feature, used in prediction clipping.
    """
    df = pd.read_csv(simulated_data, sep="\t")

    if fit_lr:
        X = df.iloc[:, [1, 3, 2]].values  # snps, lr, q
        term = te(0, 1, 2, n_splines=k)
    else:
        X = df.iloc[:, [1, 2]].values  # snps, q
        term = te(0, 1, n_splines=k)

    y = df.iloc[:, 0].values  # s_star
    model = LinearGAM(term).fit(X, y)

    data_ranges = {
        i: (X[:, i].min(), X[:, i].max()) for i in range(X.shape[1])
    }

    return model, None, data_ranges


def _predict_res(
    gam_tuple: Tuple[LinearGAM, None, Dict[int, Tuple[float, float]]],
    score_file: str,
    recomb_rate: float,
    recomb_map: str,
    quantile: float,
    fit_lr: bool
) -> List[str]:
    """
    Predict significant S* scores using a fitted pygam model (structure-free, index-based).

    Parameters
    ----------
    gam_tuple : tuple
        Fitted GAM model, None, and feature ranges.
    score_file : str
        Path to file containing observed S* scores.
    recomb_rate : float
        Fallback recombination rate (used if recomb_map is None).
    recomb_map : str
        Path to recombination map (or None).
    quantile : float
        Quantile used in significance decision.
    fit_lr : bool
        Whether to use local recombination rate as a feature.

    Returns
    -------
    res : list of str
        Tab-separated result lines.
    """
    model, _, data_ranges = gam_tuple

    if recomb_map is not None:
        recomb_map = _read_recomb_map(recomb_map)

    data, meta_data = _read_score_file(score_file, recomb_rate, recomb_map, quantile, fit_lr)

    thresholds = {}
    for s in data:
        if fit_lr:
            X = np.stack([data[s]['snps'], data[s]['lr'], data[s]['q']], axis=1)
        else:
            X = np.stack([data[s]['snps'], data[s]['q']], axis=1)

        for i in range(X.shape[1]):
            X[:, i] = np.clip(X[:, i], *data_ranges[i])

        thresholds[s] = model.predict(X)

    res = []
    for s in meta_data:
        for i in range(len(meta_data[s])):
            score = float(meta_data[s][i][4])
            null_score = round(thresholds[s][i], 6)
            sig = "True" if score > null_score else "False"

            line = "\t".join(meta_data[s][i])
            line += "\t" + str(null_score)
            line += "\t" + str(data[s]["lr"][i]) if fit_lr else "\t" + str(recomb_rate)
            line += "\t" + str(data[s]["q"][i])
            line += "\t" + sig

            res.append(line)

    return res


# @profile
def _build_gam_model_sm(
    simulated_data: str,
    k: int,
    fit_lr: bool
) -> Tuple[Any, Any, Dict[str, Tuple[float, float]]]:
    """
    Helper function for building a generalized additive model with statsmodels.
    
    This builds a model to predict expected S* scores under a demographic model
    without introgression, based on the number of SNPs, quantile, and optionally
    local recombination rate.

    Parameters
    ----------
    simulated_data : str
        Name of the file containing simulated data.
    k : int
        Dimension(s) of the bases used to represent the smooth term.
    fit_lr : bool
        If True, build a GAM that includes local recombination rate.

    Returns
    -------
    tuple
        A tuple containing (model_results, smoother, data_ranges), where:
        - model_results : statsmodels object
            Fitted statsmodels GAM model.
        - smoother : BSplines
            BSplines smoother used for basis expansion.
        - data_ranges : dict
            Dictionary with min/max values of predictors for clipping during prediction.

    Notes
    -----
    The model stores data ranges to ensure prediction data stays within the
    training data bounds, which is necessary for statsmodels spline implementation.
    """
    data = pd.read_csv(simulated_data, sep="\t")

    s_star = data.iloc[:, 0]
    snps = data.iloc[:, 1]
    q = data.iloc[:, 2]

    data_ranges = {}
    data_ranges["snps"] = (np.min(snps), np.max(snps))
    data_ranges["q"] = (np.min(q), np.max(q))

    if fit_lr:
        lr = data.iloc[:, 3]
        data_ranges["lr"] = (np.min(lr), np.max(lr))
        X = np.column_stack([snps, lr, q])
        smoother = BSplines(X, df=[k, k, k], degree=[3, 3, 3])
    else:
        X = np.column_stack([snps, q])
        smoother = BSplines(X, df=[k, k], degree=[3, 3])

    family = sm.families.Gaussian()
    gam_model = GLMGam(s_star, smoother=smoother, family=family)
    results = gam_model.fit()

    return (results, smoother, data_ranges)


# @profile
def _predict_res_sm(
    gam_tuple: Tuple[Any, Any, Dict[str, Tuple[float, float]]],
    score_file: str,
    recomb_rate: float,
    recomb_map: str,
    quantile: float,
    fit_lr: bool
) -> List[List[Any]]:
    """
    Helper function to predict significant S* scores using statsmodels GAM.

    Uses the fitted GAM to predict expected S* scores under a null model,
    then compares observed S* scores to these predictions to determine
    significance (observed > expected).

    Parameters
    ----------
    gam_tuple : tuple
        Tuple containing (model_results, smoother, data_ranges).
    score_file : str
        Name of the file containing S* scores.
    recomb_rate : float
        Uniform recombination rate assumed across the genome.
    recomb_map : str
        Name of the file containing recombination maps.
    quantile : float
        Quantile for estimating significant S* scores.
    fit_lr : bool
        If True, model includes local recombination rate.

    Returns
    -------
    res : list of list
        List containing results for output with columns:
        chrom, start, end, sample, S*_score, expected_S*_score,
        local_recomb_rate, quantile, significant.
    """
    gam_results, smoother, data_ranges = gam_tuple

    if recomb_map != None:
        recomb_map = _read_recomb_map(recomb_map)
    data, meta_data = _read_score_file(
        score_file, recomb_rate, recomb_map, quantile, fit_lr
    )

    thresholds = dict()
    for s in data.keys():
        snps_clipped = np.clip(
            data[s]["snps"], data_ranges["snps"][0], data_ranges["snps"][1]
        )
        q_clipped = np.clip(data[s]["q"], data_ranges["q"][0], data_ranges["q"][1])

        if fit_lr:
            lr_clipped = np.clip(
                data[s]["lr"], data_ranges["lr"][0], data_ranges["lr"][1]
            )
            X_pred = np.column_stack([snps_clipped, lr_clipped, q_clipped])
        else:
            X_pred = np.column_stack([snps_clipped, q_clipped])

        X_trans = smoother.transform(X_pred)
        thresholds[s] = gam_results.predict(exog=X_trans)

    res = []
    for s in meta_data.keys():
        for i in range(len(meta_data[s])):
            score = float(meta_data[s][i][4])
            null_score = round(thresholds[s][i], 6)
            if score > null_score:
                sig = "True"
            else:
                sig = "False"

            line = "\t".join(meta_data[s][i])
            line += "\t" + str(null_score)
            line += "\t" + str(data[s]["lr"][i])
            line += "\t" + str(data[s]["q"][i])
            line += "\t" + sig

            res.append(line)
    return res


# @profile
def _read_score_file(
    score_file: str,
    recomb_rate: float,
    recomb_map: str,
    quantile: float,
    fit_lr: bool
) -> Tuple[Dict, Dict]:
    """
    Helper function for reading files from `sstar score`.

    Parameters
    ----------
    score_file : str
        Name of the file containing S* scores.
    recomb_rate : float
        Uniform recombination rate assumed across the genome.
    recomb_map : str
        Name of the file containing recombination maps.
    quantile : float
        Quantile for estimating significant S* scores.
    fit_lr : bool
        If True, building a GAM with local recombination rate.

    Returns
    -------
    data : dict
        Dictionary containing data for predicting significant S* scores.
    meta_data : dict
        Dictionary containing meta data for output.
    """
    data = dict()
    meta_data = dict()
    with open(score_file, "r") as f:
        f.readline()
        for line in f.readlines():
            line = line.rstrip()
            element = line.split("\t")
            if element[6] == "NA":
                continue
            else:
                chr_name = element[0]
                win_start = element[1]
                win_end = element[2]
                sample = element[3]
                score = element[4]
                total_snp_num = element[5]
                if fit_lr:
                    if recomb_map != None:
                        key = chr_name + ":" + win_start + "-" + win_end
                        if key in recomb_map.keys():
                            local_recomb_rate = recomb_map[key]
                        else:
                            continue
                    else:
                        local_recomb_rate = recomb_rate
                else:
                    local_recomb_rate = np.nan

                if sample not in data.keys():
                    data[sample] = dict()
                    data[sample]["snps"] = list()
                    data[sample]["lr"] = list()
                    data[sample]["q"] = list()
                if sample not in meta_data.keys():
                    meta_data[sample] = list()
                data[sample]["snps"].append(float(total_snp_num))
                data[sample]["lr"].append(local_recomb_rate)
                data[sample]["q"].append(quantile)
                meta_data[sample].append([chr_name, win_start, win_end, sample, score])

    return data, meta_data


# @profile
def _read_recomb_map(recomb_map_file: str) -> Dict:
    """
    Helper function for reading recombination maps from files.

    Parameters
    ----------
    recomb_map_file : str
        Name of the file containing recombination maps.

    Returns
    -------
    recomb_map : dict
        Dictionary containing recombination maps.
    """
    recomb_map = dict()
    with open(recomb_map_file, "r") as f:
        for line in f.readlines():
            line = line.rstrip()
            element = line.split("\t")
            key = element[0] + ":" + element[1] + "-" + element[2]
            if key not in recomb_map.keys():
                recomb_map[key] = float(element[3])

    return recomb_map
