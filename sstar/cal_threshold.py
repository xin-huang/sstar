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

import numpy as np
import pandas as pd
import re
import rpy2.robjects as ro
from typing import Any, Dict, List, Optional, Tuple
from rpy2.robjects import FloatVector, Formula, pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects.conversion import localconverter


def cal_threshold(
    simulated_data: str,
    score_file: str,
    recomb_rate: float,
    recomb_map: Optional[str],
    quantile: float,
    output: str,
    k: int,
    phased: bool = False,
) -> None:
    """
    Description:
        Calculate significant S* scores from simulated data.

    Arguments:
        simulated_data str: Name of the file containing simulated data.
        score_file str: Name of the file containing S* scores.
        recomb_rate float: Uniform recombination rate assumed across the genome.
        recomb_map str: Name of the file containing recombination maps.
        quantile float: Quantile for estimating significant S* scores.
        output str: Name of the output file.
        k int: Dimension(s) of the bases used to represent the smooth term.
        phased bool: If True, expect score_file to contain haplotype sample IDs (*_<hap-index>, e.g. *_1).
    """

    # FORCE legacy + correct behavior:
    # If a recombination map is provided, we use lr (this matches the expected file)
    fit_lr = recomb_map is not None

    gam = _build_gam_model(simulated_data, k, fit_lr)
    res = _predict_res(
        gam=gam,
        score_file=score_file,
        recomb_rate=recomb_rate,
        recomb_map=recomb_map,
        quantile=quantile,
        fit_lr=fit_lr,
        phased=phased,
    )

    header = (
        "chrom\tstart\tend\tsample\tS*_score\texpected_S*_score\t"
        "local_recomb_rate\tquantile\tsignificant"
    )
    with open(output, "w") as o:
        o.write(header + "\n")
        o.write("\n".join(res) + "\n")


def _build_gam_model(simulated_data: str, k: int, fit_lr: bool) -> Any:
    """
    Description:
        Helper function for building a generalized additive model with R package MGCV.

    Arguments:
        simulated_data str: Name of the file containing simulated data.
        k int: Dimension(s) of the bases used to represent the smooth term.
        fit_lr bool: If True, build a gam with local recombination rate.

    Returns:
        gam rpy2.robjects: Generalized additive model for prediction.
    """
    mgcv = importr("mgcv")

    data = pd.read_csv(simulated_data, sep="\t")

    s_star = FloatVector(data.iloc[:, 0].values)
    snps = FloatVector(data.iloc[:, 1].values)
    q = FloatVector(data.iloc[:, 2].values)

    if fit_lr:
        lr = FloatVector(data.iloc[:, 3].values)
        fmla_str = f"s_star ~ te(snps, lr, q, k={k})"
    else:
        lr = None
        fmla_str = f"s_star ~ te(snps, q, k={k})"

    fmla = Formula(fmla_str)

    # Put variables into the formula environment (this is what fixes "lr not found")
    fenv = fmla.environment
    fenv["s_star"] = s_star
    fenv["snps"] = snps
    fenv["q"] = q
    if fit_lr:
        fenv["lr"] = lr

    gam = mgcv.gam(fmla)
    return gam


def _predict_res(
    gam: Any,
    score_file: str,
    recomb_rate: float,
    recomb_map: Optional[str],
    quantile: float,
    fit_lr: bool,
    phased: bool,
) -> List[str]:
    """
    Description:
        Helper function to predict significant S* scores.

    Returns:
        res list: List containing results for output.
    """
    stats = importr("stats")

    if recomb_map is not None:
        recomb_map = _read_recomb_map(recomb_map)

    data, meta_data = _read_score_file(
        score_file=score_file,
        recomb_rate=recomb_rate,
        recomb_map=recomb_map,
        quantile=quantile,
        fit_lr=fit_lr,
        phased=phased,
    )

    thresholds = {}
    for s in data.keys():
        if fit_lr:
            pd_df = pd.DataFrame(
                {"snps": data[s]["snps"], "lr": data[s]["lr"], "q": data[s]["q"]}
            )
        else:
            pd_df = pd.DataFrame({"snps": data[s]["snps"], "q": data[s]["q"]})

        with localconverter(ro.default_converter + pandas2ri.converter):
            r_from_pd_df = ro.conversion.py2rpy(pd_df)

        thresholds[s] = np.array(stats.predict(gam, r_from_pd_df))

    res = []
    for s in meta_data.keys():
        for i in range(len(meta_data[s])):
            score = float(meta_data[s][i][4])
            null_score = round(float(thresholds[s][i]), 6)
            sig = "True" if score > null_score else "False"

            line = "\t".join(meta_data[s][i])
            line += "\t" + str(null_score)
            line += "\t" + str(data[s]["lr"][i])
            line += "\t" + str(data[s]["q"][i])
            line += "\t" + sig
            res.append(line)

    return res


def _read_score_file(
    score_file: str,
    recomb_rate: float,
    recomb_map: Optional[Dict[str, float]],
    quantile: float,
    fit_lr: bool,
    phased: bool,
) -> Tuple[Dict[str, Dict[str, List[float]]], Dict[str, List[List[str]]]]:
    """
    Description:
        Helper function for reading files from `sstar score`.

    Returns:
        data dict: Data used for prediction.
        meta_data dict: Meta data for output.
    """
    data = {}
    meta_data = {}

    with open(score_file, "r") as f:
        header = f.readline().rstrip().split("\t")
        col = {name: i for i, name in enumerate(header)}
        for line in f:
            line = line.rstrip()
            if not line:
                continue

            element = line.split("\t")

            # Skip rows where S* SNP number is NA
            if element[col["S*_SNP_number"]] == "NA":
                continue

            chr_name = element[col["chrom"]]
            win_start = element[col["start"]]
            win_end = element[col["end"]]
            sample = element[col["sample"]]
            hap_index = element[col["hap_index"]] if "hap_index" in col else "NA"

            is_hap_row = hap_index != "NA"

            if phased and not is_hap_row:
                raise ValueError(
                    "threshold called with --phased but score file does not contain phased rows "
                    "(expected hap_index column to contain non-NA values)."
                )

            if (not phased) and is_hap_row:
                raise ValueError(
                    "threshold called without --phased but score file contains phased rows "
                    "(hap_index has non-NA values). Use phased quantiles and run threshold with --phased."
                )

            score = element[col["S*_score"]]
            total_snp_num = element[col["region_ind_SNP_number"]]

            if fit_lr:
                if recomb_map is not None:
                    key = f"{chr_name}:{win_start}-{win_end}"
                    if key in recomb_map:
                        local_recomb_rate = recomb_map[key]
                    else:
                        # Keep legacy behavior: if window not in recomb map, skip it
                        continue
                else:
                    # uniform rate (only reached if recomb_rate not in (None, 0))
                    local_recomb_rate = recomb_rate
            else:
                local_recomb_rate = np.nan

            if sample not in data:
                data[sample] = {"snps": [], "lr": [], "q": []}
            if sample not in meta_data:
                meta_data[sample] = []

            data[sample]["snps"].append(float(total_snp_num))
            data[sample]["lr"].append(local_recomb_rate)
            data[sample]["q"].append(quantile)
            meta_data[sample].append([chr_name, win_start, win_end, sample, score])

    return data, meta_data


def _read_recomb_map(recomb_map_file: str) -> Dict[str, float]:
    """
    Description:
        Helper function for reading recombination maps from files.

    Returns:
        recomb_map dict: {"chr:start-end": float(lr)}
    """
    recomb_map = {}
    with open(recomb_map_file, "r") as f:
        for line in f:
            line = line.rstrip()
            if not line.strip():
                continue
            element = line.split("\t")
            key = f"{element[0]}:{element[1]}-{element[2]}"
            if key not in recomb_map:
                recomb_map[key] = float(element[3])
    return recomb_map


if __name__ == "__main__":
    cal_threshold(
        simulated_data="../examples/data/simulated_data/gravel_asn_scale_60k.simulated.data",
        score_file="../tests/data/test.score",
        recomb_rate=1.29,
        recomb_map=None,
        quantile=0.99,
        output="threshold.out",
        k=8,
        phased=False,
    )
