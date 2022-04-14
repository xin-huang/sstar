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
import rpy2.robjects as ro
from rpy2.robjects import FloatVector, Formula, pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects.conversion import localconverter

#@profile
def cal_threshold(simulated_data, score_file, recomb_rate, recomb_map, quantile, output, k):
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
    """
    if (recomb_rate != None) and (recomb_map != None): fit_lr = True
    else: fit_lr = False

    gam = _build_gam_model(simulated_data, k, fit_lr)
    res = _predict_res(gam, score_file, recomb_rate, recomb_map, quantile, fit_lr)
    header = 'chrom\tstart\tend\tsample\tS*_score\texpected_S*_score\tlocal_recomb_rate\tquantile\tsignificant'
    with open(output, 'w') as o:
        o.write(header+"\n")
        o.write("\n".join(res)+"\n")

#@profile
def _build_gam_model(simulated_data, k, fit_lr):
    """
    Description:
        Helper function for building a generalized additive model with R package MGCV.

    Arguments:
        simulated_data str: Name of the file containing simulated data.
        k int: Dimension(s) of the bases used to represent the smooth term.
        fit_lr bool: If True, building a gam with local recombination rate.

    Returns:
        gam rpy2.robjects: Generalized additive model for prediction.
    """
    mgcv = importr('mgcv')

    data = pd.read_csv(simulated_data, sep="\t")

    s_star = FloatVector(data.iloc[:,0].values)
    snps = FloatVector(data.iloc[:,1].values)
    q = FloatVector(data.iloc[:,2].values)

    ro.globalenv['s_star'] = s_star
    ro.globalenv['snps'] = snps
    ro.globalenv['q'] = q

    if fit_lr: 
        lr = FloatVector(data.iloc[:,3].values)
        ro.globalenv['lr'] = lr
        fmla = f's_star ~ te(snps, lr, q, k={k})'
    else: fmla = f's_star ~ te(snps, q, k={k})'

    fmla = Formula(fmla)

    gam = mgcv.gam(fmla)

    return gam

#@profile
def _predict_res(gam, score_file, recomb_rate, recomb_map, quantile, fit_lr):
    """
    Description:
        Helper function to predict significant S* scores.

    Arguments:
        gam rpy2.robjects: Generalized additive model for prediction.
        score_file str: Name of the file containing S* scores.
        recomb_rate float: Uniform recombination rate assumed across the genome.
        recomb_map str: Name of the file containing recombination maps.
        quantile float: Quantile for estimating significant S* scores.
        fit_lr bool: If True, building a gam with local recombination rate.

    Returns:
        res list: List containing results for output.
    """
    stats = importr('stats')

    if recomb_map != None: recomb_map = _read_recomb_map(recomb_map)
    data, meta_data = _read_score_file(score_file, recomb_rate, recomb_map, quantile, fit_lr)

    thresholds = dict()
    for s in data.keys():
        if fit_lr: pd_df = pd.DataFrame({'snps': data[s]['snps'], 'lr': data[s]['lr'], 'q': data[s]['q']})
        else: pd_df = pd.DataFrame({'snps': data[s]['snps'], 'q': data[s]['q']})

        with localconverter(ro.default_converter + pandas2ri.converter):
            r_from_pd_df = ro.conversion.py2rpy(pd_df)
        
        thresholds[s] = np.array(stats.predict(gam, r_from_pd_df))

    res = []
    for s in meta_data.keys():
        for i in range(len(meta_data[s])):
            score = float(meta_data[s][i][4])
            null_score = round(thresholds[s][i], 6)
            if score > null_score: sig = 'True'
            else: sig = 'False'

            line = "\t".join(meta_data[s][i])
            line += "\t" + str(null_score)
            line += "\t" + str(10**data[s]['lr'][i])
            line += "\t" + str(data[s]['q'][i])
            line += "\t" + sig

            res.append(line)

    return res

#@profile
def _read_score_file(score_file, recomb_rate, recomb_map, quantile, fit_lr):
    """
    Description:
        Helper function for reading files from `sstar score`.

    Arguments:
        score_file str: Name of the file containing S* scores.
        recomb_rate float: Uniform recombination rate assumed across the genome. 
        recomb_map str: Name of the file containing recombination maps.
        quantile float: Quantile for estimating significant S* scores.
        fit_lr bool: If True, building a gam with local recombination rate.

    Returns:
        data dict: Dictionary containing data for predicting significat S* scores.
        meta_data dict: Dictionary containing meta data for output.
    """
    data = dict()
    meta_data = dict()
    with open(score_file, 'r') as f:
        f.readline()
        for line in f.readlines():
            line = line.rstrip()
            element = line.split("\t")
            if element[6] == 'NA': continue
            else:
                chr_name = element[0]
                win_start = element[1]
                win_end = element[2]
                sample = element[3]
                score = element[4]
                total_snp_num = element[5]
                if fit_lr:
                    if recomb_map != None: 
                        key = chr_name+":"+win_start+"-"+win_end
                        if key in recomb_map.keys(): local_recomb_rate = np.log10(recomb_map[key])
                        else: continue
                    else:
                        local_recomb_rate = np.log10(recomb_rate)
                else: local_recomb_rate = np.nan

                if sample not in data.keys(): 
                    data[sample] = dict()
                    data[sample]['snps'] = list()
                    data[sample]['lr'] = list()
                    data[sample]['q'] = list()
                if sample not in meta_data.keys(): meta_data[sample] = list()
                data[sample]['snps'].append(float(total_snp_num))
                data[sample]['lr'].append(local_recomb_rate)
                data[sample]['q'].append(quantile)
                meta_data[sample].append([chr_name, win_start, win_end, sample, score])

    return data, meta_data

#@profile
def _read_recomb_map(recomb_map_file):
    """
    Description:
        Helper function for reading recombination maps from files.

    Arguments:
        recomb_map_file str: Name of the file containing recombination maps.

    Returns:
        recomb_map dict: Dictionary containing recombination maps.
    """
    recomb_map = dict()
    with open(recomb_map_file, 'r') as f:
        for line in f.readlines():
            line = line.rstrip()
            element = line.split("\t")
            key = element[0]+":"+element[1]+"-"+element[2]
            if key not in recomb_map.keys(): recomb_map[key] = float(element[3])

    return recomb_map

if __name__ == '__main__':
   cal_threshold(simulated_data="../examples/data/simulated_data/gravel_asn_scale_60k.simulated.data", score_file="../tests/data/test.score", recomb_rate=1.29, recomb_map=None, quantile=0.99, output="threshold.out")
