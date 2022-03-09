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

import gzip
import numpy as np
import rpy2.robjects as ro
from rpy2.robjects import FloatVector
from rpy2.robjects.packages import importr

def cal_likelihood(sim_null, sim_alt1, sim_alt2, real_data, threshold, output, grid_size, bdw_null, bdw_alt1, bdw_alt2, step_size):
    """
    Description:
        Calculate posterior probabilities under the no introgression model, the model with introgression only from the src1 population,
        and the model with introgression only from the src2 population.

    Arguments:
        sim_null list: Names of the files containing p-values from simulated data under the no introgression model (the null model).
        sim_alt1 list: Names of the files containing p-values from simulated data under the model with introgression only from the src1 population (the alternative model 1).
        sim_alt2 list: Names of the files containing p-values from simulated data under the model with introgression only from the src2 population (the alternative model 2).
        real_data list: Names of the files containing p-values from real data.
        threshold str: Name of the file containing windows with statistically significant S* scores.
        output str: Name of the output file.
        grid_size int: Grid size for kernel density estimation.
        bdw_null float: Bandwidth for kernel density estimation with the null model.
        bdw_alt1 float: Bandwidth for kernel density estimation with the alternative model 1.
        bdw_alt2 float: Bandwidth for kernel density estimation with the alternative model 2.
        step_size float: Step size for prior probabilities.
    """
    sim_null_pvals = _read_pvals(sim_null, None)
    sim_alt1_pvals = _read_pvals(sim_alt1, None)
    sim_alt2_pvals = _read_pvals(sim_alt2, None)
    real_data_pvals = _read_pvals(real_data, threshold)

    sim_null_logits1, sim_null_logits2, sim_null_keys, sim_null_metadata = _logit(sim_null_pvals)
    sim_alt1_logits1, sim_alt1_logits2, sim_alt1_keys, sim_alt1_metadata = _logit(sim_alt1_pvals)
    sim_alt2_logits1, sim_alt2_logits2, sim_alt2_keys, sim_alt2_metadata = _logit(sim_alt2_pvals)
    real_data_logits1, real_data_logits2, real_data_keys, real_data_metadata = _logit(real_data_pvals)

    lims = [min(sim_null_logits1+sim_null_logits2+sim_alt1_logits1+sim_alt1_logits2+sim_alt2_logits1+sim_alt2_logits2)-5, 
            max(sim_null_logits1+sim_null_logits2+sim_alt1_logits1+sim_alt1_logits2+sim_alt2_logits1+sim_alt2_logits2)+5]
    lims = lims + lims

    null_density = _kde2d(px=sim_null_logits1, py=sim_null_logits2, n=grid_size, h=bdw_null, lims=lims, x=real_data_logits1, y=real_data_logits2)
    alt1_density = _kde2d(px=sim_alt1_logits1, py=sim_alt1_logits2, n=grid_size, h=bdw_alt1, lims=lims, x=real_data_logits1, y=real_data_logits2)
    alt2_density = _kde2d(px=sim_alt2_logits1, py=sim_alt2_logits2, n=grid_size, h=bdw_alt2, lims=lims, x=real_data_logits1, y=real_data_logits2)

    max_ll = dict()
    max_ll['ll'] = -np.inf
    for prior_alt1 in np.arange(1e-5, 1-1e-5+step_size, step_size):
        for prior_alt2 in np.arange(1e-5, 1-1e-5+step_size, step_size):
            ll = _ll(prior_alt1, prior_alt2, null_density, alt1_density, alt2_density, len(real_data_logits1))
            if ll['ll'] == 'NA': continue
            if ll['ll'] > max_ll['ll']: max_ll = ll

    header = "chrom\tstart\tend\tsample\thap_index\tS*_start\tS*_end\tpost_null\tpost_src1\tpost_src2"
    with open(output, 'w') as o:
        o.write(header+"\n")
        size = len(max_ll['post_null'])
        ll = round(max_ll['ll'], 8)
        for k in real_data_pvals.keys():
            if k in real_data_keys:
                index = real_data_keys.index(k)
                s_star_start, s_star_end = real_data_metadata[index].split(":")
                chr_name, win_start, win_end, sample, hap_index = k.split(":")

                src1_pval = real_data_pvals[k]['pval1']
                src2_pval = real_data_pvals[k]['pval2']
                post_null = round(max_ll['post_null'][index], 8)
                post_alt1 = round(max_ll['post_alt1'][index], 8)
                post_alt2 = round(max_ll['post_alt2'][index], 8)

                line = f'{chr_name}\t{win_start}\t{win_end}\t{sample}\t'
                line += f'{hap_index}\t{s_star_start}\t{s_star_end}\t{post_null}\t{post_alt1}\t{post_alt2}'
                o.write(line+"\n")

def _kde2d(px, py, n, h, lims, x, y):
    """
    Description:
        Helper function for 2d kernel density estimation with R packages MASS and spatstat.

    Arguments:
        px list: P-values from simulated data for kernel density estimation in dimension x.
        py list: P-values from simulated data for kernel density estimation in dimension y.
        n int: Grid size for kernel density estimation.
        h float: Bandwidth for kernel density estimation.
        lims list: Limits for kernel density estimation.
        x list: P-values from real data for interpolation in dimension x.
        y list: P-values from real data for interpolation in dimension y.

    Returns:
        density np.ndarray: Interpolated results.
    """
    base = importr('base')
    mass = importr('MASS')
    spatstat_geom = importr('spatstat.geom')

    px = FloatVector(px)
    py = FloatVector(py)
    x = FloatVector(x)
    y = FloatVector(y)
    n = FloatVector([n])
    h = FloatVector([h])
    lims = FloatVector(lims)

    kde_res = mass.kde2d(px, py, n=n, h=h, lims=lims)
    im_res = spatstat_geom.im(base.t(kde_res[2]), xcol=kde_res[0], yrow=kde_res[1])
    density = np.array(spatstat_geom.interp_im(im_res, x, y))

    return density

def _ll(prior_alt1, prior_alt2, null_density, alt1_density, alt2_density, data_size):
    """
    Description:
        Helper function for calculating likelihoods and posterior probabilities.

    Arguments:
        prior_alt1 float: Prior probability for the alternative model 1.
        prior_alt2 float: Prior probability for the alternative model 2.
        null_density np.ndarray: Interpolated results under the null model.
        alt1_density np.ndarray: Interpolated results under the alternative model 1.
        alt2_density np.ndarray: Interpolated results under the alternative model 2.
        data_size int: Size of the data for estimating the posterior probabilities.

    Returns:
        res dict: Dictionary containing likelihoods and posterior probabilities for analysis.
    """
    res = dict()
    res['ll'] = 'NA'
    res['post_null'] = 'NA'
    res['post_alt1'] = 'NA'
    res['post_alt2'] = 'NA'
    res['prior_alt1'] = 'NA'
    res['prior_alt2'] = 'NA'

    if prior_alt1 + prior_alt2 >= 1: return res
    
    prior_null = 1 - prior_alt1 - prior_alt2
    pz_null = null_density * prior_null
    pz_alt1 = alt1_density * prior_alt1
    pz_alt2 = alt2_density * prior_alt2
    pz_null[pz_null == 0] = 1 / data_size
    pz_alt1[pz_alt1 == 0] = 1 / data_size
    pz_alt2[pz_alt2 == 0] = 1 / data_size
    pz = pz_null + pz_alt1 + pz_alt2

    res['ll'] = np.sum(np.log(pz))
    res['post_null'] = pz_null / pz
    res['post_alt1'] = pz_alt1 / pz
    res['post_alt2'] = pz_alt2 / pz
    res['prior_alt1'] = prior_alt1
    res['prior_alt2'] = prior_alt2

    return res

def _read_pvals(filenames, threshold_file):
    """
    Description:
        Helper function for reading p-values from files.

    Arguments:
        filenames list: Names of the files containing p-values.
        threshold_file str: Name of the file containing windows with statistically significant S* scores.

    Returns:
        pvals dict: Dictionary containing p-values and other information for analysis.
    """
    pvals = dict()
    thresholds = dict()

    if threshold_file is not None:
        with open(threshold_file, 'r') as f: 
            f.readline()
            for line in f.readlines():
                line = line.rstrip()
                elements = line.split("\t")
                chr_name = elements[0]
                win_start = elements[1]
                win_end = elements[2]
                sample = elements[3]
                key = f'{chr_name}:{win_start}:{win_end}:{sample}'
                key1 = key + ":1"
                key2 = key + ":2"
                if elements[-1] == 'False':
                    thresholds[key1] = False
                    thresholds[key2] = False
    
    f = gzip.open(filenames[0], 'rt')
    try:
        f.readline()
    except:
        f.close()
        f = open(filenames[0], 'r')
        f.readline()

    for line in f.readlines():
        line = line.rstrip()
        elements = line.split("\t")
        chr_name = elements[0]
        win_start = elements[1]
        win_end = elements[2]
        sample = elements[3]
        pvalue = elements[4]
        hap_index = elements[6]
        s_star_start = elements[7]
        s_star_end = elements[8]
        if pvalue == 'NA': continue
        key = f'{chr_name}:{win_start}:{win_end}:{sample}:{hap_index}'
        metadata = ":".join(elements[9:11]) 
        if key in thresholds.keys(): continue
        if key not in pvals.keys(): pvals[key] = dict()
        pvals[key]['metadata'] = f'{s_star_start}:{s_star_end}'
        pvals[key]['pval1'] = pvalue
    f.close()

    f = gzip.open(filenames[1], 'rt')
    try:
        f.readline()
    except:
        f.close()
        f = open(filenames[1], 'r')
        f.readline()
    
    for line in f.readlines():
        line = line.rstrip()
        elements = line.split("\t")
        chr_name = elements[0]
        win_start = elements[1]
        win_end = elements[2]
        sample = elements[3]
        pvalue = elements[4]
        hap_index = elements[6]
        s_star_start = elements[7]
        s_star_end = elements[8]
        if pvalue == 'NA': continue
        key = f'{chr_name}:{win_start}:{win_end}:{sample}:{hap_index}'
        if key in thresholds.keys(): continue
        if key in pvals.keys(): pvals[key]['pval2'] = pvalue
    f.close()

    return pvals

def _logit(pvals):
    """
    Description:
        Helper function for converting p-values into logits.

    Arguments:
        pvals list: List containing p-values for convertion.

    Returns:
        logits1 list: List containing logits from p-values calculating with the src1 population.
        logits2 list: List containing logits from p-values calculating with the src2 population.
        keys list: List containing sample information and position information.
        metadata list: List containing statistics from haplotypes in windows.
    """
    logits1 = list()
    logits2 = list()
    keys = list()
    metadata = list()

    for k in pvals.keys():
        if 'pval1' not in pvals[k].keys(): continue
        if 'pval2' not in pvals[k].keys(): continue
        if (pvals[k]['pval1'] != 'NA') and (pvals[k]['pval2'] != 'NA'):
            pval1 = float(pvals[k]['pval1'])
            pval2 = float(pvals[k]['pval2'])
            if (pval1 != 0) and (pval1 != 1) and (pval2 != 0) and (pval2 != 1):
                logit1 = np.log(pval1/(1-pval1))
                logit2 = np.log(pval2/(1-pval2))
                logits1.append(logit1)
                logits2.append(logit2)
                keys.append(k)
                metadata.append(pvals[k]['metadata'])

    return logits1, logits2, keys, metadata
