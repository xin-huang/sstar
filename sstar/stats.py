# Apache License Version 2.0
# Copyright 2023 Xin Huang
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
import scipy.stats as sps
from scipy.spatial import distance_matrix


def cal_n_ton(tgt_gt, ploidy):
    """
    Description:
        Calculates individual frequency spetra for samples.

    Arguments:
        tgt_gt numpy.ndarray: Genotype matrix from the target population.
        ploidy int: Ploidy of the genomes.

    Returns:
        spectra numpy.ndarray: Individual frequency spectra for haplotypes.
    """
    mut_num, sample_num = tgt_gt.shape
    iv = np.ones((sample_num, 1))
    counts = (tgt_gt>0)*np.matmul(tgt_gt, iv)
    spectra = np.array([np.bincount(counts[:,idx].astype('int8'), minlength=sample_num*ploidy+1) for idx in range(sample_num)])
    # ArchIE does not count non-segragating sites
    spectra[:,0] = 0

    return spectra


def cal_dist(gt1, gt2):
    """
    Description:
        Calculates pairswise Euclidean distances between two genotype matrixes.

    Arguments:
        gt1 numpy.ndarray: Genotype matrix 1.
        gt2 numpy.ndarray: Genotype matrix 2.

    Returns:
        dists numpy.ndarray: Distances estimated.
    """
    dists = distance_matrix(np.transpose(gt2), np.transpose(gt1))
    dists.sort()

    return dists


def cal_mut_num(ref_gt, tgt_gt, mut_type):
    """
    Description:
        Calculates number of private or total mutations in a sample.

    Arguments:
        ref_gt numpy.ndarray: Genotype matrix from the reference population.
        tgt_gt numpy.ndarray: Genotype matrix from the target population.
        mut_type str: Type of mutations. Private or total.
        
    Returns:
        mut_num numpy.ndarray: Numbers of mutations.
    """
    counts = np.sum(ref_gt, axis=1)
    counts = np.reshape(counts, (counts.shape[0],1))

    if mut_type == 'private': mut_num = np.sum((tgt_gt>0)*(counts==0), axis=0)
    if mut_type == 'total': mut_num = np.sum((tgt_gt>0), axis=0) + np.sum((tgt_gt==0)*(counts>0), axis=0)

    return mut_num


def cal_sstar(tgt_gt, pos, method, match_bonus=5000, max_mismatch=5, mismatch_penalty=-10000):
    """
    Description:
        Calculates sstar scores for a given genotype matrix.

    Arguments:
        tgt_gt numpy.ndarray: Genotype matrix for individuals from the target population.
        pos numpy.ndarray: Positions for the variants.
        method str: Method to create the physical distance matrix and genotype distance matrix.
                    Supported methods: vernot2016, vernot2014, and archie.
                    vernot2016 calculates genotype distances with a single individual from the target population.
                    vernot2014 calculates genotype distances with all individuals from the target population.
                    archie uses vernot2014 to calculate genotype distances but removes singletons before calculating genotype distances.
        match_bonus int: Bonus for matching genotypes of two different variants.
        max_mismatch int: Maximum genotype distance allowed.
        mismatch_penalty int: Penalty for mismatching genotypes of two different variants.

    Returns:
        sstar_scores list: The estimated sstar scores.
        sstar_snp_nums list: Numbers of sstar SNPs.
        haplotypes list: The haplotypes used for obtaining the estimated sstar scores.
    """
    def _create_matrixes(gt, pos, idx, method):
        hap = gt[:,idx]
        pos = pos[hap!=0]

        if method == 'vernot2016':
            # Calculate genotype distance with a single individual
            gt = hap[hap!=0]
            geno_matrix = np.tile(gt, (len(pos), 1))
            gd_matrix = np.transpose(geno_matrix)-geno_matrix
        elif method == 'vernot2014':
            # Calculate genotype distance with all individuals
            gd_matrix = distance_matrix(geno_matrix, geno_matrix, p=1)
        elif method == 'archie':
            geno_matrix = gt[hap!=0]
            # Remove singletons
            idx = np.sum(geno_matrix,axis=1)!=1
            pos = pos[idx]
            geno_matrix = geno_matrix[idx]
            gd_matrix = distance_matrix(geno_matrix, geno_matrix, p=1)
        else: raise ValueError(f'Method {method} is not supported!')

        pos_matrix = np.tile(pos, (len(pos), 1))
        pd_matrix = np.transpose(pos_matrix)-pos_matrix
        pd_matrix = pd_matrix.astype('float')

        return pd_matrix, gd_matrix, pos

    def _cal_ind_sstar(pd_matrix, gd_matrix, pos, match_bonus, max_mismatch, mismatch_penalty):
        pd_matrix[pd_matrix<10] = -np.inf
        pd_matrix[(pd_matrix>=10)*(gd_matrix==0)] += match_bonus
        pd_matrix[(pd_matrix>=10)*(gd_matrix>0)*(gd_matrix<=max_mismatch)] = mismatch_penalty
        pd_matrix[(pd_matrix>=10)*(gd_matrix>max_mismatch)] = -np.inf

        snp_num = len(pos)
        max_scores = [0] * snp_num
        max_score_snps = [[]] * snp_num
        for j in range(snp_num):
            max_score = -np.inf
            snps = []
            for i in range(j):
                score = max_scores[i] + pd_matrix[j,i]
                max_score_i = max(max_score, score, pd_matrix[j,i])
                if max_score_i == max_score: continue
                elif max_score_i == score: snps = max_score_snps[i] + [pos[j]]
                elif max_score_i == pd_matrix[j,i]: snps = [pos[i], pos[j]]
                max_score = max_score_i
            max_scores[j] = max_score
            max_score_snps[j] = snps

        try:
            sstar_score = max(max_scores)
            last_snp = max_scores.index(sstar_score)
            haplotype = max_score_snps[last_snp]
        except ValueError:
            sstar_score = 0
            last_snp = None
            haplotype = None

        return sstar_score, haplotype

    mut_num, ind_num = tgt_gt.shape
    sstar_scores = []
    sstar_snp_nums = []
    haplotypes = []
    for i in range(ind_num):
        pd_matrix, gd_matrix, sub_pos = _create_matrixes(tgt_gt, pos, i, method)
        sstar_score, haplotype = _cal_ind_sstar(pd_matrix, gd_matrix, sub_pos, match_bonus, max_mismatch, mismatch_penalty)
        if sstar_score == -np.inf: sstar_score = 0
        if (haplotype is not None) and (len(haplotype) != 0): 
            sstar_snp_num = len(haplotype)
            haplotype = ",".join([str(x) for x in haplotype])
        else: 
            sstar_snp_num = 'NA'
            haplotype = 'NA'
        sstar_scores.append(sstar_score)
        sstar_snp_nums.append(sstar_snp_num)
        haplotypes.append(haplotype)

    return sstar_scores, sstar_snp_nums, haplotypes


def cal_pr(ntruth_tracts, ninferred_tracts, ntrue_positives):
    """
    """
    precision = np.nan if float(ninferred_tracts) == 0 else ntrue_positives / float(ninferred_tracts) * 100
    recall = np.nan if float(ntruth_tracts) == 0 else ntrue_positives / float(ntruth_tracts) * 100

    return precision, recall
