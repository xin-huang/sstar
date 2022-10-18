import numpy as np
import scipy as sp
from scipy.spatial import distance_matrix


def cal_n_ton(tgt_gt):
    """
    Description:
        Calculates individual frequency spetra for haplotypes.

    Arguments:
        tgt_gt numpy.ndarray: Genotype matrix from the target population.

    Returns:
        spectra numpy.ndarray: Individual frequency spectra for haplotypes.
    """
    mut_num = tgt_gt.shape[0]
    ind_num = tgt_gt.shape[1]
    ploidy = tgt_gt.shape[2]
    hap_num = ind_num * ploidy

    tgt_gt = np.reshape(tgt_gt, (mut_num, hap_num))
    iv = np.ones((hap_num, 1))
    counts = tgt_gt * np.matmul(tgt_gt, iv)

    spectra = np.array([np.bincount(counts[:,idx].astype('int64'), minlength=hap_num+1) for idx in range(hap_num)])
    # ArchIE does not count non-segragating sites
    # spectra[:,0] = 0

    return spectra


def cal_ref_dist(ref_gt, tgt_gt):
    """
    Description:
        Calculates Euclidean distances between the reference and target populations.

    Arguments:
        ref_gt numpy.ndarray: Genotype matrix from the reference population.
        tgt_gt numpy.ndarray: Genotype matrix from the target population.

    Returns:
        dist_min numpy.ndarray: Minimum distances to the reference population.
    """
    ref_dist = distance_matrix(np.transpose(tgt_gt), np.transpose(ref_gt))
    dist_min = np.min(ref_dist, axis=1)

    return dist_min


def cal_tgt_dist(tgt_gt):
    """
    Description:
        Calculates Euclidean distances within the target population.

    Arguments:
        tgt_gt numpy.ndarray: Genotype matrix from the target population.

    Returns:
        tgt_dist numpy.ndarray: Pairwise distances between haplotypes from the target population.
        dist_mean numpy.ndarray: Mean of the pairwise distances between haplotypes from the target population.
        dist_var numpy.ndarray: Variation of the pairwise distances between haplotypes from the target population.
        dist_skew numpy.ndarray: Skew of the pairwise distances between haplotypes from the target population.
        dist_kurtosis numpy.ndarray: Kurtosis of the pairwise distances between haplotypes from the target population.
    """
    tgt_dist = distance_matrix(np.transpose(tgt_gt), np.transpose(tgt_gt))
    # ArchIE sorts tgt_dist
    # tgt_dist.sort()
    dist_mean = np.mean(tgt_dist, axis=1)
    dist_var = np.var(tgt_dist, axis=1)
    dist_skew = sp.stats.skew(tgt_dist, axis=1)
    dist_kurtosis = sp.stats.kurtosis(tgt_dist, axis=1)

    return tgt_dist, dist_mean, dist_var, dist_skew, dist_kurtosis


def cal_pvt_mut_num():
    pass


def cal_sstar():
    pass
