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

import allel
import numpy as np
import os
from multiprocessing import Process, Queue
from typing import Optional, Union
from sstar.utils import read_data


# @profile
def cal_s_star(
    vcf: str,
    ref_ind_file: str,
    tgt_ind_file: str,
    anc_allele_file: Optional[str],
    output: str,
    win_len: int,
    win_step: int,
    thread: int,
    match_bonus: int,
    max_mismatch: int,
    mismatch_penalty: int,
    is_phased: bool = False,
) -> None:
    """
    Calculate S* scores from genotype data.

    Parameters
    ----------
    vcf : str
        Path to the VCF file containing genotype data.
    ref_ind_file : str
        Path to the file containing reference individual IDs.
    tgt_ind_file : str
        Path to the file containing target individual IDs.
    anc_allele_file : str or None
        Path to the ancestral allele file. If None, ancestral allele information is not used.
    output : str
        Path to the output score file.
    win_len : int
        Length of sliding windows.
    win_step : int
        Step size of sliding windows.
    thread : int
        Number of worker processes.
    match_bonus : int
        Bonus for matching genotypes between two variants.
    max_mismatch : int
        Maximum genotype distance allowed before a pair is discarded.
    mismatch_penalty : int
        Penalty for mismatching genotypes between two variants.
    is_phased : bool, optional
        If True, calculate scores on phased haplotypes; otherwise calculate scores on dosages. Default: `False`.
    """

    # --- Load data in 3D GenotypeArray form (V x N x ploidy) -----------------
    ref_data, ref_samples, tgt_data, tgt_samples, src_data, src_samples = read_data(
        vcf_file=vcf,
        ref_ind_file=ref_ind_file,
        tgt_ind_file=tgt_ind_file,
        src_ind_file=None,
        anc_allele_file=anc_allele_file,
        is_phased=is_phased,
    )

    # --- Convert GT to 2D: phased haplotypes or genotype dosage -------------
    for data_dict in (ref_data, tgt_data, src_data):
        if data_dict is None:
            continue
        for d in data_dict.values():
            GT = d["GT"]  # allel.GenotypeArray (V, N, P)
            A = GT.values  # ndarray (V, N, P)
            if is_phased:
                v, n, p = A.shape
                d["GT"] = A.reshape(v, n * p)  # (V, N*P) haplotypes (0/1)
            else:
                d["GT"] = A.sum(axis=2)  # (V, N) dosage (0/1/2)
            # ensure POS is a numpy array
            if not isinstance(d["POS"], np.ndarray):
                d["POS"] = np.asarray(d["POS"])

    # --- Remove variants observed in the reference populations --------------
    # original logic (Xin Huang):
    #   # Remove variants observed in the reference populations
    #   # Assume 1 is the alt allele
    #   # (original used ref_data[c]["GT"].is_hom_ref(); in 2D dosage, hom-ref == 0)
    for c, d_tgt in tgt_data.items():
        ref_gt = ref_data[c]["GT"]  # 2D: (V_ref, N_ref)
        variants_not_in_ref = np.all(ref_gt == 0, axis=1)
        for key in ("POS", "REF", "ALT", "GT"):
            if key in d_tgt:
                d_tgt[key] = d_tgt[key][variants_not_in_ref]

    if thread > 1:
        thread = min(os.cpu_count() - 1, len(tgt_samples), thread)

    _cal_score(
        ref_data,
        tgt_data,
        tgt_samples,
        win_len,
        win_step,
        output,
        thread,
        match_bonus,
        max_mismatch,
        mismatch_penalty,
        is_phased=is_phased,
    )


# @profile
def _cal_score(
    ref_data: dict,
    tgt_data: dict,
    samples: list,
    win_len: int,
    win_step: int,
    output: str,
    thread: int,
    match_bonus: int,
    max_mismatch: int,
    mismatch_penalty: int,
    is_phased: bool = False,
) -> None:
    """
    Manage multiprocessing for S* score calculation.

    Parameters
    ----------
    ref_data : dict
        Reference genotype data by chromosome.
    tgt_data : dict
        Target genotype data by chromosome.
    samples : list
        Target sample IDs.
    win_len : int
        Length of sliding windows.
    win_step : int
        Step size of sliding windows.
    output : str
        Path to the output score file.
    thread : int
        Number of worker processes.
    match_bonus : int
        Bonus for matching genotypes between two variants.
    max_mismatch : int
        Maximum genotype distance allowed before a pair is discarded.
    mismatch_penalty : int
        Penalty for mismatching genotypes between two variants.
    is_phased : bool, optional
        If True, process each haplotype separately; otherwise process each individual as dosages. Default: `False`.
    """

    try:
        from pytest_cov.embed import cleanup_on_sigterm
    except ImportError:
        pass
    else:
        cleanup_on_sigterm()

    res = []
    # Use multiprocessing to calculate S* scores in different samples
    in_queue, out_queue = Queue(), Queue()
    workers = [
        Process(
            target=_cal_score_worker,
            args=(
                in_queue,
                out_queue,
                ref_data,
                tgt_data,
                win_len,
                win_step,
                match_bonus,
                max_mismatch,
                mismatch_penalty,
                is_phased,
            ),
        )
        for _ in range(thread)
    ]

    if is_phased:
        # GT shape is V x (N * ploidy) after flattening, so we need 1 job per haplotype
        any_chr = next(iter(tgt_data.values()))
        n_cols = any_chr["GT"].shape[1]  # total haplotypes
        n_samples = len(samples)
        ploidy = n_cols // n_samples

        n_jobs = n_cols
        for hap_idx in range(n_cols):
            sample_idx = hap_idx // ploidy
            hap_no = (hap_idx % ploidy) + 1
            sample_name = samples[sample_idx]
            in_queue.put((hap_idx, sample_name, str(hap_no)))
    else:
        n_jobs = len(samples)
        # CHANGE REQUESTED BY XIN (Dec 5, 2025): use n_jobs in the loop range
        for s in range(n_jobs):
            sample_name = samples[s]
            in_queue.put((s, sample_name, "NA"))

    try:
        for worker in workers:
            worker.start()

        for _ in range(n_jobs):
            item = out_queue.get()
            if item != "":
                res.append(item)

        for worker in workers:
            worker.terminate()
    finally:
        for worker in workers:
            worker.join()

    header = "chrom\tstart\tend\tsample\thap_index\tS*_score\tregion_ind_SNP_number\tS*_SNP_number\tS*_SNPs"
    with open(output, "w") as o:
        o.write(header + "\n")
        o.write("\n".join(res))
        o.write("\n")


# @profile
def _cal_score_worker(
    in_queue: Queue,
    out_queue: Queue,
    ref_data: dict,
    tgt_data: dict,
    win_len: int,
    win_step: int,
    match_bonus: int,
    max_mismatch: int,
    mismatch_penalty: int,
    is_phased: bool = False,
) -> None:
    """
    Worker process for S* score calculation.

    Parameters
    ----------
    in_queue : multiprocessing.Queue
        Queue providing target sample or haplotype jobs.
    out_queue : multiprocessing.Queue
        Queue receiving formatted S* output.
    ref_data : dict
        Reference genotype data by chromosome.
    tgt_data : dict
        Target genotype data by chromosome.
    win_len : int
        Length of sliding windows.
    win_step : int
        Step size of sliding windows.
    match_bonus : int
        Bonus for matching genotypes between two variants.
    max_mismatch : int
        Maximum genotype distance allowed before a pair is discarded.
    mismatch_penalty : int
        Penalty for mismatching genotypes between two variants.
    is_phased : bool, optional
        If True, process each haplotype separately; otherwise process each individual as dosages. Default: `False`.
    """

    while True:
        s, sample_name, hap_index = in_queue.get()
        chr_names = tgt_data.keys()
        for c in chr_names:
            tgt_gt = tgt_data[c]["GT"]  # 2D: (V, N_tgt)
            tgt_pos = tgt_data[c]["POS"]
            ind = tgt_gt[:, s]  # 1D: (V,)

            ref_gt = ref_data[c]["GT"]  # 2D: (V, N_ref)
            ref_pos = ref_data[c]["POS"]
            # original: ref_sub_pos = ref_pos[~np.all(ref_gt.is_hom_ref(), axis=1)]
            # 2D dosage/haplotype equivalent: hom-ref == 0
            ref_sub_pos = ref_pos[~np.all(ref_gt == 0, axis=1)]
            # Assume the ref allele is 0 and the alt allele is 1
            # original: tgt_sub_gt = tgt_gt[~ind.is_hom_ref()][:,s]; tgt_sub_pos = tgt_pos[~ind.is_hom_ref()]
            # 2D equivalent: keep sites where this target sample is not hom-ref (ind != 0)
            mask = ind != 0
            tgt_sub_gt = ind[mask]
            tgt_sub_pos = tgt_pos[mask]
            res = _cal_score_ind(
                c,
                sample_name,
                hap_index,
                ref_sub_pos,
                tgt_sub_pos,
                tgt_sub_gt,
                win_step,
                win_len,
                match_bonus,
                max_mismatch,
                mismatch_penalty,
                chr_last_pos=int(tgt_pos[-1]),
            )
            out_queue.put("\n".join(res))


# @profile
def _cal_score_ind(
    chr_name: str,
    sample_name: str,
    hap_index: str,
    ref_pos: Union[list, np.ndarray],
    tgt_pos: Union[list, np.ndarray],
    tgt_gt: np.ndarray,
    win_step: int,
    win_len: int,
    match_bonus: int,
    max_mismatch: int,
    mismatch_penalty: int,
    chr_last_pos: int,
) -> list:
    """
    Calculate S* scores for one target individual or haplotype.

    Parameters
    ----------
    chr_name : str
        Chromosome name.
    sample_name : str
        Sample ID written to the score output.
    hap_index : str
        Haplotype index, or `"NA"` for unphased runs.
    ref_pos : list or numpy.ndarray
        Variant positions observed in the reference population.
    tgt_pos : list or numpy.ndarray
        Variant positions observed in the target individual or haplotype.
    tgt_gt : numpy.ndarray
        One-dimensional target dosages or haploid allele values.
    win_step : int
        Step size of sliding windows.
    win_len : int
        Length of sliding windows.
    match_bonus : int
        Bonus for matching genotypes between two variants.
    max_mismatch : int
        Maximum genotype distance allowed before a pair is discarded.
    mismatch_penalty : int
        Penalty for mismatching genotypes between two variants.

    Returns
    -------
    list
        Formatted output lines containing S* scores for each window.
    """
    res = []
    if len(tgt_gt) <= 2:
        # An individual should have at least three SNPs for S* calculation
        last_pos = 0
        win_start = 0
    else:
        win_start = (tgt_pos[0] + win_step) // win_step * win_step - win_len
        if win_start < 0:
            win_start = 0
        last_pos = chr_last_pos

    while last_pos > win_start:
        win_end = win_start + win_len
        tgt_snps = np.where((tgt_pos > win_start) & (tgt_pos <= win_end))[0]
        ref_snps = np.where((ref_pos > win_start) & (ref_pos <= win_end))[0]
        total_snps = len(np.union1d(ref_pos[ref_snps], tgt_pos[tgt_snps]))

        # Only calculate S* for a region with more than two SNPs
        tgt_snp_num = len(tgt_snps)
        if tgt_snp_num > 2:
            tgt_snps_gt = tgt_gt[tgt_snps]
            tgt_snps_pos = tgt_pos[tgt_snps]
            # Store the maximum score ended at the j-th SNP
            max_scores = [0] * tgt_snp_num
            max_score_snps = [[]] * tgt_snp_num
            for j in range(tgt_snp_num):
                max_score = -np.inf
                snps = []
                for i in range(j):
                    geno_dist = abs(tgt_snps_gt[j] - tgt_snps_gt[i])

                    phy_dist = abs(tgt_snps_pos[j] - tgt_snps_pos[i])
                    if phy_dist < 10:
                        score_ij = -np.inf
                    elif geno_dist == 0:
                        score_ij = match_bonus + phy_dist
                    elif geno_dist <= max_mismatch:
                        score_ij = mismatch_penalty
                    else:
                        score_ij = -np.inf
                    score = max_scores[i] + score_ij
                    max_score_i = max(max_score, score, score_ij)
                    if max_score_i == max_score:
                        continue
                    elif max_score_i == score:
                        snps = max_score_snps[i] + [str(tgt_snps_pos[j])]
                    elif max_score_i == score_ij:
                        snps = [str(tgt_snps_pos[i]), str(tgt_snps_pos[j])]
                    max_score = max_score_i
                max_scores[j] = max_score
                max_score_snps[j] = snps

            s_star_score = max(max_scores)
            last_snp = max_scores.index(s_star_score)
            haplotype = max_score_snps[last_snp]

            s_star_snp_num = len(haplotype)
            if len(haplotype) > 2:
                haplotype = ",".join(haplotype)
            else:
                s_star_snp_num = "NA"
                s_star_score = "NA"
                haplotype = "NA"
        else:
            s_star_snp_num = "NA"
            s_star_score = "NA"
            haplotype = "NA"

        line = (
            f"{chr_name}\t{win_start}\t{win_end}\t{sample_name}\t"
            f"{hap_index}\t{s_star_score}\t{total_snps}\t{s_star_snp_num}\t{haplotype}"
        )

        win_start += win_step
        res.append(line)

    return res
