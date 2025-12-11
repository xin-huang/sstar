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
from sstar.utils import read_data


# @profile
def cal_s_star(
    vcf,
    ref_ind_file,
    tgt_ind_file,
    anc_allele_file,
    output,
    win_len,
    win_step,
    thread,
    match_bonus,
    max_mismatch,
    mismatch_penalty,
    is_phased: bool = False,  # added but currently unused; kept for future phased support
):
    """
    Description:
        Calculate S* scores.

    Arguments:
        vcf str: Name of the VCF file containing genotype data.
        ref_ind_file str: Name of the file containing sample information from the reference population.
        tgt_ind_file str: Name of the file containing sample information from the target population.
        anc_allele_file str: Name of the file containing ancestral allele information.
        output str: Name of the output file.
        win_len int: Length of sliding windows.
        win_step int: Step size of sliding windows.
        thread int: Number of threads.
        match_bonus int: Bonus for matching genotypes of two different variants.
        max_mismatch int: Maximum genotype distance allowed.
        mismatch_penalty int: Penalty for mismatching genotypes of two different variants.
        is_phased bool: If True, flatten haplotypes; else convert to dosages.
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
        is_phased=is_phased,  # threaded but currently unused; kept for future phased support
    )


# @profile
def _cal_score(
    ref_data,
    tgt_data,
    samples,
    win_len,
    win_step,
    output,
    thread,
    match_bonus,
    max_mismatch,
    mismatch_penalty,
    is_phased: bool = False,  # threaded but currently unused; kept for future phased support
):
    """
    Description:
        Helper function to manage worker functions to calculate S* scores with multiprocessing.

    Arguments:
        ref_data dict: Dictionary containing genotype data from the reference population.
        tgt_data dict: Dictionary containing genotype data from the target population.
        samples list: Sample information.
        win_len int: Length of sliding windows.
        win_step int: Step size of sliding windows.
        output str: Name of the output file.
        thread int: Number of threads.
        match_bonus int: Bonus for matching genotypes of two different variants.
        max_mismatch int: Maximum genotype distance allowed.
        mismatch_penalty int: Penalty for mismatching genotypes of two different variants.
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
                is_phased,  # threaded but currently unused
            ),
        )
        for _ in range(thread)
    ]
    if is_phased:
    # GT shape is V x (N * ploidy) after flattening, so we need 1 job per haplotype
        any_chr = next(iter(tgt_data.values()))
        n_cols = any_chr["GT"].shape[1]      # total haplotypes
        n_samples = len(samples)
        ploidy = n_cols // n_samples

        n_jobs = n_cols
        for hap_idx in range(n_cols):
            sample_idx = hap_idx // ploidy
            hap_no = (hap_idx % ploidy) + 1
            sample_name = f"{samples[sample_idx]}_hap{hap_no}"
            in_queue.put((hap_idx, sample_name))
    else:
        n_jobs = len(samples)
        for s in range(len(samples)):
            in_queue.put((s, samples[s]))

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

    header = "chrom\tstart\tend\tsample\tS*_score\tregion_ind_SNP_number\tS*_SNP_number\tS*_SNPs"
    with open(output, "w") as o:
        o.write(header + "\n")
        o.write("\n".join(res))
        o.write("\n")


# @profile
def _cal_score_worker(
    in_queue,
    out_queue,
    ref_data,
    tgt_data,
    win_len,
    win_step,
    match_bonus,
    max_mismatch,
    mismatch_penalty,
    is_phased: bool = False,  # threaded but currently unused; kept for future phased support
):
    """
    Description:
        Worker function to calculate S* scores with multiprocessing.

    Arguments:
        in_queue multiprocessing.Queue: multiprocessing.Queue instance to receive parameters from the manager.
        out_queue multiprocessing.Queue: multiprocessing.Queue instance to send results back to the manager.
        ref_data dict: Dictionary containing genotype data from the reference population.
        tgt_data dict: Dictionary containing genotype data from the target population.
        win_len int: Length of sliding windows.
        win_step int: Step size of sliding windows.
        match_bonus int: Bonus for matching genotypes of two different variants.
        max_mismatch int: Maximum genotype distance allowed.
        mismatch_penalty int: Penalty for mismatching genotypes of two different variants.
    """

    while True:
        s, sample_name = in_queue.get()
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
                ref_sub_pos,
                tgt_sub_pos,
                tgt_sub_gt,
                win_step,
                win_len,
                match_bonus,
                max_mismatch,
                mismatch_penalty,
            )
            out_queue.put("\n".join(res))


# @profile
def _cal_score_ind(
    chr_name,
    sample_name,
    ref_pos,
    tgt_pos,
    tgt_gt,
    win_step,
    win_len,
    match_bonus,
    max_mismatch,
    mismatch_penalty,
):
    """
    Description:
        Helper function for calculating S* in a single individual.

    Arguments:
        chr_name str: Name of the chromosome.
        sample_name str: Name of the sample.
        ref_pos list: List of position for variants in reference populations.
        tgt_pos list: List of position for variants in target populations.
        tgt_gt: Genotype data (here: 1D array of dosages/haploid alleles) from target populations.
        win_step int: Length of sliding windows.
        win_len int: Step size of sliding windows.
        match_bonus int: Bonus for matching genotypes of two different variants.
        max_mismatch int: Maximum genotype distance allowed.
        mismatch_penalty int: Penalty for mismatching genotypes of two different variants.

    Returns:
        res list: Results of calculated S*
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
        last_pos = tgt_pos[-1]

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
            f"{s_star_score}\t{total_snps}\t{s_star_snp_num}\t{haplotype}"
        )

        win_start += win_step
        res.append(line)

    return res
