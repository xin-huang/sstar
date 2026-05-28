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
import math
import numpy as np
from typing import Optional, Tuple, Union


# @profile
def parse_ind_file(filename: str) -> list:
    """
    Read individual IDs from a text file.

    Parameters
    ----------
    filename : str
        Path to the file containing one individual ID per line.

    Returns
    -------
    list
        Individual IDs read from the file.
    """

    with open(filename, "r") as f:
        samples = [l.rstrip() for l in f.readlines()]

    if len(samples) == 0:
        raise Exception(f"No sample is found in {filename}! Please check your data.")

    return samples


# @profile
def read_geno_data(
    vcf: str,
    ind: list,
    anc_allele_file: Optional[str],
    filter_missing: bool,
) -> dict:
    """
    Read genotype data for selected individuals from a VCF file.

    Parameters
    ----------
    vcf : str
        Path to the VCF file containing genotype data.
    ind : list
        Individual IDs to read from the VCF file.
    anc_allele_file : str or None
        Path to the ancestral allele file. If None, ancestral allele filtering is not applied.
    filter_missing : bool
        If True, remove variants that are missing in all selected individuals.

    Returns
    -------
    dict
        Genotype data grouped by chromosome.
    """

    vcf = allel.read_vcf(vcf, alt_number=1, samples=ind)
    gt = vcf["calldata/GT"]
    chr_names = np.unique(vcf["variants/CHROM"])
    samples = vcf["samples"]
    pos = vcf["variants/POS"]
    ref = vcf["variants/REF"]
    alt = vcf["variants/ALT"]

    if anc_allele_file != None:
        anc_allele = read_anc_allele(anc_allele_file)
    data = dict()
    for c in chr_names:
        if c not in data.keys():
            data[c] = dict()
            data[c]["POS"] = pos
            data[c]["REF"] = ref
            data[c]["ALT"] = alt
            data[c]["GT"] = gt
        index = np.where(vcf["variants/CHROM"] == c)
        data = filter_data(data, c, index)
        # Remove missing data
        if filter_missing:
            index = data[c]["GT"].count_missing(axis=1) == len(samples)
            data = filter_data(data, c, ~index)
        if anc_allele_file != None:
            data = check_anc_allele(data, anc_allele, c)

    return data


# @profile
def filter_data(data: dict, c: str, index: np.ndarray) -> dict:
    """
    Filter genotype data for one chromosome by variant index.

    Parameters
    ----------
    data : dict
        Genotype data grouped by chromosome.
    c : str
        Chromosome name to filter.
    index : numpy.ndarray
        Boolean or integer index selecting variants to keep.

    Returns
    -------
    dict
        Filtered genotype data.
    """

    data[c]["POS"] = data[c]["POS"][index]
    data[c]["REF"] = data[c]["REF"][index]
    data[c]["ALT"] = data[c]["ALT"][index]
    data[c]["GT"] = allel.GenotypeArray(data[c]["GT"][index])

    return data


# @profile
def read_data(
    vcf_file: str,
    ref_ind_file: Optional[str],
    tgt_ind_file: Optional[str],
    src_ind_file: Optional[str],
    anc_allele_file: Optional[str],
) -> Tuple[
    Optional[dict],
    Optional[list],
    Optional[dict],
    Optional[list],
    Optional[dict],
    Optional[list],
]:
    """
    Read and harmonize reference, target, and source genotype data.

    Parameters
    ----------
    vcf_file : str
        Path to the VCF file containing genotype data.
    ref_ind_file : str or None
        Path to the file containing reference individual IDs. If None, reference data is not loaded.
    tgt_ind_file : str or None
        Path to the file containing target individual IDs. If None, target data is not loaded.
    src_ind_file : str or None
        Path to the file containing source individual IDs. If None, source data is not loaded.
    anc_allele_file : str or None
        Path to the ancestral allele file. If None, ancestral allele filtering is not applied.

    Returns
    -------
    tuple
        Tuple containing reference data, reference samples, target data, target samples, source data, and source samples.
    """
    # --- Initialize group structure and containers ----------------------------
    # Each group defines its individual file and whether data is read as phased.
    groups = {
        "ref": (ref_ind_file, True),
        "tgt": (tgt_ind_file, True),
        "src": (src_ind_file, True),
    }
    data, samples = {"ref": None, "tgt": None, "src": None}, {
        "ref": None,
        "tgt": None,
        "src": None,
    }

    # --- Load genotype data for each population group -------------------------
    # Parse sample lists and load genotype data from VCF if input file provided.
    for k, (f, filter_missing) in groups.items():
        if f is not None:
            s = parse_ind_file(f)
            samples[k] = s
            # Keep src less stringently filtered for missing (as in sstar)
            data[k] = read_geno_data(vcf_file, s, anc_allele_file, filter_missing)

    # --- Remove variants fixed as hom-ALT in BOTH ref and tgt -----------------
    # For each shared chromosome, remove variants fixed in both populations.
    # The same sites are removed from the source population for consistency.
    if data["ref"] and data["tgt"]:
        for c in data["ref"].keys() & data["tgt"].keys():
            ref_fixed = np.sum(data["ref"][c]["GT"].is_hom_alt(), axis=1) == len(
                samples["ref"]
            )
            tgt_fixed = np.sum(data["tgt"][c]["GT"].is_hom_alt(), axis=1) == len(
                samples["tgt"]
            )
            keep = ~(ref_fixed & tgt_fixed)
            fixed_pos = data["ref"][c]["POS"][~keep]
            data["ref"] = filter_data(data["ref"], c, keep)
            data["tgt"] = filter_data(data["tgt"], c, keep)
            if data["src"] and (c in data["src"]):
                src_keep = ~np.in1d(data["src"][c]["POS"], fixed_pos)
                data["src"] = filter_data(data["src"], c, src_keep)

    # --- Return loaded and processed data -------------------------------------
    return (
        data["ref"],
        samples["ref"],
        data["tgt"],
        samples["tgt"],
        data["src"],
        samples["src"],
    )


# @profile
def get_ref_alt_allele(
    ref: Union[list, np.ndarray],
    alt: Union[list, np.ndarray],
    pos: Union[list, np.ndarray],
) -> Tuple[dict, dict]:
    """
    Index REF and ALT alleles by genomic position.

    Parameters
    ----------
    ref : list or numpy.ndarray
        REF alleles.
    alt : list or numpy.ndarray
        ALT alleles.
    pos : list or numpy.ndarray
        Genomic positions corresponding to `ref` and `alt`.

    Returns
    -------
    tuple
        Tuple containing REF alleles and ALT alleles keyed by position.
    """

    ref_allele = dict()
    alt_allele = dict()

    for i in range(len(pos)):
        r = ref[i]
        a = alt[i]
        p = pos[i]
        ref_allele[p] = r
        alt_allele[p] = a

    return ref_allele, alt_allele


# @profile
def read_anc_allele(anc_allele_file: str) -> dict:
    """
    Read ancestral allele information from a BED-like file.

    Parameters
    ----------
    anc_allele_file : str
        Path to the file containing ancestral allele information.

    Returns
    -------
    dict
        Ancestral alleles keyed by chromosome and position.
    """

    anc_allele = dict()
    with open(anc_allele_file, "r") as f:
        for line in f.readlines():
            e = line.rstrip().split()
            if e[0] not in anc_allele:
                anc_allele[e[0]] = dict()
            anc_allele[e[0]][int(e[2])] = e[3]

    if not anc_allele:
        raise Exception(f"No ancestral allele is found! Please check your data.")

    return anc_allele


# @profile
def check_anc_allele(data: dict, anc_allele: dict, c: str) -> dict:
    """
    Filter and polarize genotype data using ancestral allele information.

    Variants are retained only when the ancestral allele matches either the REF or ALT allele. Variants whose ALT allele is ancestral are flipped.

    Parameters
    ----------
    data : dict
        Genotype data grouped by chromosome.
    anc_allele : dict
        Ancestral alleles keyed by chromosome and position.
    c : str
        Chromosome name to process.

    Returns
    -------
    dict
        Filtered and polarized genotype data.
    """

    ref_allele, alt_allele = get_ref_alt_allele(
        data[c]["REF"], data[c]["ALT"], data[c]["POS"]
    )
    # Remove variants not in the ancestral allele file
    intersect_snps = np.intersect1d(list(ref_allele.keys()), list(anc_allele[c].keys()))
    # Remove variants that neither the ref allele nor the alt allele is the ancestral allele
    removed_snps = []
    # Flip variants that the alt allele is the ancestral allele
    flipped_snps = []

    for v in intersect_snps:
        if (anc_allele[c][v] != ref_allele[v]) and (anc_allele[c][v] != alt_allele[v]):
            removed_snps.append(v)
        elif anc_allele[c][v] == alt_allele[v]:
            flipped_snps.append(v)

    intersect_snps = np.in1d(data[c]["POS"], intersect_snps)
    data = filter_data(data, c, intersect_snps)

    if len(removed_snps) != 0:
        remained_snps = np.logical_not(np.in1d(data[c]["POS"], removed_snps))
        data = filter_data(data, c, remained_snps)

    is_flipped_snps = np.in1d(data[c]["POS"], flipped_snps)
    # Assume no missing data
    for i in range(len(data[c]["POS"])):
        if is_flipped_snps[i]:
            data[c]["GT"][i] = allel.GenotypeVector(abs(data[c]["GT"][i] - 1))

    return data


# @profile
def read_mapped_region_file(mapped_region: Optional[str]) -> Optional[dict]:
    """
    Read mapped genomic regions from a BED file.

    Parameters
    ----------
    mapped_region : str or None
        Path to the BED file containing mapped regions. If None, no mapped-region intervals are loaded.

    Returns
    -------
    dict or None
        Mapped-region intervals by chromosome, or None if `mapped_region` is None.
    """

    if mapped_region != None:
        mapped_intervals = dict()
        with open(mapped_region, "r") as m:
            for line in m.readlines():
                elements = line.rstrip().split("\t")
                if elements[0] not in mapped_intervals.keys():
                    mapped_intervals[elements[0]] = []
                mapped_intervals[elements[0]].append(
                    (int(elements[1]), int(elements[2]))
                )
    else:
        mapped_intervals = None

    return mapped_intervals


# @profile
def cal_matchpct(
    chr_name: str,
    mapped_intervals: Optional[dict],
    data: dict,
    src_data: dict,
    tgt_ind_index: int,
    src_ind_index: int,
    hap_index: int,
    win_start: Union[int, str],
    win_end: Union[int, str],
    sample_size: int,
) -> list:
    """
    Calculate haplotype source-match statistics for one target-source pair.

    Parameters
    ----------
    chr_name : str
        Chromosome name.
    mapped_intervals : dict or None
        Mapped-region intervals by chromosome. If None, full windows are treated as mapped.
    data : dict
        Target genotype data by chromosome.
    src_data : dict
        Source genotype data by chromosome.
    tgt_ind_index : int
        Index of the target individual.
    src_ind_index : int
        Index of the source individual.
    hap_index : int
        Index of the target haplotype.
    win_start : int or str
        Start position of the window, or `"NA"` when no window is available.
    win_end : int or str
        End position of the window, or `"NA"` when no window is available.
    sample_size : int
        Number of target individuals analyzed.

    Returns
    -------
    list
        Haplotype statistics for the target-source pair.
    """

    if (win_start == "NA") and (win_end == "NA"):
        return ["NA", "NA", "NA", "NA", "NA", "NA", "NA"]

    win_start = int(win_start)
    win_end = int(win_end)

    hap_match_pct = "NA"
    gt = data[chr_name]["GT"]
    pos = data[chr_name]["POS"]
    src_gt = src_data[chr_name]["GT"]
    src_pos = src_data[chr_name]["POS"]
    res = []

    sub_snps = np.where((pos >= win_start) & (pos <= win_end))[0]
    sub_gt = gt[sub_snps][:, tgt_ind_index]
    sub_pos = pos[sub_snps]

    sub_src_snps = np.where((src_pos >= win_start) & (src_pos <= win_end))[0]
    sub_src_gt = src_gt[sub_src_snps][:, src_ind_index]
    sub_src_pos = src_pos[sub_src_snps]
    src_hom_variants = sub_src_pos[sub_src_gt.is_hom_alt()]
    src_het_variants = sub_src_pos[sub_src_gt.is_het()]
    src_variants = np.concatenate((src_hom_variants, src_het_variants))

    hap = sub_gt[:, hap_index]
    hap_len = win_end - win_start
    hap_mapped_len = _cal_mapped_len(mapped_intervals, chr_name, win_start, win_end)
    hap_variants_num, hap_site_num, hap_match_src_allele_num, hap_sfs, hap_match_pct = (
        _cal_hap_stats(
            gt[sub_snps],
            hap,
            sub_pos,
            src_variants,
            src_hom_variants,
            src_het_variants,
            sample_size,
        )
    )

    res.append(hap_variants_num)
    res.append(hap_len)
    res.append(hap_mapped_len)
    res.append(hap_match_src_allele_num)
    res.append(hap_site_num)
    res.append(hap_sfs)
    res.append(hap_match_pct)

    return res


# @profile
def _cal_mapped_len(
    mapped_intervals: Optional[dict], chr_name: str, win_start: int, win_end: int
) -> int:
    """
    Calculate mapped length within a genomic window.

    Parameters
    ----------
    mapped_intervals : dict or None
        Mapped-region intervals by chromosome. If None, the full window length is used.
    chr_name : str
        Chromosome name.
    win_start : int
        Start position of the window.
    win_end : int
        End position of the window.

    Returns
    -------
    int
        Mapped length rounded down to the nearest kilobase.
    """

    if mapped_intervals is None or chr_name not in mapped_intervals.keys():
        mapped_len = win_end - win_start
    else:
        overlaps = [
            (idx[0], idx[1])
            for idx in mapped_intervals[chr_name]
            if (win_start >= idx[0] and win_start <= idx[1])
            or (win_end >= idx[0] and win_end <= idx[1])
            or (win_start <= idx[0] and win_end >= idx[1])
        ]
        mapped_len = 0
        for idx in overlaps:
            if (idx[0] <= win_start) and (idx[1] >= win_end):
                mapped_len += win_end - win_start
            elif (idx[0] >= win_start) and (idx[1] <= win_end):
                mapped_len += idx[1] - idx[0]
            elif (idx[0] <= win_start) and (idx[1] <= win_end):
                mapped_len += idx[1] - win_start
            else:
                mapped_len += win_end - idx[0]

    mapped_len = mapped_len // 1000 * 1000

    return mapped_len


# @profile
def _cal_hap_stats(
    gt: allel.GenotypeArray,
    hap: Optional[allel.GenotypeVector],
    pos: Union[list, np.ndarray],
    src_variants: Union[list, np.ndarray],
    src_hom_variants: Union[list, np.ndarray],
    src_het_variants: Union[list, np.ndarray],
    sample_size: int,
) -> tuple:
    """
    Calculate summary statistics for one haplotype in a window.

    Parameters
    ----------
    gt : allel.GenotypeArray
        Genotype data for all haplotypes in the window.
    hap : allel.GenotypeVector or None
        Haplotype to summarize. If None, all returned values are `"NA"`.
    pos : list or numpy.ndarray
        Variant positions in the window.
    src_variants : list or numpy.ndarray
        Variant positions observed in the source individual.
    src_hom_variants : list or numpy.ndarray
        Homozygous ALT variant positions observed in the source individual.
    src_het_variants : list or numpy.ndarray
        Heterozygous variant positions observed in the source individual.
    sample_size : int
        Number of target individuals analyzed.

    Returns
    -------
    tuple
        Haplotype variant count, site count, source-matching allele count, SFS summary, and match percent.
    """

    if hap is None:
        return "NA", "NA", "NA", "NA", "NA"
    else:
        hap_variants = pos[np.equal(hap, 1)]
        hap_variants_num = len(hap_variants)
        # Assume the alternative allele is the derived allele
        hap_shared_src_hom_site_num = len(
            np.intersect1d(hap_variants, src_hom_variants)
        )
        hap_shared_src_het_site_num = len(
            np.intersect1d(hap_variants, src_het_variants)
        )
        hap_site_num = len(np.union1d(hap_variants, src_variants))
        hap_match_src_allele_num = (
            hap_shared_src_hom_site_num + 0.5 * hap_shared_src_het_site_num
        )
        hap_shared_src_site_num = (
            hap_shared_src_hom_site_num + hap_shared_src_het_site_num
        )
        if hap_site_num != 0:
            hap_match_pct = round(hap_match_src_allele_num / hap_site_num, 6)
        else:
            hap_match_pct = "NA"

        hap_sfs = np.sum(np.sum(gt[hap == 1], axis=2), axis=1)
        if hap_sfs.size != 0:
            hap_sfs_mean = np.mean(hap_sfs)
            # See https://stackoverflow.com/questions/10825926/python-3-x-rounding-behavior
            # if not np.isnan(sfs_mean): sfs_mean = int(round(sfs_mean))
            # if not np.isnan(hap_sfs_mean): hap_sfs = int(int(py2round(hap_sfs_mean))/10*108)
            # if not np.isnan(hap_sfs_mean): hap_sfs = int(py2round(hap_sfs_mean))/(2*sample_size)
            if not np.isnan(hap_sfs_mean):
                hap_sfs = round(hap_sfs_mean / (2 * sample_size), 6)
        else:
            hap_sfs = np.nan

    return (
        hap_variants_num,
        hap_site_num,
        hap_match_src_allele_num,
        hap_sfs,
        hap_match_pct,
    )


def py2round(x: float, d: int = 0) -> float:
    """
    Round a number using Python 2 half-away-from-zero behavior.

    Parameters
    ----------
    x : float
        Value to round.
    d : int, optional
        Number of decimal places. Default: `0`.

    Returns
    -------
    float
        Rounded value.
    """
    p = 10**d
    if x > 0:
        return float(math.floor((x * p) + 0.5)) / p
    elif x < 0:
        return float(math.ceil((x * p) - 0.5)) / p
    else:
        return 0.0
