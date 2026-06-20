# Copyright 2026 Xin Huang and Andrea Koca
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

import os
from typing import Optional, Union

import allel
import numpy as np
import pandas as pd

from sstar.utils import parse_ind_file, read_geno_data


def calc_match_pct(
    x: Union[float, list, np.ndarray],
    y: Union[float, list, np.ndarray],
    denominator: int,
    P: int = 2,
) -> Union[float, str]:
    """
    Calculate window-normalized dosage match rate between two genotype arrays.

    Missing values must be encoded as -1 and contribute zero to the numerator.
    The denominator is the population-level number of variants in the window.
    Callable sites are positions where both genotype arrays have dosage values
    greater than or equal to zero.

    Parameters
    ----------
    x : float, list, or numpy.ndarray
        ALT-allele dosage values for the first genotype vector.
    y : float, list, or numpy.ndarray
        ALT-allele dosage values for the second genotype vector.
    denominator : int
        Number of population-level variants in the window. Missing genotypes
        are included in this count even though they do not contribute to the
        numerator.
    P : int, optional
        Ploidy used to normalize dosage differences. Default: 2.

    Returns
    -------
    float or str
        Window-normalized match rate in [0, 1], or "NA" if the window has no
        variants. If the window has variants but no callable pairwise sites,
        return 0.0.
    """
    if denominator == 0:
        return "NA"

    x = np.atleast_1d(np.asarray(x, dtype=float))
    y = np.atleast_1d(np.asarray(y, dtype=float))

    ok = (x >= 0) & (y >= 0)
    score_sum = np.sum(1.0 - np.abs(x[ok] - y[ok]) / float(P))

    return float(score_sum / denominator)


def _build_position_index(pos: Union[list, np.ndarray]) -> dict[int, int]:
    """
    Build a genomic-position to array-index lookup.

    Parameters
    ----------
    pos : list or numpy.ndarray
        Variant positions.

    Returns
    -------
    dict[int, int]
        Mapping from genomic position to genotype-array row index.
    """
    return {int(p): i for i, p in enumerate(pos)}


def _dosage_for_position_index(
    gt: allel.GenotypeArray,
    pos_to_i: dict[int, int],
    ind_index: int,
    positions: list,
    hap_index: Optional[int] = None,
    P: int = 2,
) -> np.ndarray:
    """
    Return ALT dosage for one individual using a precomputed position index.

    Missing genotypes and absent positions are encoded as -1. When extracting a
    single haplotype, called alleles are multiplied by ``P`` so haplotype allele
    states are represented on the same dosage scale used for diploid genotypes.

    Parameters
    ----------
    gt : allel.GenotypeArray
        Genotype array with variant, sample, and ploidy axes.
    pos_to_i : dict[int, int]
        Mapping from genomic position to genotype-array row index.
    ind_index : int
        Index of the individual to extract.
    positions : list
        Genomic positions to query.
    hap_index : int or None, optional
        Zero-based haplotype index to extract for phased targets. When None,
        return diploid ALT dosage using scikit-allel's to_n_alt().
    P : int, optional
        Ploidy multiplier for haploid allele extraction. Default: 2, which
        maps haplotype allele 0 to dosage 0 and haplotype allele 1 to dosage 2.

    Returns
    -------
    numpy.ndarray
        Dosage array aligned to positions. Missing genotypes and absent
        positions are encoded as -1.
    """
    gt_arr = np.asarray(gt)
    out = np.full(len(positions), -1.0, dtype=float)
    idx = []
    out_i = []

    for j, p in enumerate(positions):
        i = pos_to_i.get(int(p))
        if i is not None:
            idx.append(i)
            out_i.append(j)

    if not idx:
        return out

    idx = np.asarray(idx, dtype=int)
    out_i = np.asarray(out_i, dtype=int)
    if hap_index is None:
        sub = gt_arr[idx, ind_index : ind_index + 1, :]
        g = allel.GenotypeArray(sub)
        dos = g.to_n_alt().astype(float)
        dos[g.is_missing()] = -1
        out[out_i] = dos[:, 0]
        return out

    alleles = gt_arr[idx, ind_index, hap_index].astype(float)
    called = alleles >= 0
    dos = np.full(len(idx), -1.0, dtype=float)
    dos[called] = alleles[called] * P
    out[out_i] = dos

    return out


def _dosage_for_positions(
    gt: allel.GenotypeArray,
    pos: Union[list, np.ndarray],
    ind_index: int,
    positions: list,
    hap_index: Optional[int] = None,
    P: int = 2,
) -> np.ndarray:
    """
    Return ALT dosage for one individual at exact genomic positions.

    If a requested position is absent or the genotype is missing, the returned
    dosage is -1. When extracting a single haplotype, called alleles are
    multiplied by ``P`` so that haplotype allele states are represented on the
    same dosage scale used for diploid genotypes.

    Parameters
    ----------
    gt : allel.GenotypeArray
        Genotype array with variant, sample, and ploidy axes.
    pos : list or numpy.ndarray
        Variant positions corresponding to gt.
    ind_index : int
        Index of the individual to extract.
    positions : list
        Genomic positions to query.
    hap_index : int or None, optional
        Zero-based haplotype index to extract for phased targets. When None,
        return diploid ALT dosage using scikit-allel's to_n_alt().
    P : int, optional
        Ploidy multiplier for haploid allele extraction. Default: 2, which
        maps haplotype allele 0 to dosage 0 and haplotype allele 1 to dosage 2.

    Returns
    -------
    numpy.ndarray
        Dosage array aligned to positions. Missing genotypes and absent
        positions are encoded as -1.
    """
    pos_to_i = _build_position_index(pos)
    return _dosage_for_position_index(
        gt=gt,
        pos_to_i=pos_to_i,
        ind_index=ind_index,
        positions=positions,
        hap_index=hap_index,
        P=P,
    )


def match(
    tgt_gt: allel.GenotypeArray,
    tgt_pos: Union[list, np.ndarray],
    tgt_index: int,
    src_gt: allel.GenotypeArray,
    src_pos: Union[list, np.ndarray],
    src_index: int,
    positions: list,
    denominator: int,
    P: int = 2,
    tgt_hap_index: Optional[int] = None,
) -> Union[float, str]:
    """
    Calculate pairwise match rate between one target and one source individual.

    Target and source dosages are extracted at the requested positions, then
    compared with ``calc_match_pct`` using the provided window-level
    denominator. If ``tgt_hap_index`` is provided, the target dosage is derived
    from that haplotype only.

    Parameters
    ----------
    tgt_gt : allel.GenotypeArray
        Target genotype array.
    tgt_pos : list or numpy.ndarray
        Variant positions corresponding to ``tgt_gt``.
    tgt_index : int
        Index of the target individual in ``tgt_gt``.
    src_gt : allel.GenotypeArray
        Source genotype array.
    src_pos : list or numpy.ndarray
        Variant positions corresponding to ``src_gt``.
    src_index : int
        Index of the source individual in ``src_gt``.
    positions : list
        Window-level variant positions to compare.
    denominator : int
        Number of population-level variants in the window.
    P : int, optional
        Ploidy used to normalize dosage differences. Default: 2.
    tgt_hap_index : int or None, optional
        Zero-based target haplotype index to compare. When None, compare the
        target diploid genotype dosage.

    Returns
    -------
    float or str
        Window-normalized pairwise match rate, or "NA" if the denominator is
        zero.
    """
    tgt_dos = _dosage_for_positions(
        tgt_gt,
        tgt_pos,
        tgt_index,
        positions,
        hap_index=tgt_hap_index,
        P=P,
    )
    src_dos = _dosage_for_positions(
        src_gt,
        src_pos,
        src_index,
        positions,
        P=P,
    )

    return calc_match_pct(tgt_dos, src_dos, denominator=denominator, P=P)


def _base_sample_name(sample: str, samples: list[str]) -> tuple[str, Optional[int]]:
    """
    Resolve a target tract label to a base sample and optional haplotype index.

    Labels that exactly match an entry in ``samples`` are treated as unphased
    diploid samples. Labels ending in ``_1`` or ``_2`` are treated as phased
    haplotype labels when the prefix matches a sample name.

    Parameters
    ----------
    sample : str
        Target tract sample label.
    samples : list of str
        Valid target sample names.

    Returns
    -------
    tuple[str, int or None]
        Base sample name and zero-based haplotype index. The haplotype index is
        None for unphased diploid samples.

    Raises
    ------
    ValueError
        If ``sample`` cannot be resolved against ``samples``.
    """
    if sample in samples:
        return sample, None

    base, sep, suffix = sample.rpartition("_")
    if sep and suffix in {"1", "2"} and base in samples:
        return base, int(suffix) - 1

    raise ValueError(f"Target sample {sample} is not present in the target list.")


def _resolve_chrom(chrom: object, data: dict) -> str:
    """
    Resolve a BED chromosome label against VCF chromosome keys.

    If the label is not found exactly, this function tries the corresponding
    representation with or without a ``chr`` prefix.

    Parameters
    ----------
    chrom : object
        Chromosome label from the tract file.
    data : dict
        Genotype data keyed by chromosome label.

    Returns
    -------
    str
        Chromosome key present in ``data``.

    Raises
    ------
    ValueError
        If neither the exact label nor the alternate ``chr``-prefixed form is
        present in ``data``.
    """
    chrom = str(chrom)
    if chrom in data:
        return chrom

    if chrom.startswith("chr"):
        alt_chrom = chrom[3:]
    else:
        alt_chrom = f"chr{chrom}"

    if alt_chrom in data:
        return alt_chrom

    available = ", ".join(str(c) for c in data.keys())
    raise ValueError(
        f"Chromosome {chrom} is not present in the VCF data. "
        f"Available chromosomes: {available}"
    )


def _validate_disjoint_samples(tgt_samples: list[str], src_samples: list[str]) -> None:
    """
    Raise an error if target and source sample lists overlap.

    Parameters
    ----------
    tgt_samples : list of str
        Target sample names.
    src_samples : list of str
        Source sample names.

    Raises
    ------
    ValueError
        If any sample appears in both target and source sample lists.
    """
    overlap = set(tgt_samples) & set(src_samples)
    if overlap:
        overlap_text = ", ".join(sorted(overlap))
        raise ValueError(
            f"Target and source sample lists must not overlap: {overlap_text}"
        )


def _validate_sorted_positions(data: dict) -> None:
    """
    Raise an error if positions are not sorted within any chromosome.

    Parameters
    ----------
    data : dict
        Genotype data keyed by chromosome.

    Raises
    ------
    ValueError
        If variant positions decrease within a chromosome.
    """
    for chrom, chrom_data in data.items():
        pos = np.asarray(chrom_data["POS"])
        if np.any(pos[:-1] > pos[1:]):
            raise ValueError(f"VCF positions must be sorted within chromosome {chrom}.")


def _positions_in_bed_interval(pos: np.ndarray, start: int, end: int) -> np.ndarray:
    """
    Return VCF positions satisfying ``start < POS <= end``.

    Parameters
    ----------
    pos : numpy.ndarray
        Sorted VCF positions for one chromosome.
    start : int
        BED-style tract start.
    end : int
        BED-style tract end.

    Returns
    -------
    numpy.ndarray
        Positions inside the interval.
    """
    left = np.searchsorted(pos, start, side="right")
    right = np.searchsorted(pos, end, side="right")
    return pos[left:right]


def run_match(
    vcf_file: str,
    tgt_ind_file: str,
    src_ind_file: str,
    tract_file: str,
    output_file: str,
    ploidy: int = 2,
) -> None:
    """
    Calculate aggregated source match rates for inferred tracts.

    The tract file must contain at least four columns: chromosome, start, end,
    and target sample label. Tract coordinates are interpreted as BED-style
    0-based half-open intervals, so a 1-based VCF position belongs to a tract
    when ``start < POS <= end``. Target and source samples are loaded from the
    same VCF, so window variants are the VCF positions in the tract. For each
    source individual, the pairwise score is normalized by this window-level
    variant count, and the final tract score is the mean across all source
    individuals.

    Parameters
    ----------
    vcf_file : str
        Input VCF file containing target and source samples.
    tgt_ind_file : str
        File containing one target sample name per line.
    src_ind_file : str
        File containing one source sample name per line.
    tract_file : str
        BED-like tract file with at least chrom, start, end, and sample
        columns. Additional columns are ignored.
    output_file : str
        Output TSV path. The output columns are chrom, start, end, sample, and
        match_rate.
    ploidy : int, optional
        Ploidy used to normalize dosage differences and scale haplotype
        dosages. Default: 2.
    """
    tgt_samples = parse_ind_file(tgt_ind_file)
    src_samples = parse_ind_file(src_ind_file)
    _validate_disjoint_samples(tgt_samples, src_samples)

    analysis_samples = tgt_samples + src_samples
    data, returned_samples = read_geno_data(
        vcf_file, analysis_samples, None, filter_missing=False
    )
    sample_to_index = {str(sample): i for i, sample in enumerate(returned_samples)}
    src_indices = [sample_to_index[sample] for sample in src_samples]

    _validate_sorted_positions(data)

    pos_by_chrom = {
        chrom: np.asarray(chrom_data["POS"]) for chrom, chrom_data in data.items()
    }
    position_index_by_chrom = {
        chrom: _build_position_index(pos) for chrom, pos in pos_by_chrom.items()
    }

    tracts = pd.read_csv(tract_file, sep="\t", header=None)
    if tracts.shape[1] < 4:
        raise ValueError(
            "The tract file must contain at least four columns: "
            "chrom, start, end, and sample."
        )
    tracts = tracts.iloc[:, :4]
    tracts.columns = ["chrom", "start", "end", "sample"]

    rows = []
    for _, tract in tracts.iterrows():
        chrom = tract["chrom"]
        start = int(tract["start"])
        end = int(tract["end"])
        tract_sample = str(tract["sample"])
        base_sample, tgt_hap_index = _base_sample_name(tract_sample, tgt_samples)
        tgt_index = sample_to_index[base_sample]

        resolved_chrom = _resolve_chrom(chrom, data)
        gt = data[resolved_chrom]["GT"]
        pos = pos_by_chrom[resolved_chrom]
        pos_to_i = position_index_by_chrom[resolved_chrom]

        # Tract coordinates are BED-style 0-based half-open [start, end),
        # while VCF POS values are 1-based; therefore a VCF position belongs
        # to the tract when start < POS <= end.
        region_pos = _positions_in_bed_interval(pos, start, end)
        positions = region_pos.astype(int).tolist()
        window_variant_count = len(positions)

        if window_variant_count == 0:
            rate = "NA"
        else:
            tgt_dos = _dosage_for_position_index(
                gt=gt,
                pos_to_i=pos_to_i,
                ind_index=tgt_index,
                positions=positions,
                hap_index=tgt_hap_index,
                P=ploidy,
            )
            source_rates = []
            for src_index in src_indices:
                src_dos = _dosage_for_position_index(
                    gt=gt,
                    pos_to_i=pos_to_i,
                    ind_index=src_index,
                    positions=positions,
                    P=ploidy,
                )
                source_rates.append(
                    calc_match_pct(
                        tgt_dos,
                        src_dos,
                        denominator=window_variant_count,
                        P=ploidy,
                    )
                )
            rate = float(np.mean(source_rates))

        rows.append(
            {
                "chrom": chrom,
                "start": start,
                "end": end,
                "sample": tract_sample,
                "match_rate": rate,
            }
        )

    output_dir = os.path.dirname(output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    columns = ["chrom", "start", "end", "sample", "match_rate"]
    pd.DataFrame(rows, columns=columns).to_csv(
        output_file, sep="\t", index=False, header=False
    )
