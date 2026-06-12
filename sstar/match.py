# Copyright 2026 Xin Huang
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
    P: int = 2,
) -> Union[float, str]:
    """
    Calculate dosage-based match rate between two genotype dosage arrays.

    Missing values must be encoded as -1 and are ignored.

    Parameters
    ----------
    x : float, list, or numpy.ndarray
        ALT-allele dosage values for the first genotype vector.
    y : float, list, or numpy.ndarray
        ALT-allele dosage values for the second genotype vector.
    P : int, optional
        Ploidy used to normalize dosage differences. Default: 2.

    Returns
    -------
    float or str
        Mean match rate in [0, 1], or "NA" if there are no callable sites.
    """
    x = np.atleast_1d(np.asarray(x, dtype=float))
    y = np.atleast_1d(np.asarray(y, dtype=float))

    ok = (x >= 0) & (y >= 0)
    if not np.any(ok):
        return "NA"

    return float(np.mean(1.0 - np.abs(x[ok] - y[ok]) / float(P)))


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
    dosage is -1.

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
        Haplotype index to extract for phased targets. When None, return
        diploid ALT dosage using scikit-allel's to_n_alt().
    P : int, optional
        Ploidy multiplier for haploid allele extraction. Default: 2.

    Returns
    -------
    numpy.ndarray
        Dosage array aligned to positions.
    """
    pos = np.asarray(pos)
    gt_arr = np.asarray(gt)
    pos_to_i = {int(p): i for i, p in enumerate(pos)}
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


def match(
    tgt_gt: allel.GenotypeArray,
    tgt_pos: Union[list, np.ndarray],
    tgt_index: int,
    src_gt: allel.GenotypeArray,
    src_pos: Union[list, np.ndarray],
    src_index: int,
    positions: list,
    P: int = 2,
    tgt_hap_index: Optional[int] = None,
) -> Union[float, str]:
    """
    Calculate pairwise match rate between one target and one source individual.
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

    return calc_match_pct(tgt_dos, src_dos, P)


def _base_sample_name(sample: str, samples: list[str]) -> tuple[str, Optional[int]]:
    """
    Resolve a target tract label to a base sample and optional haplotype index.
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


def run_match(
    vcf_file: str,
    tgt_ind_file: str,
    src_ind_file: str,
    tract_file: str,
    output_file: str,
    ploidy: int = 2,
) -> None:
    """
    Calculate source match rates for inferred tracts and write them to a TSV.
    """
    tgt_samples = parse_ind_file(tgt_ind_file)
    src_samples = parse_ind_file(src_ind_file)

    tgt_data = read_geno_data(vcf_file, tgt_samples, None, filter_missing=False)
    src_data = read_geno_data(vcf_file, src_samples, None, filter_missing=False)

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
        tgt_index = tgt_samples.index(base_sample)

        tgt_chrom = _resolve_chrom(chrom, tgt_data)
        src_chrom = _resolve_chrom(chrom, src_data)
        tgt_gt = tgt_data[tgt_chrom]["GT"]
        tgt_pos = np.asarray(tgt_data[tgt_chrom]["POS"])
        src_gt = src_data[src_chrom]["GT"]
        src_pos = np.asarray(src_data[src_chrom]["POS"])

        # Tract coordinates are BED-style 0-based half-open [start, end),
        # while VCF POS values are 1-based; therefore a VCF position belongs
        # to the tract when start < POS <= end.
        tgt_region_pos = tgt_pos[(tgt_pos > start) & (tgt_pos <= end)]
        src_region_pos = src_pos[(src_pos > start) & (src_pos <= end)]
        positions = np.union1d(tgt_region_pos, src_region_pos).astype(int).tolist()

        for src_index, src_sample in enumerate(src_samples):
            rate = match(
                tgt_gt=tgt_gt,
                tgt_pos=tgt_pos,
                tgt_index=tgt_index,
                src_gt=src_gt,
                src_pos=src_pos,
                src_index=src_index,
                positions=positions,
                tgt_hap_index=tgt_hap_index,
                P=ploidy,
            )
            rows.append(
                {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "sample": tract_sample,
                    "match_rate": rate,
                    "src_sample": src_sample,
                }
            )

    output_dir = os.path.dirname(output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    columns = ["chrom", "start", "end", "sample", "match_rate", "src_sample"]
    pd.DataFrame(rows, columns=columns).to_csv(output_file, sep="\t", index=False)
