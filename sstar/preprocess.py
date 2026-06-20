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
import pandas as pd
from sstar.mp_manager import mp_manager
from sstar.generators import WindowDataGenerator
from sstar.feature_vector_preprocessor import FeatureVectorPreprocessor
from sstar.utils import parse_ind_file


def preprocess(
    vcf_file: str,
    chr_name: str,
    ref_ind_file: str,
    tgt_ind_file: str,
    win_len: int,
    win_step: int,
    output_file: str,
    match_bonus: int,
    max_mismatch: int,
    mismatch_penalty: int,
    nprocess: int = 1,
    ploidy: int = 2,
    is_phased: bool = True,
    anc_allele_file: str = None,
) -> None:
    """
    Preprocess genomic data to generate feature vectors for machine learning models.

    This function orchestrates the preprocessing pipeline by initializing a genomic data generator
    and a feature vector preprocessor. It utilizes multiprocessing to efficiently process large
    genomic datasets, generating feature vectors based on the specified configuration.

    Parameters
    ----------
    vcf_file : str
        Path to the VCF file containing genomic variants.
    chr_name : str
        Name of the chromosome to process.
    ref_ind_file : str
        Path to the file listing reference individual identifiers.
    tgt_ind_file : str
        Path to the file listing target individual identifiers.
    win_len : int
        Length of the sliding window for analysis, in base pairs.
    win_step : int
        Step size for the sliding window, in base pairs.
    output_file : str
        Output feature file path.
    match_bonus : int
        S* match bonus.
    max_mismatch : int
        S* maximum mismatches.
    mismatch_penalty : int
        S* mismatch penalty.
    nprocess : int, optional
        Number of worker processes to use for parallel processing. Default: 1.
    ploidy : int, optional
        Ploidy of the samples, typically 2 for diploid organisms. Default: 2.
    is_phased : bool, optional
        Indicates whether the genomic data is phased. Default: True.
    anc_allele_file : str, optional
        Path to the file containing ancestral allele information. Default: None.

    Raises
    ------
    ValueError
        If any of the provided parameters are invalid, such as negative window lengths or step sizes.
    """
    if nprocess <= 0:
        raise ValueError("Number of processes must be greater than 0.")

    generator = WindowDataGenerator(
        vcf_file=vcf_file,
        ref_ind_file=ref_ind_file,
        tgt_ind_file=tgt_ind_file,
        anc_allele_file=anc_allele_file,
        ploidy=ploidy,
        is_phased=is_phased,
        chr_name=chr_name,
        win_len=win_len,
        win_step=win_step,
    )

    ref_samples = parse_ind_file(ref_ind_file)
    tgt_samples = parse_ind_file(tgt_ind_file)

    preprocessor = FeatureVectorPreprocessor(
        ref_samples=ref_samples,
        tgt_samples=tgt_samples,
        match_bonus=match_bonus,
        max_mismatch=max_mismatch,
        mismatch_penalty=mismatch_penalty,
    )
    res = mp_manager(job=preprocessor, data_generator=generator, nprocess=nprocess)

    if res == "error":
        raise SystemExit("Some errors occurred, stopping the program ...")

    res.sort(key=lambda x: (x["Chromosome"], x["Start"], x["End"]))

    output_dir = os.path.dirname(output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    pd.DataFrame(res).to_csv(output_file, sep="\t", index=False, na_rep="NA")
