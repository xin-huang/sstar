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

import numpy as np
from typing import Any
from sstar.sstar import Sstar


class FeatureVectorPreprocessor:
    """
    A preprocessor subclass for generating feature vectors from genomic data.

    This class extends DataPreprocessor to include additional functionality for creating
    feature vectors based on genomic variants, reference and target individual genotypes,
    and window-based genomic statistics.
    """

    def __init__(
        self,
        ref_samples: list[str],
        tgt_samples: list[str],
        match_bonus: int,
        max_mismatch: int,
        mismatch_penalty: int,
    ):
        """
        Initializes a new instance of FeatureVectorsPreprocessor with specific parameters.

        Parameters
        ----------
        ref_samples : list[str]
            Reference individual identifiers.
        tgt_samples : list[str]
            Target individual identifiers.
        match_bonus : int
            S* match bonus.
        max_mismatch : int
            S* maximum mismatches.
        mismatch_penalty : int
            S* mismatch penalty.
        """
        if len(ref_samples) == 0:
            raise ValueError("No reference sample is provided.")
        if len(tgt_samples) == 0:
            raise ValueError("No target sample is provided.")

        self.match_bonus = match_bonus
        self.max_mismatch = max_mismatch
        self.mismatch_penalty = mismatch_penalty
        self.samples = {
            "ref": ref_samples,
            "tgt": tgt_samples,
        }

    def run(
        self,
        chr_name: str,
        start: int,
        end: int,
        ploidy: int,
        is_phased: bool,
        ref_gts: np.ndarray,
        tgt_gts: np.ndarray,
        pos: np.ndarray,
    ) -> list[dict[str, Any]]:
        """
        Executes the feature vector generation process for a specified genomic window.

        Parameters
        ----------
        chr_name : str
            Name of the chromosome.
        start : int
            Start position of the genomic window.
        end : int
            End position of the genomic window.
        ploidy : int
            Ploidy of the samples, typically 2 for diploid organisms.
        is_phased : bool
            Indicates whether the genomic data is phased.
        ref_gts : np.ndarray
            Genotype array for the reference individuals.
        tgt_gts : np.ndarray
            Genotype array for the target individuals.
        pos : np.ndarray
            Array of variant positions within the genomic window.

        Returns
        -------
        list
            A list of dictionaries containing the formatted feature vectors for the genomic window.
        """
        params = {
            "ref_gts": ref_gts,
            "tgt_gts": tgt_gts,
            "pos": pos,
        }

        items = {
            "chr_name": chr_name,
            "start": start,
            "end": end,
        }

        stat_params = {}
        stat_params.update(**params)
        stat_params.update(
            match_bonus=self.match_bonus,
            max_mismatch=self.max_mismatch,
            mismatch_penalty=self.mismatch_penalty,
        )
        sstar_res = Sstar.compute(**stat_params)
        items.update(sstar_res)

        fmtted_res = self._fmt_res(
            res=items,
            ploidy=ploidy,
            is_phased=is_phased,
        )

        return fmtted_res

    def _fmt_res(
        self,
        res: dict[str, Any],
        ploidy: int,
        is_phased: bool,
    ) -> list[dict[str, Any]]:
        """
        Formats the result dictionaries into a pandas DataFrame with appropriate headers.

        Parameters
        ----------
        res : dict[str, Any]
            A dictionary representing the results for a genomic window.
        ploidy : int
            The ploidy of the samples being processed.
        is_phased : bool
            Indicates whether the genomic data is phased.

        Returns
        -------
        list
            A list of dictionaries containing the formatted results with one row per sample and one column per feature.
        """
        if is_phased:
            num_samples = len(self.samples["tgt"]) * ploidy
        else:
            num_samples = len(self.samples["tgt"])

        base_dict = {
            "Chromosome": res["chr_name"],
            "Start": res["start"],
            "End": res["end"],
        }
        sample_dicts = [base_dict.copy() for _ in range(num_samples)]
        for i, sample_dict in enumerate(sample_dicts):
            sample = (
                f'{self.samples["tgt"][i//ploidy]}_{i % ploidy+1}'
                if is_phased
                else self.samples["tgt"][i]
            )
            sample_dict["Sample"] = sample
            sample_dict["S*_score"] = res["sstar"][i]
            sample_dict["Region_ind_SNP_number"] = res["region_ind_SNP_number"][i]

        return sample_dicts
