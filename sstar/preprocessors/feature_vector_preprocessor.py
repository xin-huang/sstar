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

import yaml
import numpy as np
from typing import Any
from gaishi.configs import FeatureConfig
from gaishi.preprocessors import GenericPreprocessor
from gaishi.registries.stat_registry import STAT_REGISTRY
from gaishi.utils import parse_ind_file


class FeatureVectorPreprocessor(GenericPreprocessor):
    """
    A preprocessor subclass for generating feature vectors from genomic data.

    This class extends DataPreprocessor to include additional functionality for creating
    feature vectors based on genomic variants, reference and target individual genotypes,
    and window-based genomic statistics.

    """

    def __init__(self, ref_ind_file: str, tgt_ind_file: str, feature_config_file: str):
        """
        Initializes a new instance of FeatureVectorsPreprocessor with specific parameters.

        Parameters:
        -----------
        ref_ind_file : str
            Path to the file listing reference individual identifiers.
        tgt_ind_file : str
            Path to the file listing target individual identifiers.
        feature_config_file : str
            Path to the configuration file specifying the features to be computed.

        Raises
        ------
        FileNotFoundError
            If the feature configuration file is not found.
        ValueError
            If the feature configuration file is incorrectly formatted or does not contain any features.

        """
        try:
            with open(feature_config_file, "r") as f:
                config_dict = yaml.safe_load(f)
        except FileNotFoundError:
            raise FileNotFoundError(
                f"Feature configuration file {feature_config_file} not found."
            )
        except yaml.YAMLError as exc:
            raise ValueError(f"Error parsing feature configuration: {exc}")

        self.feature_config = FeatureConfig(**config_dict)
        ref_samples = parse_ind_file(ref_ind_file)
        tgt_samples = parse_ind_file(tgt_ind_file)
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
            "is_phased": is_phased,
            "ploidy": ploidy,
            "pos": pos,
        }

        items = dict()
        items["chr_name"] = chr_name
        items["start"] = start
        items["end"] = end

        for stat_name in self.feature_config.root.keys():
            if not self.feature_config.root[stat_name]:
                continue
            stat_cls = STAT_REGISTRY.get(stat_name)
            stat_params = {}
            stat_params.update(**params)
            if isinstance(self.feature_config.root[stat_name], dict):
                stat_params.update(**self.feature_config.root[stat_name])
            items.update(stat_cls.compute(**stat_params))

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
                f'{self.samples["tgt"][i//ploidy]}_{i%ploidy+1}'
                if is_phased
                else self.samples["tgt"][i]
            )
            sample_dict["Sample"] = sample

            if "sstar" in res.keys():
                sample_dict["Sstar"] = res["sstar"][i]
            if "private_mutation" in res.keys():
                sample_dict["Private_mutation"] = res["private_mutation"][i]
            if "spectrum" in res.keys():
                for j, value in enumerate(res["spectrum"][i]):
                    sample_dict[f"{j}_ton"] = value

            for pop in ["ref", "tgt"]:
                if f"{pop}_dist" in self.feature_config.root.keys():
                    sample_dict.update(
                        self._fmt_dist_res(
                            row=res,
                            idx=i,
                            is_phased=is_phased,
                            ploidy=ploidy,
                            pop=pop,
                        )
                    )

        return sample_dicts

    def _fmt_dist_res(
        self,
        row: dict[str, Any],
        idx: int,
        is_phased: bool,
        ploidy: int,
        pop: str,
    ) -> dict[str, float]:
        """
        Formats distance-related results for a single sample and population (reference or target).

        Parameters
        ----------
        row : dict[str, Any]
            A dictionary containing results for a genomic window.
        idx : int
            The index of the sample in the target list.
        is_phased : bool
            Indicates whether the genomic data is phased.
        ploidy : int
            The ploidy of the samples being processed.
        pop : str
            Indicates the population type ('Ref' or 'Tgt') for which distances are being formatted.

        Returns
        -------
        dict[str, float]
            A dictionary with keys for each distance-related feature and values for the current sample.

        """
        items = dict()

        for feature in [
            "Minimum",
            "Maximum",
            "Mean",
            "Median",
            "Variance",
            "Skew",
            "Kurtosis",
        ]:
            if f"{feature}_{pop}_dist" in row.keys():
                items[f"{feature}_{pop}_dist"] = row[f"{feature}_{pop}_dist"][idx]

        if f"All_{pop}_dist" in row.keys():
            if is_phased:
                for j in range(len(row[f"All_{pop}_dist"][idx])):
                    sample = self.samples[pop][int(j / ploidy)]
                    items[f"{pop}_dist_{sample}_{j%ploidy+1}"] = row[f"All_{pop}_dist"][
                        idx
                    ][j]
            else:
                for j in range(len(row[f"All_{pop}_dist"][idx])):
                    sample = self.samples[pop][j]
                    items[f"{pop}_dist_{sample}"] = row[f"All_{pop}_dist"][idx][j]

        return items
