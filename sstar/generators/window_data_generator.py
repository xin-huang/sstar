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


from sstar.utils import read_data, split_genome
from sstar.generators import GenericGenerator


class WindowDataGenerator(GenericGenerator):
    """
    Generates genomic data for each specified window from VCF and other related files.
    """

    def __init__(
        self,
        vcf_file: str,
        chr_name: str,
        ref_ind_file: str,
        tgt_ind_file: str,
        win_len: int,
        win_step: int,
        anc_allele_file: str = None,
        src_ind_file: str = None,
        ploidy: int = 2,
        is_phased: bool = True,
    ):
        """
        Initializes a new instance of WindowDataGenerator.

        Parameters
        ----------
        vcf_file : str
            The path to the VCF file containing variant data.
        chr_name : str
            The name of the chromosome to process data for.
        ref_ind_file : str
            The path to the file containing identifiers for reference individuals.
        tgt_ind_file : str
            The path to the file containing identifiers for target individuals.
        win_len : int
            The length of each window in base pairs.
        win_step : int
            The step size between windows in base pairs.
        anc_allele_file : str, optional
            The path to the file containing ancestral allele information.
            Default: None.
        src_ind_file : str, optional
            The path to the file containing identifiers for source individuals.
            Default: None.
        ploidy : int, optional
            The ploidy of the genome. Default: 2.
        is_phased : bool, optional
            Specifies whether the genotype data is phased. Default: True.

        Raises
        ------
        ValueError
            If `win_len` is less than or equal to 0, if `win_step` is negative,
            if `ploidy` is less than or equal to 0, or if `chr_name` is not in the VCF file.
        """
        if win_len <= 0:
            raise ValueError("win_len must be greater than 0.")

        if win_step < 0:
            raise ValueError("win_step must be non-negative.")

        if ploidy <= 0:
            raise ValueError("ploidy must be greater than 0.")

        self.ploidy = ploidy
        self.is_phased = is_phased

        ref_data, _ref_samples, tgt_data, _tgt_samples, src_data, _src_samples = (
            read_data(
                vcf_file,
                ref_ind_file,
                tgt_ind_file,
                anc_allele_file,
                is_phased,
                src_ind_file,
            )
        )

        if chr_name not in tgt_data:
            raise ValueError(f"{chr_name} is not present in the VCF file.")

        pos = tgt_data[chr_name]["POS"]

        windows = split_genome(
            pos=pos,
            chr_name=chr_name,
            polymorphism_size=win_len,
            step_size=win_step,
        )

        self.data = []
        for w in range(len(windows)):
            chr_name = windows[w][0]
            start = windows[w][1][0]
            end = windows[w][1][1]
            ref_gts = ref_data[chr_name]["GT"]
            tgt_gts = tgt_data[chr_name]["GT"]
            src_gts = None
            if src_data is not None and chr_name in src_data:
                src_gts = src_data[chr_name]["GT"]
            idx = (pos >= start) & (pos <= end)
            sub_ref_gts = ref_gts[idx]
            sub_tgt_gts = tgt_gts[idx]
            sub_src_gts = src_gts[idx] if src_gts is not None else None
            sub_pos = pos[idx]

            d = {
                "chr_name": chr_name,
                "start": start,
                "end": end,
                "ploidy": self.ploidy,
                "is_phased": self.is_phased,
                "ref_gts": sub_ref_gts,
                "tgt_gts": sub_tgt_gts,
                "pos": sub_pos,
            }
            if sub_src_gts is not None:
                d["src_gts"] = sub_src_gts

            self.data.append(d)

    def get(self):
        """
        Yields genomic data for each window.

        Yields
        ------
        dict
            A dictionary containing chromosome name, start and end positions,
            ploidy and phase information, reference, target, optional source
            genotypes, and positions for each window.
        """
        for d in self.data:
            yield d

    def __len__(self):
        return len(self.data)
