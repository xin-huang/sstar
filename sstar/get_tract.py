# Copyright 2025 Xin Huang
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


from typing import List, Optional
import pandas as pd


def get_tract(
    threshold_file: str,
    match_pct_files: Optional[List[str]],
    output_prefix: str,
    diff: float,
) -> None:
    """
    Gets candidate introgressed fragments based on threshold and match percentage data.

    Parameters
    ----------
    threshold_file : str
        Path to the file containing threshold values.
    match_pct_files : list of str or None
        List of file paths containing match percentages to source populations.
        If None, only threshold filtering is applied.
        Must contain either 1 or 2 file paths if provided.
    output_prefix : str
        Prefix for the output file(s).
    diff : float
        Cut-off used to determine source assignment when two source files are provided.
    """
    threshold_df = pd.read_csv(threshold_file, sep="\t")

    required_cols = {"chrom", "start", "end", "sample"}
    if not required_cols.issubset(threshold_df.columns):
        raise ValueError(f"`{threshold_file}` must contain columns: {required_cols}")

    significant_regions = threshold_df[threshold_df["significant"] == True]

    if match_pct_files is None:
        _output_bed(significant_regions, output_prefix + ".bed")
    elif len(match_pct_files) == 1:
        src_df = pd.read_csv(match_pct_files[0], sep="\t")
        merged_df = pd.merge(
            src_df, significant_regions, on=["chrom", "start", "end", "sample"]
        )
        _output_bed(merged_df, output_prefix + ".bed", col="match_rate")
    elif len(match_pct_files) == 2:
        src1_df = pd.read_csv(match_pct_files[0], sep="\t")
        src2_df = pd.read_csv(match_pct_files[1], sep="\t")

        src1_sig_df = pd.merge(
            src1_df, significant_regions, on=["chrom", "start", "end", "sample"]
        )
        src2_sig_df = pd.merge(
            src2_df, significant_regions, on=["chrom", "start", "end", "sample"]
        )
        merged_df = pd.merge(
            src1_sig_df, src2_sig_df, on=["chrom", "start", "end", "sample"]
        )

        assigned_to_src1 = merged_df[
            merged_df["match_rate_x"] - merged_df["match_rate_y"] > diff
        ]
        assigned_to_src2 = merged_df[
            merged_df["match_rate_x"] - merged_df["match_rate_y"] < -diff
        ]

        _output_bed(assigned_to_src1, output_prefix + ".src1.bed", col="match_rate_x")
        _output_bed(assigned_to_src2, output_prefix + ".src2.bed", col="match_rate_y")
    else:
        raise ValueError(
            "`match_pct_files` must be None, a list of 1 file path, or a list of 2 file paths."
        )


def _output_bed(df: pd.DataFrame, output: str, col: Optional[str] = None) -> None:
    """
    Outputs candidate introgressed fragments in BED format.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing candidate introgressed fragments.
    output : str
        Path to the output BED file.
    col : str or None, optional
        Name of the column containing source match rates to include in the output.
        If None, only basic BED fields are written.
    """
    df = df.sort_values(by=["sample", "chrom", "start", "end"])
    if col is None:
        colnames = ["chrom", "start", "end", "sample"]
    else:
        colnames = ["chrom", "start", "end", "sample", col]
    df.to_csv(output, columns=colnames, sep="\t", header=False, index=False)
