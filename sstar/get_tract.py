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

import pandas as pd


def get_tract(threshold_file, match_pct_files, output_prefix, diff):
    """
    Description:
        Gets candidate introgressed fragments.

    Arguments:
        threshold_file str: Name of the file containing thresholds.
        match_pct_files list: List containing names of files containing source match rates.
        output_prefix str: Prefix of the output files.
        diff float: Cut-off to determine the origin of the candidate introgressed fragments.
    """
    threshold_df = pd.read_csv(threshold_file, sep="\t")
    threshold_df = threshold_df[threshold_df['significant'] == True]

    if match_pct_files is None:
        _output_bed(threshold_df, output_prefix+'.bed')
    elif len(match_pct_files) == 1:
        src_df = pd.read_csv(match_pct_files[0], sep="\t")
        merge_df = pd.merge(src_df, threshold_df, on=['chrom', 'start', 'end', 'sample'])
        _output_bed(merge_df, output_prefix+'.bed', col='match_rate')
    else:
        src1_df = pd.read_csv(match_pct_files[0], sep="\t")
        src2_df = pd.read_csv(match_pct_files[1], sep="\t")
        src1_sig_df = pd.merge(src1_df, threshold_df, on=['chrom', 'start', 'end', 'sample'])
        src2_sig_df = pd.merge(src2_df, threshold_df, on=['chrom', 'start', 'end', 'sample'])
        df = pd.merge(src1_sig_df, src2_sig_df, on=['chrom', 'start', 'end', 'sample'])
        src1_df = df[df['match_rate_x']-df['match_rate_y']>diff]
        src2_df = df[df['match_rate_x']-df['match_rate_y']<diff]
        _output_bed(src1_df, output_prefix+'.src1.bed', col='match_rate_x')
        _output_bed(src2_df, output_prefix+'.src2.bed', col='match_rate_y')


def _output_bed(df, output, col=None):
    """
    Description:
        Helper fuction for outputing candidate introgressed fragments in BED format.

    Arguments:
        df pandas.dataframe: Dataframe containing candidate introgressed fragments.
        output str: Name of the output file.
        col str: Name of the column containing source match rates in the dataframe.
    """
    df = df.sort_values(by=['sample', 'chrom', 'start', 'end'])
    if col is None:
        colnames = ['chrom', 'start', 'end', 'sample']
    else: colnames = ['chrom', 'start', 'end', 'sample', col]
    df.to_csv(output, columns=colnames, sep="\t", header=False, index=False)
