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
import pybedtools

def get_tract(threshold_file, match_pct_files, output_prefix, diff):
    """
    """
    threshold_df = pd.read_csv(threshold_file, sep="\t")
    threshold_df = threshold_df[threshold_df['significant'] == True]
    cols = ['chrom', 'start', 'end']

    if match_pct_files is None:
        pybedtools.BedTool.from_dataframe(threshold_df).sort().merge().moveto(output_prefix+'.bed')
    else:
        src1_df = pd.read_csv(match_pct_files[0], sep="\t")
        src2_df = pd.read_csv(match_pct_files[1], sep="\t")
        src1_sig_df = pd.merge(src1_df, threshold_df, on=['chrom', 'start', 'end', 'sample'])
        src2_sig_df = pd.merge(src2_df, threshold_df, on=['chrom', 'start', 'end', 'sample'])
        df = pd.merge(src1_sig_df, src2_sig_df, on=['chrom', 'start', 'end', 'sample'])
        src1_df = df[df['match_pct_x']-df['match_pct_y']>diff]
        src2_df = df[df['match_pct_x']-df['match_pct_y']<diff]
        pybedtools.BedTool.from_dataframe(src1_df).sort().merge().moveto(output_prefix+'.src1.bed')
        pybedtools.BedTool.from_dataframe(src2_df).sort().merge().moveto(output_prefix+'.src2.bed')
