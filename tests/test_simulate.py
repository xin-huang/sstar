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


import pytest
import pandas as pd
from sstar.simulate import simulate


def test_simulate(tmp_path):
    output_dir = tmp_path
    output_prefix = "check.simulate"
    output_file = output_dir / f"{output_prefix}.training.features.tsv"

    simulate(
        demo_model_file="tests/data/ArchIE_3D19_wo_intro.yaml",
        nrep=100,
        nref=50,
        ntgt=50,
        ref_id="Reference",
        tgt_id="Target",
        ploidy=2,
        seq_len=50000,
        mut_rate=1.25e-8,
        rec_rate=1.0e-8,
        feature_config_file="tests/data/sstar.yaml",
        output_dir=output_dir,
        output_prefix=output_prefix,
        nfeature=10000,
        is_phased=False,
        is_shuffled=False,
        keep_sim_data=False,
        seed=12345,
        nprocess=2,
    )

    df = pd.read_csv(output_file, sep="\t")
    df_expected = pd.read_csv("tests/exp_results/test.simulate.expected.tsv", sep="\t")

    pd.testing.assert_frame_equal(
        df,
        df_expected,
        check_dtype=False,
        check_exact=False,
        rtol=1e-6,
        atol=1e-8,
    )
