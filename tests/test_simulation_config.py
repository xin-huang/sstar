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

from pathlib import Path

from sstar.configs.simulation_config import SimulationConfig


def test_simulation_config_normalizes_output_file(tmp_path):
    rel_out = tmp_path / "nested" / "sim.features.tsv"
    cfg = SimulationConfig(
        nrep=1,
        nref=2,
        ntgt=2,
        ref_id="REF",
        tgt_id="TGT",
        ploidy=2,
        is_phased=True,
        seq_len=1000,
        mut_rate=1e-8,
        rec_rate=1e-8,
        nprocess=1,
        force_balanced=False,
        output_file=rel_out,
        keep_sim_data=False,
        seed=1,
        feature_config_file=tmp_path / "feature.yml",
        is_shuffled=True,
        nfeature=10,
    )

    assert cfg.output_file == Path(rel_out).expanduser().resolve()
