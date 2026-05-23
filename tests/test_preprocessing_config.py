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

from sstar.configs.preprocessing_config import PreprocessingConfig


def test_preprocess_config_normalizes_output_file(tmp_path):
    rel_out = tmp_path / "nested" / "sample.features"
    cfg = PreprocessingConfig(
        vcf_file=tmp_path / "input.vcf",
        chr_name="chr1",
        ref_ind_file=tmp_path / "ref.txt",
        tgt_ind_file=tmp_path / "tgt.txt",
        output_file=rel_out,
        win_len=1000,
        win_step=500,
        feature_config_file=tmp_path / "feature.yml",
    )

    assert cfg.output_file == Path(rel_out).expanduser().resolve()
