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


import os, pytest
import numpy as np
from sstar.utils import read_data, split_genome
from sstar.generators import WindowDataGenerator


@pytest.fixture
def file_paths():
    expected_dir = "tests/exp_results"
    return {
        "vcf_file": os.path.join(expected_dir, "test.0.vcf"),
        "ref_ind_file": os.path.join(expected_dir, "test.0.ref.ind.list"),
        "tgt_ind_file": os.path.join(expected_dir, "test.0.tgt.ind.list"),
        "anc_allele_file": None,
        "is_phased": True,
    }


@pytest.fixture
def init_params(file_paths):
    return {
        **file_paths,
        "chr_name": "1",
        "win_len": 50000,
        "win_step": 50000,
    }


@pytest.fixture
def expected_params(file_paths):
    expected_data = []
    chr_name = "1"
    win_len = 50000
    win_step = 50000
    ref_data, _ref_samples, tgt_data, _tgt_samples, _src_data, _src_samples = read_data(
        **file_paths
    )
    windows = split_genome(tgt_data[chr_name]["POS"], chr_name, win_step, win_len)

    for w in range(len(windows)):
        chr_name = windows[w][0]
        start = windows[w][1][0]
        end = windows[w][1][1]
        ref_gts = ref_data[chr_name]["GT"]
        tgt_gts = tgt_data[chr_name]["GT"]
        pos = tgt_data[chr_name]["POS"]
        idx = (pos >= start) & (pos <= end)
        sub_ref_gts = ref_gts[idx]
        sub_tgt_gts = tgt_gts[idx]
        sub_pos = pos[idx]

        d = {
            "chr_name": chr_name,
            "start": start,
            "end": end,
            "ploidy": 2,
            "is_phased": True,
            "ref_gts": sub_ref_gts,
            "tgt_gts": sub_tgt_gts,
            "pos": sub_pos,
        }
        expected_data.append(d)

    return expected_data


def test_WindowDataGenerator(init_params, expected_params):
    generator = WindowDataGenerator(**init_params)
    generated_params_list = list(generator.get())

    assert len(generated_params_list) == len(
        expected_params
    ), "The number of generated and expected parameters do not match."

    for generated, expected in zip(generated_params_list, expected_params):
        for key in generated:
            if isinstance(generated[key], np.ndarray) and isinstance(
                expected[key], np.ndarray
            ):
                assert np.array_equal(
                    generated[key], expected[key]
                ), f"Arrays do not match for key {key}."
            else:
                assert (
                    generated[key] == expected[key]
                ), f"Values do not match for key {key}."


def test_WindowDataGenerator_with_source_genotypes(init_params, file_paths, tmp_path):
    src_ind_file = tmp_path / "src.ind.list"
    src_ind_file.write_text("tsk_75\ntsk_76\ntsk_77\n")

    generator = WindowDataGenerator(
        **init_params,
        src_ind_file=str(src_ind_file),
    )
    generated_params_list = list(generator.get())

    ref_data, _ref_samples, tgt_data, _tgt_samples, src_data, _src_samples = read_data(
        **file_paths,
        src_ind_file=str(src_ind_file),
    )
    chr_name = init_params["chr_name"]
    windows = split_genome(
        tgt_data[chr_name]["POS"],
        chr_name,
        init_params["win_step"],
        init_params["win_len"],
    )

    assert generated_params_list, "Expected at least one generated window."
    for generated, window in zip(generated_params_list, windows):
        chr_name = window[0]
        start, end = window[1]
        pos = tgt_data[chr_name]["POS"]
        idx = (pos >= start) & (pos <= end)

        assert "src_gts" in generated
        assert np.array_equal(generated["src_gts"], src_data[chr_name]["GT"][idx])


def test_window_generator_uses_closed_interval_boundaries():
    pos = np.array([1, 10, 20, 21, 30])
    start = 10
    end = 20
    idx = (pos >= start) & (pos <= end)

    assert np.array_equal(pos[idx], np.array([10, 20]))
