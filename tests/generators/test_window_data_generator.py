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
    ref_data, ref_samples, tgt_data, tgt_samples = read_data(**file_paths)
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


def test_window_generator_uses_closed_interval_boundaries():
    pos = np.array([1, 10, 20, 21, 30])
    start = 10
    end = 20
    idx = (pos >= start) & (pos <= end)

    assert np.array_equal(pos[idx], np.array([10, 20]))


def test_window_generator_scopes_alignment_to_requested_chromosome(tmp_path):
    vcf_file = tmp_path / "multi_chrom_no_shared.vcf"
    ref_ind_file = tmp_path / "ref.ind.list"
    tgt_ind_file = tmp_path / "tgt.ind.list"

    vcf_file.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tref1\ttgt1",
                "1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|0",
                "1\t200\t.\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|1",
                "2\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t.|.",
                "2\t200\t.\tA\tG\t.\tPASS\t.\tGT\t.|.\t0|1",
            ]
        )
        + "\n"
    )
    ref_ind_file.write_text("ref1\n")
    tgt_ind_file.write_text("tgt1\n")

    generator = WindowDataGenerator(
        vcf_file=str(vcf_file),
        chr_name="1",
        ref_ind_file=str(ref_ind_file),
        tgt_ind_file=str(tgt_ind_file),
        win_len=1000,
        win_step=1000,
    )

    generated_windows = list(generator.get())
    assert len(generated_windows) == 1
    assert generated_windows[0]["chr_name"] == "1"
    assert np.array_equal(generated_windows[0]["pos"], np.array([100, 200]))
