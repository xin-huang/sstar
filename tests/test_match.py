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

import allel
import numpy as np
import pandas as pd
import pytest

from sstar.__main__ import _sstar_cli_parser
from sstar.match import (
    _base_sample_name,
    _dosage_for_positions,
    _resolve_chrom,
    calc_match_pct,
    match,
    run_match,
)


def _write_match_input_files(tmp_path, tract_text):
    vcf_file = tmp_path / "mini.vcf"
    tgt_ind_file = tmp_path / "tgt.ind.list"
    src_ind_file = tmp_path / "src.ind.list"
    tract_file = tmp_path / "tract.bed"
    output_file = tmp_path / "nested" / "match.tsv"

    vcf_file.write_text(
        "##fileformat=VCFv4.2\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tind1\tsrc1\tsrc2\n"
        "21\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|0\t1|1\n"
        "21\t200\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t1|1\t0|1\n"
        "21\t300\t.\tA\tG\t.\tPASS\t.\tGT\t1|1\t1|1\t.|.\n"
    )
    tgt_ind_file.write_text("ind1\n")
    src_ind_file.write_text("src1\nsrc2\n")
    tract_file.write_text(tract_text)

    return vcf_file, tgt_ind_file, src_ind_file, tract_file, output_file


def test_calc_match_pct():
    assert calc_match_pct([0, 1, 2], [0, 1, 2], P=2) == 1.0
    assert calc_match_pct([0], [1], P=2) == 0.5
    assert calc_match_pct([0], [2], P=2) == 0.0
    assert calc_match_pct([0, 1, 2], [0, 2, 2], P=2) == pytest.approx(
        (1.0 + 0.5 + 1.0) / 3
    )
    assert calc_match_pct([0, -1], [0, 2], P=2) == 1.0
    assert calc_match_pct([-1], [-1], P=2) == "NA"


def test_dosage_for_positions():
    gt = allel.GenotypeArray(
        np.array(
            [
                [[0, 0]],
                [[0, 1]],
                [[1, 1]],
                [[-1, -1]],
            ]
        )
    )
    pos = np.array([100, 200, 300, 400])

    dosage = _dosage_for_positions(
        gt=gt,
        pos=pos,
        ind_index=0,
        positions=[300, 100, 250, 400, 200],
    )

    np.testing.assert_array_equal(dosage, np.array([2.0, 0.0, -1.0, -1.0, 1.0]))


def test_match():
    tgt_gt = allel.GenotypeArray(
        np.array(
            [
                [[0, 0]],
                [[0, 1]],
                [[1, 1]],
            ]
        )
    )
    src_gt = allel.GenotypeArray(
        np.array(
            [
                [[0, 0]],
                [[1, 1]],
                [[1, 1]],
            ]
        )
    )
    pos = np.array([100, 200, 300])

    assert match(tgt_gt, pos, 0, src_gt, pos, 0, [100, 200, 300], P=2) == pytest.approx(
        (1.0 + 0.5 + 1.0) / 3
    )


def test_base_sample_name_maps_phased_labels():
    samples = ["ind1", "ind2"]

    assert _base_sample_name("ind1", samples) == "ind1"
    assert _base_sample_name("ind1_1", samples) == "ind1"
    assert _base_sample_name("ind1_2", samples) == "ind1"

    with pytest.raises(ValueError):
        _base_sample_name("missing_1", samples)


def test_resolve_chrom_handles_chr_prefix():
    data = {"21": {}, "chr22": {}}

    assert _resolve_chrom("21", data) == "21"
    assert _resolve_chrom("chr21", data) == "21"
    assert _resolve_chrom("22", data) == "chr22"

    with pytest.raises(ValueError):
        _resolve_chrom("23", data)


def test_run_match(tmp_path):
    vcf_file, tgt_ind_file, src_ind_file, tract_file, output_file = _write_match_input_files(
        tmp_path, "21\t0\t250\tind1_1\n"
    )

    run_match(
        vcf_file=str(vcf_file),
        tgt_ind_file=str(tgt_ind_file),
        src_ind_file=str(src_ind_file),
        tract_file=str(tract_file),
        output_file=str(output_file),
        ploidy=2,
    )

    df = pd.read_csv(output_file, sep="\t")
    assert df.columns.tolist() == [
        "chrom",
        "start",
        "end",
        "sample",
        "match_rate",
        "src_sample",
    ]
    assert df["sample"].tolist() == ["ind1_1", "ind1_1"]

    rates = dict(zip(df["src_sample"], df["match_rate"]))
    assert rates["src1"] == pytest.approx(0.75)
    assert rates["src2"] == pytest.approx(0.5)


def test_run_match_ignores_extra_tract_columns(tmp_path):
    vcf_file, tgt_ind_file, src_ind_file, tract_file, output_file = _write_match_input_files(
        tmp_path, "21\t0\t250\tind1_1\textra_value\n"
    )

    run_match(
        vcf_file=str(vcf_file),
        tgt_ind_file=str(tgt_ind_file),
        src_ind_file=str(src_ind_file),
        tract_file=str(tract_file),
        output_file=str(output_file),
        ploidy=2,
    )

    df = pd.read_csv(output_file, sep="\t")
    assert df.columns.tolist() == [
        "chrom",
        "start",
        "end",
        "sample",
        "match_rate",
        "src_sample",
    ]
    assert df["sample"].tolist() == ["ind1_1", "ind1_1"]


def test_match_parser_accepts_required_args(tmp_path):
    vcf_file, tgt_ind_file, src_ind_file, tract_file, output_file = _write_match_input_files(
        tmp_path, "21\t0\t250\tind1_1\n"
    )
    parser = _sstar_cli_parser()

    args = parser.parse_args(
        [
            "match",
            "--vcf",
            str(vcf_file),
            "--tgt",
            str(tgt_ind_file),
            "--src",
            str(src_ind_file),
            "--tract-file",
            str(tract_file),
            "--output",
            str(output_file),
        ]
    )

    assert args.vcf == str(vcf_file)
    assert args.tgt == str(tgt_ind_file)
    assert args.src == str(src_ind_file)
    assert args.tract_file == str(tract_file)
    assert args.output == str(output_file)
    assert args.ploidy == 2
