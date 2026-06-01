# Copyright 2026 Xin Huang and Andrea Koca
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

import argparse

import pandas as pd
import pytest

from sstar.assign import assign_source_by_match_rate
from sstar.parsers.assign_parser import add_assign_parser


def _write_match_rate_file(path, rows):
    pd.DataFrame(
        rows, columns=["chrom", "start", "end", "sample", "match_rate"]
    ).to_csv(path, sep="\t", index=False, header=False)


def _read_bed(path):
    return pd.read_csv(path, sep="\t", header=None, keep_default_na=False)


def test_assign_source_by_unique_maximum_match_rate(tmp_path):
    src1_file = tmp_path / "src1.match.tsv"
    src2_file = tmp_path / "src2.match.tsv"
    src3_file = tmp_path / "src3.match.tsv"
    output_prefix = tmp_path / "out"

    rows = [
        ["1", 100, 200, "ind1", 0.8],
        ["1", 300, 400, "ind1", 0.2],
        ["1", 500, 600, "ind1", 0.5],
    ]
    _write_match_rate_file(src1_file, rows)
    _write_match_rate_file(
        src2_file,
        [
            ["1", 100, 200, "ind1", 0.5],
            ["1", 300, 400, "ind1", 0.9],
            ["1", 500, 600, "ind1", 0.5],
        ],
    )
    _write_match_rate_file(
        src3_file,
        [
            ["1", 100, 200, "ind1", 0.1],
            ["1", 300, 400, "ind1", 0.1],
            ["1", 500, 600, "ind1", 0.5],
        ],
    )

    assign_source_by_match_rate(
        match_rate_files=[str(src1_file), str(src2_file), str(src3_file)],
        source_names=["src1", "src2", "src3"],
        output_prefix=str(output_prefix),
    )

    src1_bed = _read_bed(tmp_path / "out.src1.inferred.tracts.bed")
    src2_bed = _read_bed(tmp_path / "out.src2.inferred.tracts.bed")

    assert src1_bed.values.tolist() == [[1, 100, 200, "ind1", 0.8]]
    assert src2_bed.values.tolist() == [[1, 300, 400, "ind1", 0.9]]
    assert (tmp_path / "out.src3.inferred.tracts.bed").read_text() == ""


def test_assign_source_uses_intersection_of_tract_keys(tmp_path):
    src1_file = tmp_path / "src1.match.tsv"
    src2_file = tmp_path / "src2.match.tsv"
    output_prefix = tmp_path / "out"

    _write_match_rate_file(
        src1_file,
        [
            ["1", 100, 200, "ind1", 0.8],
            ["1", 300, 400, "ind1", 0.9],
        ],
    )
    _write_match_rate_file(
        src2_file,
        [
            ["1", 100, 200, "ind1", 0.7],
            ["1", 500, 600, "ind1", 1.0],
        ],
    )

    assign_source_by_match_rate(
        match_rate_files=[str(src1_file), str(src2_file)],
        source_names=["src1", "src2"],
        output_prefix=str(output_prefix),
    )

    src1_bed = _read_bed(tmp_path / "out.src1.inferred.tracts.bed")
    assert src1_bed.values.tolist() == [[1, 100, 200, "ind1", 0.8]]
    assert (tmp_path / "out.src2.inferred.tracts.bed").read_text() == ""


def test_assign_source_ignores_na_and_rejects_ties(tmp_path):
    src1_file = tmp_path / "src1.match.tsv"
    src2_file = tmp_path / "src2.match.tsv"
    src3_file = tmp_path / "src3.match.tsv"
    output_prefix = tmp_path / "out"

    _write_match_rate_file(
        src1_file,
        [
            ["1", 100, 200, "ind1", "NA"],
            ["1", 300, 400, "ind1", 0.5],
            ["1", 500, 600, "ind1", "NA"],
        ],
    )
    _write_match_rate_file(
        src2_file,
        [
            ["1", 100, 200, "ind1", 0.4],
            ["1", 300, 400, "ind1", 0.5],
            ["1", 500, 600, "ind1", "NA"],
        ],
    )
    _write_match_rate_file(
        src3_file,
        [
            ["1", 100, 200, "ind1", "NA"],
            ["1", 300, 400, "ind1", 0.1],
            ["1", 500, 600, "ind1", "NA"],
        ],
    )

    assign_source_by_match_rate(
        match_rate_files=[str(src1_file), str(src2_file), str(src3_file)],
        source_names=["src1", "src2", "src3"],
        output_prefix=str(output_prefix),
    )

    src2_bed = _read_bed(tmp_path / "out.src2.inferred.tracts.bed")
    assert src2_bed.values.tolist() == [[1, 100, 200, "ind1", 0.4]]
    assert (tmp_path / "out.src1.inferred.tracts.bed").read_text() == ""
    assert (tmp_path / "out.src3.inferred.tracts.bed").read_text() == ""


def test_assign_source_validates_inputs(tmp_path):
    src1_file = tmp_path / "src1.match.tsv"
    src2_file = tmp_path / "src2.match.tsv"
    _write_match_rate_file(src1_file, [["1", 100, 200, "ind1", 0.8]])
    _write_match_rate_file(src2_file, [["1", 100, 200, "ind1", 0.7]])

    with pytest.raises(ValueError, match="At least two"):
        assign_source_by_match_rate([str(src1_file)], ["src1"], str(tmp_path / "out"))

    with pytest.raises(ValueError, match="must match source names"):
        assign_source_by_match_rate(
            [str(src1_file), str(src2_file)], ["src1"], str(tmp_path / "out")
        )

    with pytest.raises(ValueError, match="must be unique"):
        assign_source_by_match_rate(
            [str(src1_file), str(src2_file)],
            ["src1", "src1"],
            str(tmp_path / "out"),
        )

    with pytest.raises(ValueError, match="reserved column names"):
        assign_source_by_match_rate(
            [str(src1_file), str(src2_file)],
            ["src1", "sample"],
            str(tmp_path / "out"),
        )


def test_assign_parser_accepts_required_args(tmp_path):
    src1_file = tmp_path / "src1.match.tsv"
    src2_file = tmp_path / "src2.match.tsv"
    _write_match_rate_file(src1_file, [["1", 100, 200, "ind1", 0.8]])
    _write_match_rate_file(src2_file, [["1", 100, 200, "ind1", 0.7]])

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="subparsers")
    add_assign_parser(subparsers)

    args = parser.parse_args(
        [
            "assign",
            "--match-rate",
            str(src1_file),
            str(src2_file),
            "--source-name",
            "src1",
            "src2",
            "--output-prefix",
            str(tmp_path / "out"),
        ]
    )

    assert args.match_rate == [str(src1_file), str(src2_file)]
    assert args.source_name == ["src1", "src2"]
    assert args.output_prefix == str(tmp_path / "out")
    assert hasattr(args, "runner")
    assert callable(args.runner)
