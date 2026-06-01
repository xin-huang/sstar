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

from sstar.match import run_match
from sstar.parsers.argument_validation import existed_file, positive_int


def _run_match(args: argparse.Namespace) -> None:
    """
    Execute the match process with specified parameters.

    Parameters
    ----------
    args : argparse.Namespace
        A namespace object obtained from argparse, containing specified parameters.
    """
    run_match(
        vcf_file=args.vcf,
        tgt_ind_file=args.tgt,
        src_ind_file=args.src,
        tract_file=args.tract_file,
        output_file=args.output,
        ploidy=args.ploidy,
    )


def add_match_parser(subparsers: argparse.ArgumentParser) -> None:
    """
    Initializes and configures the command-line interface parser
    for the match subcommand.

    Parameters
    ----------
    subparsers : argparse.ArgumentParser
        A command-line interface parser to be configured.
    """
    parser = subparsers.add_parser(
        "match", help="Calculate source match rates for inferred tracts."
    )
    parser.add_argument(
        "--vcf",
        type=existed_file,
        required=True,
        help="Path to the VCF file.",
    )
    parser.add_argument(
        "--tgt",
        type=existed_file,
        required=True,
        help="Path to the target individual list.",
    )
    parser.add_argument(
        "--src",
        type=existed_file,
        required=True,
        help="Path to the source individual list.",
    )
    parser.add_argument(
        "--tract-file",
        type=existed_file,
        required=True,
        help="Path to the inferred tract BED file from `sstar2 infer`.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to the output match-rate BED file.",
    )
    parser.add_argument(
        "--ploidy",
        type=positive_int,
        default=2,
        help="Ploidy used to normalize dosage differences. Default: %(default)s.",
    )
    parser.set_defaults(runner=_run_match)
