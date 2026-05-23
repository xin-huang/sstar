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

import argparse
from sstar.parsers.argument_validation import (
    existed_file,
    non_negative_int,
    non_positive_int,
    positive_int,
)
from sstar.infer import infer


def _run_infer(args: argparse.Namespace) -> None:
    """
    Execute the infer process with specified parameters.

    Parameters
    ----------
    args : argparse.Namespace
        A namespace object obtained from argparse, containing specified parameters.
    """
    infer(
        model=args.model,
        config=args.config,
        feat_file=args.feat_file,
        pred_file=args.pred_file,
        tract_file=args.tract_file,
        match_bonus=args.match_bonus,
        max_mismatch=args.max_mismatch,
        mismatch_penalty=args.mismatch_penalty,
    )


def add_infer_parser(subparsers: argparse.ArgumentParser) -> None:
    """
    Initializes and configures the command-line interface parser
    for the infer subcommand.

    Parameters
    ----------
    subparsers : argparse.ArgumentParser
        A command-line interface parser to be configured.
    """
    parser = subparsers.add_parser(
        "infer", help="Run the infer command based on specified parameters."
    )
    parser.add_argument(
        "--model",
        type=existed_file,
        required=True,
        help="Path to the model file.",
    )
    parser.add_argument(
        "--config",
        type=existed_file,
        required=True,
        help="Path to the config file.",
    )
    parser.add_argument(
        "--feat-file", required=True, help="Path to the feature TSV file."
    )
    parser.add_argument(
        "--pred-file", required=True, help="Path to the prediction TSV file."
    )
    parser.add_argument(
        "--tract-file", required=True, help="Path to the output BED file."
    )
    parser.add_argument(
        "--match-bonus", type=positive_int, default=5000, help="S* match bonus."
    )
    parser.add_argument(
        "--max-mismatch",
        type=non_negative_int,
        default=5,
        help="S* maximum mismatches.",
    )
    parser.add_argument(
        "--mismatch-penalty",
        type=non_positive_int,
        default=-10000,
        help="S* mismatch penalty.",
    )
    parser.set_defaults(runner=_run_infer)
