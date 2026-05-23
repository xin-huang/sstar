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
from sstar.parsers.argument_validation import existed_file

# from sstar.train import train


def _run_train(args: argparse.Namespace) -> None:
    """
    Execute the train process with specified parameters.

    Parameters
    ----------
    args : argparse.Namespace
        A namespace object obtained from argparse, containing specified parameters.
    """
    pass
    # train(
    #    demes=args.demes,
    #    config=args.config,
    #    output=args.output,
    # )


def add_train_parser(subparsers: argparse.ArgumentParser) -> None:
    """
    Initializes and configures the command-line interface parser
    for the train subcommand.

    Parameters
    ----------
    subparsers : argparse.ArgumentParser
        A command-line interface parser to be configured.
    """
    parser = subparsers.add_parser(
        "train", help="Run the train command based on specified parameters."
    )
    parser.add_argument(
        "--demes",
        type=existed_file,
        required=True,
        help="Path to the demes file.",
    )
    parser.add_argument(
        "--config",
        type=existed_file,
        required=True,
        help="Path to the config file.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to the output file.",
    )
    parser.set_defaults(runner=_run_train)
