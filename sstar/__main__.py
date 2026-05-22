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

# import sstar.models
# import sstar.stats
from sstar import __version__
from sstar.parsers.train_parser import add_train_parser
from sstar.parsers.infer_parser import add_infer_parser


def _set_sigpipe_handler() -> None:
    """
    Sets the signal handler for SIGPIPE signals on POSIX systems.
    """
    import os, signal

    if os.name == "posix":
        # Set signal handler for SIGPIPE to quietly kill the program.
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)


def _sstar_cli_parser() -> argparse.ArgumentParser:
    """
    Initializes and configures the command-line interface parser
    for sstar2.

    Returns
    -------
    top_parser : argparse.ArgumentParser
        A configured command-line interface parser.
    """
    top_parser = argparse.ArgumentParser(
        description="sstar2: S*-based Archaic Introgression Detection with Machine Learning"
    )
    top_parser.add_argument("--version", action="version", version=f"{__version__}")
    subparsers = top_parser.add_subparsers(dest="subparsers")
    subparsers.required = True

    add_train_parser(subparsers)
    add_infer_parser(subparsers)

    return top_parser


def main(arg_list: list = None) -> None:
    """
    Main entry for sstar2.

    Parameters
    ----------
    arg_list : list, optional
        A list containing arguments for sstar2. Default: None.
    """
    _set_sigpipe_handler()
    parser = _sstar_cli_parser()
    args = parser.parse_args(arg_list)
    args.runner(args)
