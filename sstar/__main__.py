# Copyright 2025 Xin Huang
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
import os
import signal
import sys
from typing import Optional, List

# Import all parser modules
from sstar.parsers.score_parser import add_score_parser
from sstar.parsers.quantile_parser import add_quantile_parser
from sstar.parsers.threshold_parser import add_threshold_parser
from sstar.parsers.match_rate_parser import add_match_rate_parser
from sstar.parsers.tract_parser import add_tract_parser

def set_sigpipe_handler() -> None:
    """
    Sets the signal handler for SIGPIPE signals on POSIX systems.
    """
    if os.name == "posix":
        # Set signal handler for SIGPIPE to quietly kill the program.
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)

def sstar_cli_parser() -> argparse.ArgumentParser:
    """
    Initializes and configures the command-line interface parser for sstar.
    
    Returns
    -------
    top_parser : argparse.ArgumentParser
        A configured command-line interface parser.
    """
    top_parser = argparse.ArgumentParser()
    subparsers = top_parser.add_subparsers(dest="subcommand")
    subparsers.required = True
    
    # Add each subcommand parser
    add_score_parser(subparsers)
    add_quantile_parser(subparsers)
    add_threshold_parser(subparsers)
    add_match_rate_parser(subparsers)
    add_tract_parser(subparsers)
    
    return top_parser

def main(arg_list: Optional[List[str]] = None) -> None:
    """
    Main entry for sstar.

    Parameters
    ----------
    arg_list : list of str, optional
        List containing arguments for sstar. If None, defaults to sys.argv[1:].

    Returns
    -------
    None
    """
    set_sigpipe_handler()
    parser = sstar_cli_parser()
    args = parser.parse_args(arg_list)
    args.runner(args)

if __name__ == "__main__":
    main()
