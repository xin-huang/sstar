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

import argparse, os


def positive_int(value: str) -> int:
    """
    Validates if the provided string represents a positive integer.

    Parameters
    ----------
    value : str
        The value to validate.

    Returns
    -------
    int
        The validated positive integer.

    Raises
    ------
    argparse.ArgumentTypeError
        If the value is not a valid integer or positive integer.

    """
    if value is not None:
        try:
            ivalue = int(value)
        except ValueError:
            raise argparse.ArgumentTypeError(f"{value} is not a valid integer")
        if ivalue <= 0:
            raise argparse.ArgumentTypeError(f"{value} is not a positive integer")
    return ivalue


def positive_number(value: str) -> float:
    """
    Validates if the provided string represents a positive number.

    Parameters
    ----------
    value : str
        The value to validate.

    Returns
    -------
    float
        The validated positive number.

    Raises
    ------
    argparse.ArgumentTypeError
        If the value is not a valid number or positive number.

    """
    if value is not None:
        try:
            fvalue = float(value)
        except ValueError:
            raise argparse.ArgumentTypeError(f"{value} is not a valid number")
        if fvalue <= 0:
            raise argparse.ArgumentTypeError(f"{value} is not a positive number")
    return fvalue


def existed_file(value: str) -> str:
    """
    Validates if the provided string is a path to an existing file.

    Parameters
    ----------
    value : str
        The path to validate.

    Returns
    -------
    str
        The validated file path.

    Raises
    ------
    argparse.ArgumentTypeError
        If the file does not exist.

    """
    if value is not None:
        if not os.path.isfile(value):
            raise argparse.ArgumentTypeError(f"{value} is not found")
    return value
