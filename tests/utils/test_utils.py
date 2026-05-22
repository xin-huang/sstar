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

import pytest
import numpy as np
from sstar.utils import split_genome


@pytest.fixture
def pos():
    return np.array(
        [261, 592, 619, 693, 1008, 1125, 1307, 1424, 1436, 1491, 1743, 1853, 1917, 2075]
    )


def test_split_genome_windows(pos):
    windows = split_genome(
        pos=pos, chr_name="1", polymorphism_size=50000, step_size=10000
    )

    assert windows == [("1", [1, 50000])]


def test_split_genome_windows_starts_at_one():
    pos = np.array([1, 2, 3, 4, 5])
    windows = split_genome(pos=pos, chr_name="1", polymorphism_size=3, step_size=2)

    assert windows[0] == ("1", [1, 3])


def test_split_genome_polymorphisms(pos):
    polymorphisms = split_genome(
        pos=pos, chr_name="1", polymorphism_size=10, step_size=1, window_based=False
    )

    expected_polymorphisms = [
        ("1", [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]),
        ("1", [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]),
        ("1", [2, 3, 4, 5, 6, 7, 8, 9, 10, 11]),
        ("1", [3, 4, 5, 6, 7, 8, 9, 10, 11, 12]),
        ("1", [4, 5, 6, 7, 8, 9, 10, 11, 12, 13]),
    ]

    for actual, expected in zip(polymorphisms, expected_polymorphisms):
        assert actual[0] == expected[0]
        assert np.array_equal(actual[1], expected[1])


def test_split_genome_random_polymorphisms(pos):
    random_polymorphisms = split_genome(
        pos=pos,
        chr_name="1",
        polymorphism_size=10,
        window_based=False,
        random_polymorphisms=True,
        seed=12345,
    )

    assert random_polymorphisms == [("1", [0, 3, 6, 7, 8, 9, 10, 11, 12, 13])]


def test_split_genome_negative_step_size(pos):
    with pytest.raises(ValueError) as exc_info:
        split_genome(
            pos=pos,
            chr_name="1",
            polymorphism_size=192,
            step_size=-100,
            window_based=False,
        )

    assert (
        str(exc_info.value)
        == "`step_size` and `polymorphism_size` must be positive integers."
    )


def test_split_genome_negative_polymorphism_size(pos):
    with pytest.raises(ValueError) as exc_info:
        split_genome(
            pos=pos,
            chr_name="1",
            polymorphism_size=-192,
            step_size=100,
            window_based=False,
        )

    assert (
        str(exc_info.value)
        == "`step_size` and `polymorphism_size` must be positive integers."
    )


def test_split_genome_empty_positions(pos):
    with pytest.raises(ValueError) as exc_info:
        split_genome(
            pos=[],
            chr_name="1",
            polymorphism_size=192,
            step_size=100,
            window_based=False,
        )

    assert str(exc_info.value) == "`pos` array must not be empty."


def test_split_genome_no_windows_created(pos):
    with pytest.raises(ValueError) as exc_info:
        split_genome(
            pos=pos,
            chr_name="1",
            polymorphism_size=192,
            step_size=100,
            window_based=False,
        )

    assert (
        str(exc_info.value)
        == "No windows could be created with the given number of polymorphisms and step size."
    )


def test_split_genome_no_windows_created_random_polymorphisms(pos):
    with pytest.raises(ValueError) as exc_info:
        split_genome(
            pos=pos,
            chr_name="1",
            polymorphism_size=192,
            window_based=False,
            random_polymorphisms=True,
        )

    assert (
        str(exc_info.value)
        == "No windows could be created with the given number of polymorphisms."
    )
