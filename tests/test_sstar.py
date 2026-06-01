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
from sstar.sstar import Sstar


def test_sstar_all_zero_returns_zero():
    ref = np.zeros((3, 2), dtype=int)
    tgt = np.zeros((3, 2), dtype=int)  # (n_sites=3, n_samples=2)
    pos = np.array([10, 50, 100], dtype=int)

    out = Sstar.compute(ref_gts=ref, tgt_gts=tgt, pos=pos)["sstar"]
    assert out == [np.nan, np.nan]


def test_sstar_two_matching_loci_far_apart_gives_bonus_plus_distance():
    ref = np.zeros((3, 2), dtype=int)
    tgt = np.array(
        [
            [1, 1],  # site 0
            [2, 0],  # site 1
            [1, 1],  # site 2
        ],
        dtype=int,
    )  # shape (n_sites=2, n_samples=2)

    pos = np.array([0, 50, 100], dtype=int)
    match_bonus = 5000
    expected = 100 + match_bonus

    out = Sstar.compute(ref_gts=ref, tgt_gts=tgt, pos=pos, match_bonus=match_bonus)[
        "sstar"
    ]
    assert len(out) == 2
    assert out[0] == expected
    assert np.isnan(out[1])


def test_sstar_close_positions_or_too_few_loci_yield_nan():
    ref = np.zeros((3, 2), dtype=int)

    tgt_close = np.array(
        [
            [1, 1],
            [1, 1],
            [1, 1],
        ],
        dtype=int,
    )
    pos_close = np.array([0, 5, 8], dtype=int)  # all pair distances < 10

    out_close = Sstar.compute(
        ref_gts=ref,
        tgt_gts=tgt_close,
        pos=pos_close,
    )["sstar"]

    assert len(out_close) == 2
    assert np.all(np.isnan(out_close))

    ref = np.zeros((2, 2), dtype=int)

    tgt_too_few = np.array(
        [
            [1, 0],
            [0, 1],
        ],
        dtype=int,
    )
    pos_any = np.array([0, 100], dtype=int)

    out_too_few = Sstar.compute(
        ref_gts=ref,
        tgt_gts=tgt_too_few,
        pos=pos_any,
    )["sstar"]

    assert len(out_too_few) == 2
    assert np.all(np.isnan(out_too_few))


def test_sstar_missing_params():
    params = {}

    with pytest.raises(ValueError):
        Sstar.compute(**params)
