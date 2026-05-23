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
from sstar.generators import RandomNumberGenerator


@pytest.fixture
def expected_params():
    return [
        {"rep": 0, "seed": 843828735},
        {"rep": 1, "seed": 914636142},
        {"rep": 2, "seed": 1228959103},
        {"rep": 3, "seed": 1840268611},
        {"rep": 4, "seed": 974319581},
        {"rep": 5, "seed": 819844195},
        {"rep": 6, "seed": 220395239},
        {"rep": 7, "seed": 941243410},
        {"rep": 8, "seed": 942612052},
        {"rep": 9, "seed": 2109339755},
    ]


def test_RandomNumberGenerator(expected_params):
    generator = RandomNumberGenerator(nrep=2)
    generated_params_list = list(generator.get())

    assert generated_params_list == [
        {"rep": 0, "seed": None},
        {"rep": 1, "seed": None},
    ], "Generated parameters do not match the expected parameters."

    generator = RandomNumberGenerator(nrep=10, seed=123)
    generated_params_list = list(generator.get())

    assert (
        generated_params_list == expected_params
    ), "Generated parameters do not match the expected parameters."

    generator = RandomNumberGenerator(nrep=1, seed=123)
    generated_params_list = list(generator.get())

    assert generated_params_list == [
        {"rep": 0, "seed": 123}
    ], "Generated parameters do not match the expected parameters."
