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


import numpy as np
from sstar.mp_manager import mp_manager
from sstar.generators import RandomNumberGenerator


class SimpleJob:
    def __init__(self, item):
        self.item = item

    def run(self, rep, seed):
        return [(rep, seed, self.item)]


class FailureJob:
    def __init__(self, item):
        self.item = item

    def run(self, rep, seed):
        raise Exception("Simulating failure by stopping.")


def test_mp_manager():
    nprocess = 2
    nrep = 5
    seed = 2

    generator = RandomNumberGenerator(nrep=nrep)

    job = SimpleJob("Hello")

    results = mp_manager(job=job, data_generator=generator, nprocess=nprocess)
    expected_results = [(i, None, "Hello") for i in range(nrep)]
    assert sorted(results) == sorted(
        expected_results
    ), "mp_manager did not return the expected results"

    generator = RandomNumberGenerator(nrep=nrep, seed=seed)
    results = mp_manager(job=job, data_generator=generator, nprocess=nprocess)
    np.random.seed(seed)
    seed_list = np.random.randint(1, 2**31, nrep)
    expected_results = [(i, seed_list[i], "Hello") for i in range(nrep)]
    assert sorted(results) == sorted(
        expected_results
    ), "mp_manager did not return the expected results"

    start_rep = 10
    generator = RandomNumberGenerator(start_rep=start_rep, nrep=nrep, seed=seed)
    results = mp_manager(job=job, data_generator=generator, nprocess=nprocess)
    np.random.seed(seed)
    seed_list = np.random.randint(1, 2**31, start_rep + nrep)
    expected_results = [
        (i, seed_list[i], "Hello") for i in range(start_rep, start_rep + nrep)
    ]
    assert sorted(results) == sorted(
        expected_results
    ), "mp_manager did not return the expected results"

    generator = RandomNumberGenerator(nrep=1, seed=seed)
    results = mp_manager(job=job, data_generator=generator, nprocess=nprocess)
    expected_results = [(0, 2, "Hello")]
    assert sorted(results) == sorted(
        expected_results
    ), "mp_manager did not return the expected results"


def test_mp_manager_failure(capfd):
    nprocess = 2
    generator = RandomNumberGenerator(nrep=5)

    job = FailureJob("Hello")
    mp_manager(job=job, data_generator=generator, nprocess=nprocess)

    # Use capfd to capture stdout and stderr
    captured = capfd.readouterr()

    # Assertions to verify expected output and behavior
    assert "Simulating failure by stopping." in captured.err
    assert "did not complete successfully. Initiating shutdown." in captured.out
    assert "All workers are terminated." in captured.out
