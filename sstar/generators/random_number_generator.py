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
from gaishi.generators import GenericGenerator


class RandomNumberGenerator(GenericGenerator):
    """
    Generates random numbers based on specified parameters.

    """

    def __init__(self, nrep: int, start_rep: int = 0, seed: int = None):
        """
        Initializes a new instance of RandomNumberGenerator.

        Parameters
        ----------
        nrep : int
            The total number of repetitions to generate random numbers for.
        start_rep : int, optional
            The starting replicate number. Default is 0.
        seed : int, optional
            A global seed for the random number generator to ensure reproducibility across
            multiple runs. If not provided, random seeds will be generated without a fixed seed.
            Default is None.

        Raises
        ------
        ValueError
            If `nrep` is less than or equal to 0, or `start_rep` is negative.

        """
        if nrep <= 0:
            raise ValueError("nrep must be greater than 0.")

        if start_rep < 0:
            raise ValueError("start_rep must not be negative.")

        end_rep = start_rep + nrep
        self.rep_list = list(range(start_rep, end_rep))

        if seed is None:
            self.seed_list = [None] * len(self.rep_list)
        elif end_rep > 1:
            np.random.seed(seed)
            self.seed_list = np.random.randint(1, 2**31, end_rep).tolist()
            self.seed_list = self.seed_list[start_rep:end_rep]
        else:
            self.seed_list = [seed]

    def get(self):
        """
        Yields dictionaries containing 'rep' and 'seed' values for each repetition.

        Yields
        ------
        dict
            A dictionary with two keys: 'rep' and 'seed'. 'rep' corresponds to the replicate number,
            and 'seed' corresponds to the random seed for that replicate. Each call to `get()`
            yields a new dictionary for the next replicate until all replicates are exhausted.

        """
        for rep, seed in zip(self.rep_list, self.seed_list):
            yield {"rep": rep, "seed": seed}

    def __len__(self):
        return len(self.rep_list)
