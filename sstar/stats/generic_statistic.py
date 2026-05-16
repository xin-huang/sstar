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

from abc import ABC, abstractmethod
from typing import Dict, Any


class GenericStatistic(ABC):
    """
    Generic class for all statistics.

    This class provides a generic interface for implementing specific statistical measures
    from genotype matrices, typically representing different populations or samples.
    """

    @staticmethod
    @abstractmethod
    def compute(self, **kwargs) -> Dict[str, Any]:
        """
        Computes the statistic based on the input genotype data.

        Parameters
        ----------
        **kwargs : dict
            Additional keyword arguments specific to the statistic being implemented.

        Returns
        -------
        dict
            A dictionary containing the results of the statistic computation.
        """
        pass

    @staticmethod
    def require(kwargs: Dict[str, Any], *names: str) -> tuple[Any]:
        """
        Fetch required keyword-only arguments from `kwargs`.

        Parameters
        ----------
        kwargs : dict
            Source mapping (typically the `**kwargs` passed to `compute`).
        *names : str
            One or more required keys to retrieve from `kwargs`.

        Returns
        -------
        tuple of Any
            Values corresponding to `names`, in the same order.
            If a single name is provided, the return value is a 1-tuple.

        Raises
        ------
        ValueError
            If any of the required keys are missing from `kwargs`.
        """
        if missing := [n for n in names if n not in kwargs]:
            raise ValueError(f"Missing required arguments: {missing}")

        return tuple(kwargs[n] for n in names)
