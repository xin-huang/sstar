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


class GenericGenerator(ABC):
    """
    Abstract base class for generating data.

    This class defines a common interface for data generation. Subclasses
    must implement the get method to generate data according to specific
    requirements or configurations provided via keyword arguments.

    """

    @abstractmethod
    def get(self, **kwargs):
        """
        Generates data based on the provided keyword arguments.

        Subclasses should implement this method to generate and return data
        according to the requirements described by the keyword arguments.

        Parameters:
        **kwargs: Arbitrary keyword arguments specific to the data generation
        implementation in subclasses.

        Returns:
        The generated data, the format and type of which are determined by the
        subclass implementation.

        """
        pass
