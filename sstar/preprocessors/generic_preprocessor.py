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


class GenericPreprocessor(ABC):
    """
    Abstract base class for preprocessing genomic data.

    This class defines a common interface for various data preprocessing operations,
    such as filtering, normalization, and transformation of genomic data.

    """

    @abstractmethod
    def run(self, **kwargs):
        """
        Abstract method to run the preprocessing operations.

        Subclasses must implement this method to perform specific preprocessing tasks
        based on the initialized parameters and any additional keyword arguments.

        Parameters:
        -----------
        **kwargs : dict
            Additional keyword arguments that may be required for specific preprocessing operations.

        """
        pass
