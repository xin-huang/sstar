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
import pandas as pd
from abc import ABC, abstractmethod


class MlModel(ABC):
    """
    Abstract base class for machine learning models.

    This class defines a standard interface for training and inferring from
    machine learning models. It ensures that all derived model classes implement
    the `train` and `infer` methods.
    """

    @abstractmethod
    def train(self, **kwargs):
        """
        Train the machine learning model.

        Parameters
        ----------
        **kwargs : dict
            Additional keyword arguments specific to the model being used.
        """
        pass

    @abstractmethod
    def infer(self, **kwargs):
        """
        Perform inference using the trained model on the provided data.

        Parameters
        ----------
        **kwargs : dict
            Additional keyword arguments specific to the model being used.
        """
        pass
