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


import os
import yaml
import sstar.quantile_regression as qr
from sstar.configs import GlobalConfig
from sstar.preprocess import preprocess
from sstar.utils import UniqueKeyLoader


def infer(
    model: str,
    config: str,
    feat_file: str,
    pred_file: str,
    tract_file: str,
) -> None:
    """
    Run feature preprocessing and model-based inference from a YAML configuration.

    Parameters
    ----------
    model : str
        Path to the trained model file to be used for inference.
    config : str
        Path to the inference configuration YAML file.
    output : str
        Path where inference outputs.
    """
    try:
        with open(config, "r") as f:
            config_dict = yaml.load(f, Loader=UniqueKeyLoader)
    except FileNotFoundError:
        raise FileNotFoundError(f"Configuration file '{config}' not found.")
    except yaml.YAMLError as e:
        raise ValueError(f"Error parsing YAML configuration file '{config}': {e}")

    global_config = GlobalConfig(**config_dict)

    preprocess(
        output_file=feat_file,
        **global_config.preprocessing.model_dump(),
    )

    qr.infer(
        data=feat_file,
        model=model,
        output=pred_file,
        bed_file=tract_file,
    )
