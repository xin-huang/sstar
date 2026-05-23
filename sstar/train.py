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
from sstar.simulate import simulate
from sstar.utils import UniqueKeyLoader


def train(
    demes: str,
    config: str,
    output: str | None,
    only_simulation: bool = False,
) -> None:
    """
    Run simulation and model training from YAML configuration.

    Parameters
    ----------
    demes : str
        Path to the demography (demes) YAML file used for simulation.
    config : str
        Path to the sstar configuration YAML file.
    output : str, optional
        Output path or directory passed to the model's `train` method, used
        to store the trained model. Required when `only_simulation=False`.
        Default: None.
    only_simulation : bool, optional
        If True, run simulation checks/simulation only and skip model
        training. Default: False.
    """
    try:
        with open(config, "r") as f:
            config_dict = yaml.load(f, Loader=UniqueKeyLoader)
    except FileNotFoundError:
        raise FileNotFoundError(f"Configuration file '{config}' not found.")
    except yaml.YAMLError as e:
        raise ValueError(f"Error parsing YAML configuration file '{config}': {e}")

    global_config = GlobalConfig(**config_dict)

    output_dir = os.path.dirname(str(output))
    output_prefix = os.path.splitext(os.path.basename(str(output)))[0]
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    data = os.path.join(output_dir, f"{output_prefix}.training.features.tsv")

    if not os.path.exists(data):
        print("Training data is not found. Performing simulation ...")
        simulate(
            demo_model_file=demes,
            output_dir=output_dir,
            output_prefix=output_prefix,
            **global_config.simulation.model_dump(),
        )

    if only_simulation:
        return

    qr.train(
        data=data,
        output=output,
        **global_config.model.params,
    )
