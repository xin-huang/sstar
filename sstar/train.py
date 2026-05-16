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
from typing import Optional
from sstar.configs import GlobalConfig
from sstar.registries.model_registry import MODEL_REGISTRY
from sstar.simulate import simulate
from sstar.utils import UniqueKeyLoader, filter_model_params_for_method


def train(
    demes: str,
    config: str,
    output: Optional[str] = None,
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
    data = f"{global_config.simulation.output_dir}/{global_config.simulation.output_prefix}.tsv"
    if not os.path.exists(data):
        print("Training data is not found. Perform simulation.")
        simulation_params = filter_model_params_for_method(
            simulate, global_config.simulation.model_dump()
        )
        simulate(
            demo_model_file=demes,
            **simulation_params,
        )

    if only_simulation:
        return
    if output is None:
        raise ValueError("`output` is required unless `only_simulation=True`.")

    model_name = global_config.model.name
    model_params = global_config.model.params
    model_cls = MODEL_REGISTRY.get(model_name)
    model_params = filter_model_params_for_method(model_cls.train, model_params)
    model_cls.train(
        data=data,
        output=output,
        **model_params,
    )
