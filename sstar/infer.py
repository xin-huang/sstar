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

import yaml
from sstar.configs import GlobalConfig
from sstar.registries.model_registry import MODEL_REGISTRY
from sstar.preprocess import preprocess_feature_vectors
from sstar.preprocess import preprocess_genotype_matrices
from sstar.utils import UniqueKeyLoader, filter_model_params_for_method


def infer(
    model: str,
    config: str,
    output: str,
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
    if global_config.preprocess.process_type == "feature_vector":
        preprocess_feature_vectors(
            **global_config.preprocess.model_dump(),
        )

        data = f"{global_config.preprocess.output_dir}/{global_config.preprocess.output_prefix}.features"
    elif global_config.preprocess.process_type == "genotype_matrix":
        preprocess_genotype_matrices(
            **global_config.preprocess.model_dump(),
        )

        data = f"{global_config.preprocess.output_dir}/{global_config.preprocess.output_prefix}.h5"
    else:
        raise ValueError("")

    model_name = global_config.model.name
    model_params = global_config.model.params
    model_cls = MODEL_REGISTRY.get(model_name)
    model_params = filter_model_params_for_method(model_cls.infer, model_params)
    model_cls.infer(
        data=data,
        model=model,
        output=output,
        **model_params,
    )
