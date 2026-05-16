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

import inspect
from typing import Any, Literal, Optional

from pydantic import BaseModel, ConfigDict, Field
from pydantic import model_validator
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.linear_model import LogisticRegression


class LrParams(BaseModel):
    """
    Parameter schema for logistic regression model configuration.
    """

    model_config = ConfigDict(extra="allow")

    is_scaled: bool = False

    @model_validator(mode="after")
    def validate_known_keys(self) -> "LrParams":
        """
        Validate keys against LogisticRegression signature plus gaishi-specific keys.
        """
        allowed_keys = set(inspect.signature(LogisticRegression).parameters)
        allowed_keys.add("is_scaled")
        unknown_keys = set(self.model_extra or {}).difference(allowed_keys)
        if unknown_keys:
            unknown_str = ", ".join(sorted(unknown_keys))
            raise ValueError(f"Unknown logistic_regression params: {unknown_str}")
        return self


class EtcParams(BaseModel):
    """
    Parameter schema for extra-trees classifier model configuration.
    """

    model_config = ConfigDict(extra="allow")

    is_scaled: bool = False

    @model_validator(mode="after")
    def validate_known_keys(self) -> "EtcParams":
        """
        Validate keys against ExtraTreesClassifier signature plus gaishi-specific keys.
        """
        allowed_keys = set(inspect.signature(ExtraTreesClassifier).parameters)
        allowed_keys.add("is_scaled")
        unknown_keys = set(self.model_extra or {}).difference(allowed_keys)
        if unknown_keys:
            unknown_str = ", ".join(sorted(unknown_keys))
            raise ValueError(f"Unknown extra_trees_classifier params: {unknown_str}")
        return self


class UnetParams(BaseModel):
    """
    Strict parameter schema for the UNet++ model.
    """

    model_config = ConfigDict(extra="forbid")

    trained_model_file: Optional[str] = None
    add_rnn: bool = False
    site_weighting: bool = False
    batch_size: int = Field(default=32, gt=0)
    n_early: int = Field(default=10, gt=0)
    n_epochs: int = Field(default=100, gt=0)
    learning_rate: float = Field(default=0.001, gt=0)
    min_delta: float = Field(default=1e-4, gt=0)
    val_prop: float = Field(default=0.05, gt=0, lt=1)
    seed: Optional[int] = None
    device: Optional[str] = None
    num_workers: int = Field(default=0, ge=0)
    train_drop_last: Optional[bool] = None
    val_drop_last: Optional[bool] = None
    persistent_workers: Optional[bool] = None
    prefetch_factor: Optional[int] = Field(default=None, gt=0)
    use_amp: bool = False
    recent_window: int = Field(default=500, ge=1, le=1000)


class ModelConfig(BaseModel):
    name: Literal["logistic_regression", "extra_trees_classifier", "unet++"]
    params: dict[str, Any] = Field(default_factory=dict)

    model_config = ConfigDict(extra="forbid")

    @model_validator(mode="after")
    def validate_params_by_model_name(self) -> "ModelConfig":
        """
        Validate model params using strict per-model schemas.
        """
        schema_map = {
            "logistic_regression": LrParams,
            "extra_trees_classifier": EtcParams,
            "unet++": UnetParams,
        }
        params_model = schema_map[self.name](**self.params)
        self.params = params_model.model_dump()
        extra = params_model.model_extra or {}
        self.params.update(extra)
        return self
