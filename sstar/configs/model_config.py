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
from collections.abc import Mapping
from typing import Any

from pydantic import BaseModel, ConfigDict, Field
from pydantic import field_validator, model_validator
from sklearn.ensemble import GradientBoostingRegressor


class ModelConfig(BaseModel):
    params: dict[str, Any] = Field(default_factory=dict)

    model_config = ConfigDict(extra="forbid")

    @field_validator("params", mode="before")
    @classmethod
    def validate_params_mapping(cls, value: Any) -> dict[str, Any]:
        """
        Validate that `params` is provided as a mapping.
        """
        if value is None:
            return {}
        if not isinstance(value, Mapping):
            raise TypeError("`params` must be a mapping (e.g., dict[str, Any]).")
        return dict(value)

    @model_validator(mode="after")
    def validate_known_keys(self):
        """
        Validate `params` against the GradientBoostingRegressor signature.
        """
        allowed_keys = set(inspect.signature(GradientBoostingRegressor).parameters)
        unknown_keys = set(self.params).difference(allowed_keys)

        if unknown_keys:
            unknown_str = ", ".join(sorted(unknown_keys))
            raise ValueError(f"Unknown GradientBoostingRegressor params: {unknown_str}")

        return self
