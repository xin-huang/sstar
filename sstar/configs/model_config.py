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
from typing import Any

from pydantic import BaseModel, ConfigDict, Field
from pydantic import model_validator
from sklearn.ensemble import GradientBoostingRegressor


class ModelConfig(BaseModel):
    params: dict[str, Any] = Field(default_factory=dict)

    model_config = ConfigDict(extra="forbid")

    @model_validator(mode="after")
    def validate_known_keys(self):
        """
        Validate keys against GradientBoostingRegressor signature plus sstar-specific keys.
        """
        allowed_keys = set(inspect.signature(GradientBoostingRegressor).parameters)
        unknown_keys = set(self.model_extra or {}).difference(allowed_keys)
        if unknown_keys:
            unknown_str = ", ".join(sorted(unknown_keys))
            raise ValueError(f"Unknown GradientBoostingRegressor params: {unknown_str}")
        return self
