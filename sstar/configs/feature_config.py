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

from typing import Dict, Union, Set
from pydantic import RootModel, field_validator

SUPPORTED_FEATURES: Set[str] = {
    "ref_dist",
    "tgt_dist",
    "spectrum",
    "private_mutation",
    "sstar",
}

SUPPORTED_DIST_STATS: Set[str] = {
    "all",
    "minimum",
    "maximum",
    "mean",
    "median",
    "variance",
    "skew",
    "kurtosis",
}

SUPPORTED_SSTAR_PARAMS: Set[str] = {
    "match_bonus",
    "max_mismatch",
    "mismatch_penalty",
}


class FeatureConfig(
    RootModel[
        Dict[
            str,
            Union[bool, Dict[str, Union[bool, int]]],
        ]
    ]
):
    """
    Schema and validators for feature configuration.

    This model validates a dictionary-based feature configuration used to
    control which features and sub-options are enabled for downstream
    computations. It enforces allowed feature names and constrains the
    types/keys of their parameters.

    Structure
    ---------
    - `ref_dist` / `tgt_dist` : dict
        Mapping from distance statistic names to booleans. Allowed stats
        are in `SUPPORTED_DIST_STATS`. At least one must be `True`.
    - `spectrum` / `num_private` : bool
        Feature toggles.
    - `sstar` : dict
        Mapping from S* parameter names to integer values. Allowed params
        are in `SUPPORTED_SSTAR_PARAMS`.

    Notes
    -----
    The model is defined as a `RootModel` where the underlying `root`
    is a `Dict[str, Union[bool, Dict[str, Union[bool, int]]]]` to allow
    both boolean toggles and nested option dictionaries at the top level.
    """

    @field_validator("root")
    def check_valid_feature_types(
        cls, v: Dict[str, Union[bool, Dict[str, Union[bool, int]]]]
    ):
        """
        Validates top-level features and dispatch to per-feature checks.

        Parameters
        ----------
        v : dict
            The parsed feature configuration mapping feature names to either
            booleans (for simple toggles) or dictionaries (for nested options).

        Returns
        -------
        dict
            The validated configuration (unchanged).

        Raises
        ------
        ValueError
            If a feature name is unsupported, if a feature has an unexpected
            type (e.g., boolean where a dict is required), or if a feature's
            options fail their specific validation rules.
        """
        for feat_name, params in v.items():
            if feat_name not in SUPPORTED_FEATURES:
                raise ValueError(
                    f"Unsupported feature: {feat_name!r}. "
                    f"Allowed: {sorted(SUPPORTED_FEATURES)}"
                )

            if feat_name in {"ref_dist", "tgt_dist"}:
                if not isinstance(params, dict):
                    raise ValueError(
                        f"{feat_name} must be a mapping of stats to bools."
                    )
                cls.valid_dist(feat_name, params)

            elif feat_name in {"spectrum", "private_mutation"}:
                if not isinstance(params, bool):
                    raise ValueError(f"{feat_name} must be a boolean.")

            elif feat_name == "sstar":
                if not isinstance(params, dict):
                    raise ValueError("sstar must be a mapping of params to integers.")
                cls.valid_sstar(feat_name, params)

        return v

    @staticmethod
    def valid_dist(feat_name: str, params: Dict[str, Union[bool, int]]):
        """
        Validates distance-statistics options for `ref_dist`/`tgt_dist`.

        Parameters
        ----------
        feat_name : str
            The feature name being validated (`"ref_dist"` or `"tgt_dist"`).
        params : dict
            Mapping from statistic names to booleans.

        Raises
        ------
        ValueError
            If unknown statistic keys are present, if any value is not a bool,
            if no statistic is enabled, or if `all=True` conflicts with other
            enabled statistics.
        """
        keys = set(params.keys())
        if unknown := keys - SUPPORTED_DIST_STATS:
            raise ValueError(
                f"{feat_name}: unknown dist stats {sorted(unknown)}. "
                f"Allowed: {sorted(SUPPORTED_DIST_STATS)}"
            )

        for k, v in params.items():
            if not isinstance(v, bool):
                raise ValueError(
                    f"{feat_name}: value for '{k}' must be bool, got {type(v).__name__}."
                )

        if not any(bool(x) for x in params.values()):
            raise ValueError(f"{feat_name}: at least one stat must be True.")

    @staticmethod
    def valid_sstar(feat_name: str, params: Dict[str, Union[bool, int]]):
        """
        Validates S* (sstar) parameter options.

        Parameters
        ----------
        feat_name : str
            The feature name being validated (`"sstar"`).
        params : dict
            Mapping from S* parameter names to integer values.

        Raises
        ------
        ValueError
            If unknown parameter keys are present, if a value is not an int,
            or if numeric constraints are violated (e.g., negative `max mismatch`,
            positive `mismatch penalty`).
        """
        keys = set(params.keys())
        if unknown := keys - SUPPORTED_SSTAR_PARAMS:
            raise ValueError(
                f"{feat_name}: unknown params {sorted(unknown)}. "
                f"Allowed: {sorted(SUPPORTED_SSTAR_PARAMS)}"
            )

        for k, v in params.items():
            if type(v) is not int:
                raise ValueError(
                    f"{feat_name}: value for '{k}' must be int, got {type(v).__name__}."
                )

        if "match_bonus" in params and params["match_bonus"] <= 0:
            raise ValueError(f"{feat_name}: 'match_bonus' must be > 0.")
        if "max_mismatch" in params and params["max_mismatch"] < 0:
            raise ValueError(f"{feat_name}: 'max_mismatch' must be >= 0.")
        if "mismatch_penalty" in params and params["mismatch_penalty"] > 0:
            raise ValueError(f"{feat_name}: 'mismatch_penalty' should be <= 0.")
