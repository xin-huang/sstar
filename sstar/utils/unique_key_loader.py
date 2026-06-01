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

from typing import Any, Dict

import yaml


class UniqueKeyLoader(yaml.SafeLoader):
    """YAML SafeLoader that raises an error on duplicate mapping keys."""

    pass


def no_duplicate_keys(
    loader: UniqueKeyLoader,
    node: yaml.nodes.MappingNode,
    deep: bool = False,
) -> Dict[Any, Any]:
    """
    Construct a mapping from a YAML node, rejecting duplicate keys.

    Parameters
    ----------
    loader : UniqueKeyLoader
        The YAML loader instance (based on SafeLoader).
    node : yaml.nodes.MappingNode
        The mapping node being constructed.
    deep : bool, optional
        Whether to construct nested objects deeply (passed through to the loader). Default: False.

    Returns
    -------
    dict
        A Python dictionary built from the mapping node.

    Raises
    ------
    ValueError
        If a duplicate key is encountered in the same mapping.
    """
    mapping: Dict[Any, Any] = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        if key in mapping:
            raise ValueError(f"Duplicate key in YAML: {key!r}")
        value = loader.construct_object(value_node, deep=deep)
        mapping[key] = value
    return mapping


UniqueKeyLoader.add_constructor(
    yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
    no_duplicate_keys,
)
