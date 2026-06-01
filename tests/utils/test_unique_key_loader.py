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

import pytest
import yaml
from sstar.utils import UniqueKeyLoader


def test_unique_key_loader_single_mapping_no_duplicates():
    yaml_str = """
    force_balanced: true
    nrep: 10
    label: "test"
    """
    data = yaml.load(yaml_str, Loader=UniqueKeyLoader)

    assert isinstance(data, dict)
    assert data["force_balanced"] is True
    assert data["nrep"] == 10
    assert data["label"] == "test"


def test_unique_key_loader_nested_mapping_no_duplicates():
    yaml_str = """
    simulation:
      nrep: 10
      nref: 20
      force_balanced: false
    model:
      name: "extra_trees_classifier"
      params:
        n_estimators: 500
        n_jobs: -1
    """
    data = yaml.load(yaml_str, Loader=UniqueKeyLoader)

    assert isinstance(data, dict)
    assert data["simulation"]["nrep"] == 10
    assert data["simulation"]["force_balanced"] is False
    assert data["model"]["name"] == "extra_trees_classifier"
    assert data["model"]["params"]["n_estimators"] == 500


def test_unique_key_loader_raises_on_duplicate_top_level_key():
    yaml_str = """
    force_balanced: false
    force_balanced: true
    """
    with pytest.raises(ValueError) as excinfo:
        yaml.load(yaml_str, Loader=UniqueKeyLoader)

    msg = str(excinfo.value)
    assert "Duplicate key in YAML" in msg
    assert "'force_balanced'" in msg


def test_unique_key_loader_raises_on_duplicate_nested_key():
    yaml_str = """
    simulation:
      nrep: 10
      nrep: 20
      force_balanced: true
    """
    with pytest.raises(ValueError) as excinfo:
        yaml.load(yaml_str, Loader=UniqueKeyLoader)

    msg = str(excinfo.value)
    assert "Duplicate key in YAML" in msg
    assert "'nrep'" in msg
