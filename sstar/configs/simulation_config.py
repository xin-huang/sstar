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

from pathlib import Path
from pydantic import BaseModel, ConfigDict
from pydantic import Field, field_validator


class SimulationConfig(BaseModel):
    """Configuration for simulating training data."""

    model_config = ConfigDict(extra="forbid")

    # Replicates and sample sizes
    nrep: int = Field(..., gt=0, description="Number of replicates")
    nref: int = Field(..., gt=0, description="Number of reference samples")
    ntgt: int = Field(..., gt=0, description="Number of target samples")

    # Population IDs
    ref_id: str = Field(..., description="Reference population label")
    tgt_id: str = Field(..., description="Target population label")
    src_id: str = Field(..., description="Source population label")

    # Genomic parameters
    ploidy: int = Field(..., gt=0, description="Ploidy of samples")
    is_phased: bool = Field(..., description="Whether haplotypes are phased")
    seq_len: int = Field(..., gt=0, description="Sequence length (bp)")
    mut_rate: float = Field(
        ..., gt=0.0, description="Mutation rate per bp per generation"
    )
    rec_rate: float = Field(
        ..., ge=0.0, description="Recombination rate per bp per generation"
    )

    # Parallelization
    nprocess: int = Field(1, gt=0, description="Number of processes for simulation")

    force_balanced: bool = Field(
        False,
        description="Whether to enforce class balance between intro/non-intro",
    )

    # Output
    output_prefix: str = Field(..., description="Filename prefix for outputs")
    output_dir: Path = Field(..., description="Directory for all outputs")
    keep_sim_data: bool = Field(
        False,
        description="Whether to keep raw simulation data (trees, msprime/demes outputs, etc.)",
    )

    # Randomness
    seed: int = Field(..., description="Base random seed")

    feature_config_file: Path = Field(
        ..., description="Path to feature config YAML/JSON"
    )

    is_shuffled: bool = Field(
        True,
        description="Whether to shuffle feature rows (e.g. before saving/training)",
    )

    # Features
    nfeature: int = Field(..., gt=0, description="Number of features to sample/output")

    @field_validator("output_dir")
    @classmethod
    def _ensure_output_dir(cls, p: Path) -> Path:
        # Do not create here; just normalize path
        return p.expanduser().resolve()
