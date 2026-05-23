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
from typing import Optional
from pydantic import BaseModel, ConfigDict, Field, field_validator


class PreprocessingConfig(BaseModel):
    """Configuration for preprocessing data."""

    model_config = ConfigDict(extra="forbid")

    # Input
    vcf_file: Path = Field(..., description="Path to input VCF/BCF file")
    chr_name: str = Field(..., description="Chromosome name to process")
    ref_ind_file: Path = Field(..., description="Path to reference individual ID list")
    tgt_ind_file: Path = Field(..., description="Path to target individual ID list")

    # Output
    output_file: Path = Field(..., description="Path to output feature file")

    # Parallelization / data properties
    nprocess: int = Field(
        1,
        gt=0,
        description="Number of processes to use for preprocessing",
    )
    ploidy: int = Field(2, gt=0, description="Ploidy of samples")
    is_phased: bool = Field(
        True,
        description="Whether the VCF genotypes are phased",
    )

    # Optional ancestral allele information
    anc_allele_file: Optional[Path] = Field(
        None,
        description="Optional path to ancestral allele file (if available)",
    )

    # Windowing
    win_len: int = Field(..., gt=0, description="Window length in bp")
    win_step: int = Field(..., gt=0, description="Window step size in bp")

    # Features
    feature_config_file: Path = Field(
        ..., description="Path to feature configuration YAML/JSON"
    )

    @field_validator("output_file")
    @classmethod
    def _ensure_output_file(cls, p: Path) -> Path:
        """Normalize output_file to an absolute path without creating parent dir."""
        return p.expanduser().resolve()
