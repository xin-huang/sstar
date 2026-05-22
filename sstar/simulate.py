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

import os, random, shutil

import pandas as pd
from sstar import mp_manager
from sstar.generators import RandomNumberGenerator
from sstar.msprime_simulator import MsprimeSimulator


def simulate(
    demo_model_file: str,
    nrep: int,
    nref: int,
    ntgt: int,
    ref_id: str,
    tgt_id: str,
    src_id: str,
    ploidy: int,
    seq_len: int,
    mut_rate: float,
    rec_rate: float,
    feature_config_file: str,
    output_prefix: str,
    output_dir: str,
    nfeature: int,
    is_phased: bool,
    is_shuffled: bool,
    keep_sim_data: bool,
    seed: int,
    nprocess: int,
) -> None:
    """
    Simulate genomic data and generate feature vectors for training.

    This function orchestrates simulation and S* feature-vector generation for training.

    Parameters
    ----------
    demo_model_file : str
        Path to the demographic model file for simulations.
    nrep : int
        Number of replicates to simulate.
    nref : int
        Number of reference individuals in the simulation.
    ntgt : int
        Number of target individuals in the simulation.
    ref_id : str
        Identifier for the reference population.
    tgt_id : str
        Identifier for the target population.
    src_id : str
        Identifier for the source population.
    ploidy : int
        Ploidy of the individuals.
    seq_len : int
        Length of the sequence to simulate.
    mut_rate : float
        Mutation rate per base pair per generation.
    rec_rate : float
        Recombination rate per base pair per generation.
    feature_config_file : str
        Path to the YAML configuration file specifying features to compute.
    output_prefix : str
        Prefix for output files.
    output_dir : str
        Directory to save output files.
    nfeature : int
        Total number of features to generate.
    is_phased : bool
        Flag indicating if the data is phased.
    is_shuffled: bool
        If True, shuffles the feature vectors before output.
    keep_sim_data : bool
        If False, deletes simulation data directories after processing to save disk space.
    seed : int
        Seed for random number generation to ensure reproducibility.
    nprocess : int
        Number of processes to use for parallel execution.

    Raises
    ------
    ValueError
        If `nfeature` is less than or equal to zero or other invalid parameters are provided.
    SystemExit
        If errors occur in mp_manager.

    Notes
    -----
    The simulation process is dependent on the correctness of the `LRTrainingDataSimulator`,
    `RandomNumberGenerator`, and the feature computation classes. Ensure these components are correctly
    implemented and tested.

    """
    if nfeature <= 0:
        raise ValueError("nfeature must be positive.")

    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"{output_prefix}.tsv")

    simulator = MsprimeSimulator(
        demo_model_file=demo_model_file,
        nref=nref,
        ntgt=ntgt,
        ref_id=ref_id,
        tgt_id=tgt_id,
        src_id=src_id,
        ploidy=ploidy,
        seq_len=seq_len,
        mut_rate=mut_rate,
        rec_rate=rec_rate,
        output_prefix=output_prefix,
        output_dir=output_dir,
        is_phased=is_phased,
        feature_config_file=feature_config_file,
    )

    total_features = []
    start_rep = 0

    is_stopped = False

    while not is_stopped:
        generator = RandomNumberGenerator(start_rep=start_rep, nrep=nrep, seed=seed)
        features = mp_manager(
            job=simulator, data_generator=generator, nprocess=nprocess
        )

        if features == "error":
            raise SystemExit("Some errors occurred, stopping the program ...")

        features.sort(
            key=lambda x: (
                x["Replicate"],
                x["Chromosome"],
                x["Start"],
                x["End"],
                x["Sample"],
            )
        )

        if not keep_sim_data:
            for i in range(start_rep, start_rep + nrep):
                shutil.rmtree(os.path.join(output_dir, f"{i}"), ignore_errors=True)

        total_features.extend(features)
        is_stopped = len(total_features) >= nfeature

        print(f"Number of obtained features: {len(total_features)}")

        start_rep += nrep

    total_features = total_features[:nfeature]

    if is_shuffled:
        random.seed(seed)
        random.shuffle(total_features)
    else:
        total_features.sort(
            key=lambda x: (
                x["Replicate"],
                x["Chromosome"],
                x["Start"],
                x["End"],
                x["Sample"],
            )
        )

    pd.DataFrame(total_features).to_csv(output_file, sep="\t", index=False)
