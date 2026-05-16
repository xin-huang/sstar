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
from multiprocessing import Lock, Value
from gaishi.multiprocessing import mp_manager
from gaishi.generators import RandomNumberGenerator
from gaishi.simulators import FeatureVectorSimulator
from gaishi.simulators import GenotypeMatrixSimulator


def simulate_feature_vectors(
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
    intro_prop: float,
    non_intro_prop: float,
    output_prefix: str,
    output_dir: str,
    nfeature: int,
    is_phased: bool,
    is_shuffled: bool,
    force_balanced: bool,
    keep_sim_data: bool,
    seed: int,
    nprocess: int,
) -> None:
    """
    Simulate genomic data and generate feature vectors for training.

    This function orchestrates the entire process of simulating genomic data, labeling the data,
    and generating feature vectors. It supports multiprocessing for efficiency and allows for the
    output to be balanced between introgressed and non-introgressed classes.

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
    intro_prop : float
        Proportion threshold for labeling a window as introgressed.
    non_intro_prop : float
        Proportion threshold for labeling a window as non-introgressed.
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
    force_balanced : bool
        If True, forces the dataset to have an equal number of introgressed and non-introgressed features.
        If nfeature is odd, the extra feature will be added to the intro class.
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
    The function dynamically adjusts the number of features collected from introgressed and
    non-introgressed classes based on the `force_balanced` flag. It ensures the dataset is as balanced
    as possible while respecting the total number of requested features.

    The simulation process is dependent on the correctness of the `LRTrainingDataSimulator`,
    `RandomNumberGenerator`, and the feature computation classes. Ensure these components are correctly
    implemented and tested.

    """
    if nfeature <= 0:
        raise ValueError("nfeature must be positive.")

    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"{output_prefix}.tsv")

    simulator = FeatureVectorSimulator(
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
        intro_prop=intro_prop,
        non_intro_prop=non_intro_prop,
        feature_config_file=feature_config_file,
    )

    total_features = []
    intro_features = []
    non_intro_features = []
    start_rep = 0

    # Determine the number of features to allocate to each class.
    # If nfeature is odd, add the extra feature to the intro_features class.
    num_intro_features = nfeature // 2 + (nfeature % 2)
    num_non_intro_features = nfeature // 2

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

        features0 = [item for item in features if item["Label"] == 0]
        features1 = [item for item in features if item["Label"] == 1]

        if not keep_sim_data:
            for i in range(start_rep, start_rep + nrep):
                shutil.rmtree(os.path.join(output_dir, f"{i}"), ignore_errors=True)

        if force_balanced:
            # Calculate how many items can be added to each list without exceeding the desired count
            to_add_intro = min(num_intro_features - len(intro_features), len(features1))
            to_add_non_intro = min(
                num_non_intro_features - len(non_intro_features), len(features0)
            )

            intro_features.extend(features1[:to_add_intro])
            non_intro_features.extend(features0[:to_add_non_intro])
            is_stopped = (
                len(intro_features) >= num_intro_features
                and len(non_intro_features) >= num_non_intro_features
            )
        else:
            intro_features.extend(features1)
            non_intro_features.extend(features0)
            is_stopped = len(intro_features) + len(non_intro_features) >= nfeature

        print(
            f"Number of obtained features: {len(intro_features)+len(non_intro_features)}"
        )

        start_rep += nrep

    total_features = intro_features + non_intro_features
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


def simulate_genotype_matrices(
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
    num_polymorphisms: int,
    num_genotype_matrices: int,
    output_prefix: str,
    output_dir: str,
    output_h5: bool = True,
    is_phased: bool = True,
    is_sorted: bool = True,
    force_balanced: bool = False,
    keep_sim_data: bool = False,
    seed: int = None,
    nprocess: int = 1,
):
    """
    Simulate genomic data and generate genotype matrices for training.

    This function orchestrates the entire process of simulating genomic data, labeling the data,
    and generating feature vectors. It supports multiprocessing for efficiency and allows for the
    output to be balanced between introgressed and non-introgressed classes.

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
    num_polymorphisms : int
        Number of polymorphisms in a genotype matrix.
    num_genotype_matrices : int
        Total number of simulated genotype matrices in the training data.
    output_prefix : str
        Prefix for output files.
    output_dir : str
        Directory to save output files.
    output_h5 : bool, optional
        If True, create an HDF5 training dataset.
        If False, create a table for the training dataset.
        Default: True.
    is_phased : bool, optional
        Flag indicating if the data is phased. Default: True.
    is_sorted : bool, optional
        Indicates whether the simulated data should be sorted. Default: True.
    force_balanced : bool, optional
        If True, forces the dataset to have an equal number of introgressed and non-introgressed features.
        If nfeature is odd, the extra feature will be added to the intro class. Default: False.
    keep_sim_data : bool, optional
        If False, deletes simulation data directories after processing to save disk space.
        Default: False.
    seed : int, optional
        Seed for random number generation to ensure reproducibility. Default: None.
    nprocess : int, optional
        Number of processes to use for parallel execution. Default: 1.

    Raises
    ------
    ValueError
        If `nfeature` is less than or equal to zero or other invalid parameters are provided.
    SystemExit
        If errors occur in mp_manager.
    """
    simulator = GenotypeMatrixSimulator(
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
        num_polymorphisms=num_polymorphisms,
        num_genotype_matrices=num_genotype_matrices,
        output_prefix=output_prefix,
        output_dir=output_dir,
        output_h5=output_h5,
        is_phased=is_phased,
        is_sorted=is_sorted,
        keep_sim_data=keep_sim_data,
    )

    start_rep = 0
    nrep = nrep
    num_intro = Value("i", 0)
    num_nonintro = Value("i", 0)
    lock = Lock()

    while True:
        generator = RandomNumberGenerator(
            start_rep=start_rep,
            nrep=nrep,
            seed=seed,
        )

        mp_manager(
            job=simulator,
            data_generator=generator,
            nprocess=nprocess,
            force_balanced=force_balanced,
            nintro=num_intro,
            nnonintro=num_nonintro,
            only_intro=False,
            only_non_intro=False,
            lock=lock,
        )

        if num_intro.value + num_nonintro.value >= num_genotype_matrices:
            break
        start_rep += nrep
