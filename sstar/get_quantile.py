# Copyright 2025 Xin Huang
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


import os
import demes
import glob
import subprocess
import numpy as np
import pandas as pd
import multiprocessing
from multiprocessing import Process, Queue
from scipy.stats import norm
from scipy.stats import nbinom
from typing import List, Optional, Union


def get_quantile(
    model: str,
    ms_dir: str,
    N0: int,
    nsamp: int,
    nreps: int,
    ref_index: int,
    ref_size: int,
    tgt_index: int,
    tgt_size: int,
    mut_rate: float,
    rec_rate: float,
    seq_len: int,
    snp_num_range: List[int],
    output_dir: str,
    thread: int,
    seeds: Optional[List[int]],
) -> None:
    """
    Calculates quantiles of expected S* from ms simulations.

    Parameters
    ----------
    model : str
        Path to the file containing the demographic model for simulation.
    ms_dir : str
        Directory containing the `ms` simulation program.
    N0 : int
        Effective population size used in `ms` simulation.
    nsamp : int
        Total sample size (haploid) used in `ms` simulation.
    nreps : int
        Number of simulation replicates.
    ref_index : int
        Index (1-based) of the reference population in the demographic model.
    ref_size : int
        Sample size (haploid) of the reference population.
    tgt_index : int
        Index (1-based) of the target population in the demographic model.
    tgt_size : int
        Sample size (haploid) of the target population.
    mut_rate : float
        Mutation rate per base per generation.
    rec_rate : float
        Recombination rate per base per generation.
    seq_len : int
        Length of the simulated sequence.
    snp_num_range : list of int
        Range of SNP numbers used in `ms` simulation. Should be a list of three integers:
        [min_snp, max_snp, step].
    output_dir : str
        Path to the output directory.
    thread : int
        Number of threads to use.
    seeds : list of int, optional
        Three random seed integers used for simulation reproducibility.
    """
    if seeds is not None:
        np.random.seed(np.sum(seeds))
    output_dir = os.path.abspath(output_dir)
    if not os.path.exists(output_dir):
        subprocess.call(["mkdir", output_dir])
    _generate_mut_rec_combination(N0, nreps, mut_rate, rec_rate, seq_len, output_dir)
    _run_ms_simulation(
        model,
        ms_dir,
        N0,
        nsamp,
        nreps,
        ref_index,
        ref_size,
        tgt_index,
        tgt_size,
        seq_len,
        snp_num_range,
        output_dir,
        thread,
        seeds,
    )
    _summary(output_dir, rec_rate)


def _generate_mut_rec_combination(
    N0: int, nreps: int, mut_rate: float, rec_rate: float, seq_len: int, output_dir: str
) -> None:
    """
    Generate combinations of mutation and recombination rates for ms simulations.

    Parameters
    ----------
    N0 : int
        Effective population size used in ms simulation.
    nreps : int
        Number of replicate simulations.
    mut_rate : float
        Baseline mutation rate per base per generation.
    rec_rate : float
        Baseline recombination rate per base per generation.
    seq_len : int
        Length of the simulated sequence.
    output_dir : str
        Path to the directory where the output file will be written.
    """
    scaled_mut_rate = 4*N0*mut_rate*seq_len
    scaled_rec_rate = 4*N0*rec_rate*seq_len
    mut_rate_list = norm.rvs(loc=scaled_mut_rate, scale=0.233, size=nreps)
    rec_rate_list = nbinom.rvs(n=0.5, p=0.5/(0.5+scaled_rec_rate), size=nreps)

    rates = f'{output_dir}/rates.combination'
    with open(rates, 'w') as o:
        for i in range(len(mut_rate_list)):
            if mut_rate_list[i] < 0.001: mut_rate_list[i] = 0.001
            if rec_rate_list[i] < 0.001: rec_rate_list[i] = 0.001
            mut_rate = mut_rate_list[i]
            rec_rate = rec_rate_list[i]
            o.write(f'{mut_rate}\t{rec_rate}\n')


def _run_ms_simulation(
    model: str,
    ms_dir: str,
    N0: int,
    nsamp: int,
    nreps: int,
    ref_index: int,
    ref_size: int,
    tgt_index: int,
    tgt_size: int,
    seq_len: int,
    snp_num_range: List[int],
    output_dir: str,
    thread: int,
    seeds: List[int],
) -> None:
    """
    Helper function for running ms simulation.

    Parameters
    ----------
    model : str
        Name of file containing the demographic model for simulation.
    ms_dir : str
        Name of the directory containing the ms program.
    N0 : int
        N0 used in ms simulation.
    nsamp : int
        Sample size (haploid) used in ms simulation.
    nreps : int
        Number of replicates used in ms simulation.
    ref_index : int
        Index of the reference population in the demographic model (start from 1).
    ref_size : int
        Sample size (haploid) of the reference population.
    tgt_index : int
        Index of the target population in the demographic model (start from 1).
    tgt_size : int
        Sample size (haploid) of the target population.
    seq_len : int
        Length of simulated sequence.
    snp_num_range : list of int
        Range of SNP numbers in ms simulation; the first parameter is the minimum SNP number,
        the second parameter is the maximum SNP number, the third parameter is the step size.
    output_dir : str
        Name of the output directory.
    thread : int
        Number of threads.
    seeds : list of int
        Three random seed numbers used in ms simulation.
    """
    graph = demes.load(model)
    samples = np.zeros(len(graph.demes))
    samples[ref_index - 1] = ref_size
    samples[tgt_index - 1] = tgt_size
    ms_params = demes.to_ms(graph, N0=N0, samples=samples)
    snp_num_list = np.arange(
        snp_num_range[0], snp_num_range[1] + snp_num_range[2], snp_num_range[2]
    )

    ms_exec = os.path.abspath(ms_dir) + "/ms"
    rates = f"{output_dir}/rates.combination"
    ref_list = f"{output_dir}/sim.ref.list"
    tgt_list = f"{output_dir}/sim.tgt.list"
    ref_size = int(ref_size / 2)
    tgt_size = int(tgt_size / 2)

    if ref_index == tgt_index:
        raise Exception(
            "The reference population should be different from the target population."
        )
    elif ref_index < tgt_index:
        with open(ref_list, "w") as o:
            for i in range(ref_size):
                o.write(f"ms_{i}\n")
        with open(tgt_list, "w") as o:
            for i in range(tgt_size):
                o.write(f"ms_{i+ref_size}\n")
    else:
        with open(ref_list, "w") as o:
            for i in range(ref_size):
                o.write(f"ms_{i+tgt_size}\n")
        with open(tgt_list, "w") as o:
            for i in range(tgt_size):
                o.write(f"ms_{i}\n")

    try:
        from pytest_cov.embed import cleanup_on_sigterm
    except ImportError:
        pass
    else:
        cleanup_on_sigterm()

    in_queue, out_queue = Queue(), Queue()
    workers = [
        Process(
            target=_run_ms_simulation_worker,
            args=(
                in_queue,
                out_queue,
                output_dir,
                rates,
                ms_exec,
                nsamp,
                nreps,
                seq_len,
                ms_params,
                ref_list,
                tgt_list,
                seeds,
            ),
        )
        for ii in range(thread)
    ]

    for snp_num in snp_num_list:
        in_queue.put(snp_num)

    try:
        for worker in workers:
            worker.start()
        for snp_num in snp_num_list:
            out_queue.get()
        for worker in workers:
            worker.terminate()
    finally:
        for worker in workers:
            worker.join()


def _run_ms_simulation_worker(
    in_queue: multiprocessing.Queue,
    out_queue: multiprocessing.Queue,
    output_dir: str,
    rates: str,
    ms_exec: str,
    nsamp: int,
    nreps: int,
    seq_len: int,
    ms_params: List[str],
    ref_list: str,
    tgt_list: str,
    seeds: List[int],
) -> None:
    """
    Worker function for running ms simulation.

    Parameters
    ----------
    in_queue : multiprocessing.Queue
        Queue instance to receive parameters from the manager.
    out_queue : multiprocessing.Queue
        Queue instance to send results back to the manager.
    output_dir : str
        Name of the output directory.
    rates : str
        Path to the file containing different combinations of mutation and recombination rates.
    ms_exec : str
        Path to the ms program executable.
    nsamp : int
        Sample size (haploid) used in ms simulation.
    nreps : int
        Number of replicates used in ms simulation.
    seq_len : int
        Length of the simulated sequence.
    ms_params : list of str
        List of parameters to be passed to ms.
    ref_list : str
        Path to the file containing individuals from the reference population.
    tgt_list : str
        Path to the file containing individuals from the target population.
    seeds : list of int
        Three random seed numbers used in ms simulation.
    """
    while True:
        snp_num = in_queue.get()
        output_subdir = f"{output_dir}/{snp_num}"
        output_ms = f"{output_subdir}/sim.ms"
        output_vcf = f"{output_subdir}/sim.vcf"
        output_score = f"{output_subdir}/sim.score"
        output_quantile = f"{output_subdir}/sim.quantile"
        ms_script = f"{output_subdir}/run_ms.sh"

        if seeds is not None:
            cmd = " ".join(
                [
                    "cat",
                    rates,
                    "|",
                    ms_exec,
                    str(nsamp),
                    str(nreps),
                    "-seeds",
                    " ".join([str(s) for s in seeds]),
                    "-t",
                    "tbs",
                    "-r",
                    "tbs",
                    str(seq_len),
                    "-s",
                    str(snp_num),
                    ms_params,
                    ">",
                    output_ms,
                ]
            )
        else:
            cmd = " ".join(
                [
                    "cat",
                    rates,
                    "|",
                    ms_exec,
                    str(nsamp),
                    str(nreps),
                    "-t",
                    "tbs",
                    "-r",
                    "tbs",
                    str(seq_len),
                    "-s",
                    str(snp_num),
                    ms_params,
                    ">",
                    output_ms,
                ]
            )

        if os.path.exists(output_subdir) is False:
            subprocess.call(["mkdir", output_subdir])
        with open(ms_script, "w") as o:
            o.write(cmd + "\n")
        subprocess.call(["bash", ms_script])
        _ms2vcf(output_ms, output_vcf, nsamp, seq_len)
        subprocess.call(
            [
                "sstar",
                "score",
                "--vcf",
                output_vcf,
                "--ref",
                ref_list,
                "--tgt",
                tgt_list,
                "--output",
                output_score,
                "--win-len",
                str(seq_len),
                "--win-step",
                str(seq_len),
                "--thread",
                "1",
            ]
        )
        _cal_quantile(output_score, output_quantile, snp_num)
        out_queue.put("Finished")


def _ms2vcf(
    ms_file: str,
    vcf_file: str,
    nsamp: int,
    seq_len: int,
    ploidy: int = 2
) -> None:
    """
    Converts ms output files into the VCF format.

    Parameters
    ----------
    ms_file : str
        Name of the ms file (input).
    vcf_file : str
        Name of the VCF file (output).
    nsamp : int
        Number of haploid genomes.
    seq_len : int
        Sequence length.
    ploidy : int, optional
        Ploidy of each individual. Default is 2.
    """
    data = []
    i = -1
    header = "##fileformat=VCFv4.2\n"
    header += '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(
        ["ms_" + str(i) for i in range(int(nsamp / ploidy))]
    )

    with open(ms_file, "r") as f:
        f.readline()
        f.readline()
        for l in f.readlines():
            if l.startswith("//"):
                i += 1
                data.append({})
                data[i]["pos"] = []
                data[i]["geno"] = []
            elif l.startswith("positions"):
                data[i]["pos"] = l.rstrip().split(" ")[1:]
            elif l.startswith("0") or l.startswith("1"):
                data[i]["geno"].append(l.rstrip())

    shift = 0
    with open(vcf_file, "w") as o:
        o.write(header + "\n")
        for i in range(len(data)):
            for j in range(len(data[i]["pos"])):
                pos = int(seq_len * float(data[i]["pos"][j])) + shift
                genotypes = "".join(
                    [data[i]["geno"][k][j] for k in range(len(data[i]["geno"]))]
                )
                genotypes = "\t".join(
                    [
                        a + "|" + b
                        for a, b in zip(genotypes[0::ploidy], genotypes[1::ploidy])
                    ]
                )
                o.write(f"1\t{pos}\t.\tA\tT\t100\tPASS\t.\tGT\t{genotypes}\n")
            shift += seq_len


def _cal_quantile(
    in_file: str,
    out_file: str,
    snp_num: int
) -> None:
    """
    Helper function for calculating quantiles of expected S* with a given SNP number.

    Parameters
    ----------
    in_file : str
        Name of the input file containing S* scores.
    out_file : str
        Name of the output file.
    snp_num : int
        Number of SNPs used in the simulation.
    """
    df = pd.read_csv(in_file, sep="\t").dropna()
    quantiles = np.arange(0.5, 1, 0.005)
    mean_df = (
        df.groupby(["chrom", "start", "end"], as_index=False)["S*_score"]
        .mean()
        .dropna()
    )
    scores = np.quantile(mean_df["S*_score"], quantiles)
    with open(out_file, "w") as o:
        o.write("S*_score\tSNP_num\tquantile\n")
        for i in range(len(scores)):
            o.write(f"{scores[i]}\t{snp_num}\t{quantiles[i]}\n")


def _summary(
    output_dir: str,
    rec_rate: float
) -> None:
    """
    Helper function to summarize quantiles of expected S* from different SNP numbers.

    Parameters
    ----------
    output_dir : str
        Name of the output directory.
    rec_rate : float
        Recombination rate.
    """
    all_files = glob.glob(f"{output_dir}/*/*.quantile")
    li = []
    for filename in all_files:
        df = pd.read_csv(filename, sep="\t")
        li.append(df)

    df = pd.concat(li, ignore_index=True).round(3)
    df["log(local_recomb_rate)"] = np.log10(rec_rate)
    df.sort_values(by=["SNP_num", "quantile"]).to_csv(
        f"{output_dir}/quantile.summary.txt", sep="\t", index=False
    )
