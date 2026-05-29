# Apache License Version 2.0
# Copyright 2022 Xin Huang
#
# Licensed under the Apache License Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import os
import demes
import glob
import shutil
import subprocess
import numpy as np
import pandas as pd
from multiprocessing import Process, Queue
from scipy.stats import norm
from scipy.stats import nbinom
from typing import Optional, Sequence


# ------------------------------------------------------------------
# SAFE SUBPROCESS WRAPPER FOR "sstar score"
# ------------------------------------------------------------------
def _safe_sstar_score(cmd: list) -> None:
    """
    Run an allow-listed `sstar score` subprocess command.

    Parameters
    ----------
    cmd : list
        Command tokens for an `sstar score` invocation.
    """
    if not isinstance(cmd, list) or len(cmd) < 2:
        raise ValueError(f"Unexpected subprocess invocation: {cmd}")

    if cmd[0] != "sstar" or cmd[1] != "score":
        raise ValueError(f"Unexpected subprocess invocation: {cmd}")

    allowed_flags = {
        "--vcf",
        "--ref",
        "--tgt",
        "--output",
        "--win-len",
        "--win-step",
        "--thread",
        "--phased",
    }

    i = 2
    while i < len(cmd):
        tok = cmd[i]
        if tok not in allowed_flags:
            raise ValueError(f"Unexpected flag in sstar invocation: {tok}")

        if tok == "--phased":
            i += 1
            continue

        if i + 1 >= len(cmd):
            raise ValueError(f"Missing value for {tok}")

        i += 2

    subprocess.run(cmd, check=True)


def get_quantile(
    model: str,
    ms_dir: str,
    N0: int,
    nsamp: int,
    nreps: int,
    ref_pop: str,
    ref_size: int,
    tgt_pop: str,
    tgt_size: int,
    mut_rate: float,
    rec_rate: float,
    seq_len: int,
    snp_num_range: Sequence[int],
    output_dir: str,
    thread: int,
    seeds: Optional[list],
    quantile_step: float,
    is_phased: bool = False,
    keep_simulated_data: bool = False,
) -> None:
    """
    Calculate quantiles of expected S* scores from ms simulations.

    Parameters
    ----------
    model : str
        Path to the demographic model file used for simulation.
    ms_dir : str
        Path to the directory containing the `ms` executable.
    N0 : int
        Reference effective population size used for ms parameter scaling.
    nsamp : int
        Haploid sample size used in ms simulation.
    nreps : int
        Number of simulation replicates.
    ref_pop : str
        Name of the reference population in the demographic model.
    ref_size : int
        Haploid sample size of the reference population.
    tgt_pop : str
        Name of the target population in the demographic model.
    tgt_size : int
        Haploid sample size of the target population.
    mut_rate : float
        Mutation rate per site per generation.
    rec_rate : float
        Recombination rate per site per generation.
    seq_len : int
        Length of the simulated sequence.
    snp_num_range : sequence of int
        Three values specifying minimum SNP count, maximum SNP count, and step size.
    output_dir : str
        Path to the output directory.
    thread : int
        Number of worker processes.
    seeds : list or None
        Three random seed numbers used in ms simulation. If None, ms uses its default seeding behavior.
    quantile_step : float
        Step size between quantiles from 0.5 to less than 1.
    is_phased : bool, optional
        If True, run `sstar score` with `--phased`. Default: `False`.
    keep_simulated_data : bool, optional
        If True, keep intermediate simulation directories. Default: `False`.
    """
    _validate_quantile_step(quantile_step)
    if seeds is not None:
        np.random.seed(np.sum(seeds))
    output_dir = os.path.abspath(output_dir)
    if os.path.exists(output_dir) is False:
        os.makedirs(output_dir, exist_ok=True)  # no subprocess mkdir
    _generate_mut_rec_combination(N0, nreps, mut_rate, rec_rate, seq_len, output_dir)
    simulated_snp_dirs = _run_ms_simulation(
        model,
        ms_dir,
        N0,
        nsamp,
        nreps,
        ref_pop,
        ref_size,
        tgt_pop,
        tgt_size,
        seq_len,
        snp_num_range,
        output_dir,
        thread,
        seeds,
        quantile_step,
        is_phased,
    )
    _summary(output_dir, rec_rate)
    if not keep_simulated_data:
        _cleanup_simulated_data(output_dir, simulated_snp_dirs)


def _validate_quantile_step(quantile_step: float) -> None:
    """
    Validate the quantile step size.

    Parameters
    ----------
    quantile_step : float
        Step size between quantiles from 0.5 to less than 1.
    """
    if not np.isfinite(quantile_step) or quantile_step <= 0 or quantile_step >= 0.5:
        raise ValueError("quantile_step must be greater than 0 and less than 0.5")


def _generate_mut_rec_combination(
    N0: int, nreps: int, mut_rate: float, rec_rate: float, seq_len: int, output_dir: str
) -> None:
    """
    Generate scaled mutation-rate and recombination-rate combinations for ms.

    Parameters
    ----------
    N0 : int
        Reference effective population size used for ms parameter scaling.
    nreps : int
        Number of simulation replicates.
    mut_rate : float
        Mutation rate per site per generation.
    rec_rate : float
        Recombination rate per site per generation.
    seq_len : int
        Length of the simulated sequence.
    output_dir : str
        Path to the output directory.
    """
    scaled_mut_rate = 4 * N0 * mut_rate * seq_len
    scaled_rec_rate = 4 * N0 * rec_rate * seq_len
    mut_rate_list = norm.rvs(loc=scaled_mut_rate, scale=0.233, size=nreps)
    rec_rate_list = nbinom.rvs(n=0.5, p=0.5 / (0.5 + scaled_rec_rate), size=nreps)

    rates = f"{output_dir}/rates.combination"
    with open(rates, "w") as o:
        for i in range(len(mut_rate_list)):
            if mut_rate_list[i] < 0.001:
                mut_rate_list[i] = 0.001
            if rec_rate_list[i] < 0.001:
                rec_rate_list[i] = 0.001
            m = mut_rate_list[i]
            r = rec_rate_list[i]
            o.write(f"{m}\t{r}\n")


def _get_pop_index(graph: demes.Graph, pop_name: str) -> int:
    """
    Get the 1-based index of a deme in a demes graph.

    Parameters
    ----------
    graph : demes.Graph
        Demographic graph loaded by demes.
    pop_name : str
        Deme name to resolve.

    Returns
    -------
    int
        1-based population index for use in ms population ordering.
    """
    names = [deme.name for deme in graph.demes]
    if pop_name not in names:
        raise ValueError(
            f"Population '{pop_name}' not found in model. "
            f"Available populations: {', '.join(names)}"
        )
    return names.index(pop_name) + 1


def _run_ms_simulation(
    model: str,
    ms_dir: str,
    N0: int,
    nsamp: int,
    nreps: int,
    ref_pop: str,
    ref_size: int,
    tgt_pop: str,
    tgt_size: int,
    seq_len: int,
    snp_num_range: Sequence[int],
    output_dir: str,
    thread: int,
    seeds: Optional[list],
    quantile_step: float,
    is_phased: bool = False,
) -> list:
    """
    Run ms simulations across the requested SNP-count range.

    Parameters
    ----------
    model : str
        Path to the demographic model file used for simulation.
    ms_dir : str
        Path to the directory containing the `ms` executable.
    N0 : int
        Reference effective population size used for ms parameter scaling.
    nsamp : int
        Haploid sample size used in ms simulation.
    nreps : int
        Number of simulation replicates.
    ref_pop : str
        Name of the reference population in the demographic model.
    ref_size : int
        Haploid sample size of the reference population.
    tgt_pop : str
        Name of the target population in the demographic model.
    tgt_size : int
        Haploid sample size of the target population.
    seq_len : int
        Length of the simulated sequence.
    snp_num_range : sequence of int
        Three values specifying minimum SNP count, maximum SNP count, and step size.
    output_dir : str
        Path to the output directory.
    thread : int
        Number of worker processes.
    seeds : list or None
        Three random seed numbers used in ms simulation. If None, ms uses its default seeding behavior.
    quantile_step : float
        Step size between quantiles from 0.5 to less than 1.
    is_phased : bool, optional
        If True, run `sstar score` with `--phased`. Default: `False`.

    Returns
    -------
    list
        SNP-count directory names created by this run.
    """
    graph = demes.load(model)
    ref_index = _get_pop_index(graph, ref_pop)
    tgt_index = _get_pop_index(graph, tgt_pop)
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
                quantile_step,
                is_phased,
            ),
        )
        for _ii in range(thread)
    ]

    for snp_num in snp_num_list:
        in_queue.put(snp_num)

    try:
        for worker in workers:
            worker.start()
        for _snp_num in snp_num_list:
            out_queue.get()
        for worker in workers:
            worker.terminate()
    finally:
        for worker in workers:
            worker.join()

    return [str(snp_num) for snp_num in snp_num_list]


def _run_ms_simulation_worker(
    in_queue: Queue,
    out_queue: Queue,
    output_dir: str,
    rates: str,
    ms_exec: str,
    nsamp: int,
    nreps: int,
    seq_len: int,
    ms_params: str,
    ref_list: str,
    tgt_list: str,
    seeds: Optional[list],
    quantile_step: float,
    is_phased: bool = False,
) -> None:
    """
    Worker process for running ms simulations and score quantile calculation.

    Parameters
    ----------
    in_queue : multiprocessing.Queue
        Queue providing SNP counts to simulate.
    out_queue : multiprocessing.Queue
        Queue receiving completion messages.
    output_dir : str
        Path to the output directory.
    rates : str
        Path to the file containing scaled mutation and recombination rates.
    ms_exec : str
        Path to the `ms` executable.
    nsamp : int
        Haploid sample size used in ms simulation.
    nreps : int
        Number of simulation replicates.
    seq_len : int
        Length of the simulated sequence.
    ms_params : str
        Command-line options generated from the demographic model for ms.
    ref_list : str
        Path to the file containing simulated reference individual IDs.
    tgt_list : str
        Path to the file containing simulated target individual IDs.
    seeds : list or None
        Three random seed numbers used in ms simulation. If None, ms uses its default seeding behavior.
    quantile_step : float
        Step size between quantiles from 0.5 to less than 1.
    is_phased : bool, optional
        If True, run `sstar score` with `--phased`. Default: `False`.
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

        os.makedirs(output_subdir, exist_ok=True)  # no subprocess mkdir
        with open(ms_script, "w") as o:
            o.write(cmd + "\n")
        subprocess.call(
            ["bash", ms_script]
        )  # sourcery skip: python.lang.security.audit.dangerous-subprocess-use-audit
        _ms2vcf(output_ms, output_vcf, nsamp, seq_len)

        score_cmd = [
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
        if is_phased:
            score_cmd.append("--phased")

        _safe_sstar_score(score_cmd)  # safe wrapper (no linter error)

        _cal_quantile(output_score, output_quantile, snp_num, quantile_step)
        out_queue.put("Finished")


def _ms2vcf(
    ms_file: str, vcf_file: str, nsamp: int, seq_len: int, ploidy: int = 2
) -> None:
    """
    Convert ms output into VCF format.

    Parameters
    ----------
    ms_file : str
        Path to the ms output file.
    vcf_file : str
        Path to the output VCF file.
    nsamp : int
        Haploid sample size used in ms simulation.
    seq_len : int
        Length of the simulated sequence.
    ploidy : int, optional
        Ploidy used to combine haplotypes into genotype strings. Default: `2`.
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
    in_file: str, out_file: str, snp_num: int, quantile_step: float
) -> None:
    """
    Calculate quantiles of expected S* for one simulated SNP count.

    Parameters
    ----------
    in_file : str
        Path to the simulated S* score file.
    out_file : str
        Path to the output quantile file.
    snp_num : int
        SNP count used in the ms simulation.
    quantile_step : float
        Step size between quantiles from 0.5 to less than 1.
    """
    _validate_quantile_step(quantile_step)
    df = pd.read_csv(in_file, sep="\t")
    df["S*_score"] = pd.to_numeric(df["S*_score"], errors="coerce")
    df = df.dropna(subset=["S*_score"])

    if df.empty:
        raise RuntimeError(f"No valid S* scores found in {in_file}")

    quantiles = np.arange(0.5, 1, quantile_step)
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


def _summary(output_dir: str, rec_rate: float) -> None:
    """
    Summarize expected S* quantiles across simulated SNP counts.

    Parameters
    ----------
    output_dir : str
        Path to the output directory containing per-SNP-count quantile files.
    rec_rate : float
        Recombination rate written to the summary as `log(local_recomb_rate)`.
    """
    all_files = glob.glob(f"{output_dir}/*/*.quantile")
    li = []
    for filename in all_files:
        df = pd.read_csv(filename, sep="\t")
        li.append(df)

    df = pd.concat(li, ignore_index=True)
    df["log(local_recomb_rate)"] = np.log10(rec_rate)
    df.sort_values(by=["SNP_num", "quantile"]).to_csv(
        f"{output_dir}/quantile.summary.txt", sep="\t", index=False
    )


def _cleanup_simulated_data(output_dir: str, simulated_snp_dirs: list) -> None:
    """
    Remove intermediate files generated during ms simulations.

    Parameters
    ----------
    output_dir : str
        Path to the output directory.
    simulated_snp_dirs : list
        SNP-count directories created by this run.
    """
    generated_files = {"rates.combination", "sim.ref.list", "sim.tgt.list"}
    simulated_dir_required_files = {
        "sim.ms",
        "sim.vcf",
        "sim.score",
        "sim.quantile",
        "run_ms.sh",
    }

    for filename in generated_files:
        path = os.path.join(output_dir, filename)
        if os.path.isfile(path):
            os.remove(path)

    for dirname in simulated_snp_dirs:
        subdir = os.path.join(output_dir, dirname)
        if not os.path.isdir(subdir):
            continue
        existing_files = set(os.listdir(subdir))
        if simulated_dir_required_files.issubset(existing_files):
            shutil.rmtree(subdir)
