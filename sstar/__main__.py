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

import argparse
import os
import sys
import signal
from typing import Optional, Sequence, Type


def _set_sigpipe_handler() -> None:
    """
    Install the default SIGPIPE handler on POSIX systems.
    """
    if os.name == "posix":
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)


def _run_score(args: argparse.Namespace) -> None:
    """
    Run `sstar score` from parsed command-line arguments.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments for the `score` subcommand.
    """
    from sstar.cal_s_star import cal_s_star

    cal_s_star(
        vcf=args.vcf,
        ref_ind_file=args.ref_ind,
        tgt_ind_file=args.tgt_ind,
        anc_allele_file=args.anc_allele,
        win_len=args.win_len,
        win_step=args.win_step,
        output=args.output,
        thread=args.thread,
        match_bonus=args.match_bonus,
        max_mismatch=args.max_mismatch,
        mismatch_penalty=args.mismatch_penalty,
        is_phased=args.phased,
    )


def _run_quantile(args: argparse.Namespace) -> None:
    """
    Run `sstar quantile` from parsed command-line arguments.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments for the `quantile` subcommand.
    """
    from sstar.get_quantile import get_quantile

    get_quantile(
        model=args.model,
        ms_dir=args.ms_dir,
        N0=args.N0,
        nsamp=args.nsamp,
        nreps=args.nreps,
        ref_pop=args.ref_pop,
        ref_size=args.ref_size,
        tgt_pop=args.tgt_pop,
        tgt_size=args.tgt_size,
        mut_rate=args.mut_rate,
        rec_rate=args.rec_rate,
        seq_len=args.seq_len,
        snp_num_range=args.snp_num_range,
        output_dir=args.output_dir,
        thread=args.thread,
        seeds=args.seeds,
        is_phased=args.phased,
        keep_simulated_data=args.keep_simulated_data,
        quantile_step=args.quantile_step,
    )


def _run_threshold(args: argparse.Namespace) -> None:
    """
    Run `sstar threshold` from parsed command-line arguments.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments for the `threshold` subcommand.
    """
    from sstar.cal_threshold import cal_threshold

    cal_threshold(
        simulated_data=args.sim_data,
        score_file=args.score,
        recomb_rate=args.recomb_rate,
        recomb_map=args.recomb_map,
        quantile=args.quantile,
        output=args.output,
        k=args.k,
        phased=args.phased,
    )


def _run_match_pct(args: argparse.Namespace) -> None:
    """
    Run `sstar matchrate` from parsed command-line arguments.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments for the `matchrate` subcommand.
    """
    from sstar.cal_match_rate import cal_match_pct

    cal_match_pct(
        vcf=args.vcf,
        ref_ind_file=args.ref_ind,
        tgt_ind_file=args.tgt_ind,
        src_ind_file=args.src_ind,
        anc_allele_file=args.anc_allele,
        output=args.output,
        thread=args.thread,
        score_file=args.score,
        mapped_region_file=args.mapped_region_file,
        phased=args.phased,
    )


def _run_tract(args: argparse.Namespace) -> None:
    """
    Run `sstar tract` from parsed command-line arguments.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments for the `tract` subcommand.
    """
    from sstar.get_tract import get_tract

    get_tract(
        threshold_file=args.threshold,
        match_pct_files=args.match_pct,
        output_prefix=args.output,
        diff=args.diff,
    )


def _add_common_args(parser: argparse.ArgumentParser) -> None:
    """
    Add command-line arguments shared by multiple subcommands.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        Parser to which the shared arguments are added.
    """
    parser.add_argument(
        "--anc-allele",
        type=str,
        dest="anc_allele",
        default=None,
        help="Path to the BED format file containing ancestral allele information.",
    )
    parser.add_argument(
        "--output", type=str, required=True, help="Path to the output file."
    )
    parser.add_argument(
        "--thread",
        type=int,
        default=1,
        help="Number of threads for multiprocessing. Default: %(default)s.",
    )


def _add_mapped_args(parser: argparse.ArgumentParser) -> None:
    """
    Add the mapped-region command-line argument.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        Parser to which the mapped-region argument is added.
    """
    parser.add_argument(
        "--mapped-region",
        type=str,
        dest="mapped_region_file",
        help="Path to the BED file containing mapped regions.",
    )


def _add_window_args(parser: argparse.ArgumentParser) -> None:
    """
    Add sliding-window command-line arguments.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        Parser to which the window arguments are added.
    """
    parser.add_argument(
        "--win-len",
        type=int,
        default=50000,
        help="Length of sliding windows in base pairs. Default: %(default)s.",
    )
    parser.add_argument(
        "--win-step",
        type=int,
        default=10000,
        help="Step size for sliding windows in base pairs. Default: %(default)s.",
    )


def _add_ref_ind_args(parser: argparse.ArgumentParser) -> None:
    """
    Add the reference-individual command-line argument.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        Parser to which the reference-individual argument is added.
    """
    parser.add_argument(
        "--ref",
        type=str,
        dest="ref_ind",
        required=True,
        help="Path to the file containing reference individual IDs.",
    )


def _add_tgt_ind_args(parser: argparse.ArgumentParser) -> None:
    """
    Add the target-individual command-line argument.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        Parser to which the target-individual argument is added.
    """
    parser.add_argument(
        "--tgt",
        type=str,
        dest="tgt_ind",
        required=True,
        help="Path to the file containing target individual IDs.",
    )


def _add_src_ind_args(parser: argparse.ArgumentParser) -> None:
    """
    Add the source-individual command-line argument.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        Parser to which the source-individual argument is added.
    """
    parser.add_argument(
        "--src",
        type=str,
        dest="src_ind",
        required=True,
        help="Path to the file containing source individual IDs.",
    )


def _add_score_args(parser: argparse.ArgumentParser) -> None:
    """
    Add the S* score-file command-line argument.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        Parser to which the score-file argument is added.
    """
    parser.add_argument(
        "--score",
        type=str,
        dest="score",
        required=True,
        help="Path to the score file generated by `sstar score`.",
    )


def required_length(nmin: int, nmax: int) -> Type[argparse.Action]:
    """
    Create an argparse action that enforces a bounded number of values.

    Parameters
    ----------
    nmin : int
        Minimum allowed number of values.
    nmax : int
        Maximum allowed number of values.

    Returns
    -------
    type[argparse.Action]
        Custom argparse action class that validates the number of values.
    """

    class RequiredLength(argparse.Action):
        def __call__(
            self,
            parser: argparse.ArgumentParser,
            args: argparse.Namespace,
            values: Sequence[str],
            option_string: Optional[str] = None,
        ) -> None:
            """
            Validate the number of parsed values and store them on the namespace.

            Parameters
            ----------
            parser : argparse.ArgumentParser
                Parser that owns the action.
            args : argparse.Namespace
                Namespace where parsed values are stored.
            values : sequence of str
                Values parsed for this option.
            option_string : str, optional
                Option string used on the command line. Default: `None`.
            """
            if not nmin <= len(values) <= nmax:
                raise argparse.ArgumentTypeError(
                    f"argument '{self.dest}' requires between {nmin} and {nmax} arguments"
                )
            setattr(args, self.dest, values)

    return RequiredLength


def _s_star_cli_parser() -> argparse.ArgumentParser:
    """
    Build the top-level S* command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser configured with all S* subcommands.
    """
    top_parser = argparse.ArgumentParser()
    subparsers = top_parser.add_subparsers(dest="subcommand")
    subparsers.required = True

    # score
    parser = subparsers.add_parser("score")
    parser.add_argument(
        "--vcf",
        type=str,
        required=True,
        help="Path to the VCF file containing genotype data.",
    )
    _add_ref_ind_args(parser)
    _add_tgt_ind_args(parser)
    _add_common_args(parser)
    _add_window_args(parser)
    parser.add_argument(
        "--match-bonus",
        type=int,
        default=5000,
        help="Bonus for matching genotypes between two variants. Default: %(default)s.",
    )
    parser.add_argument(
        "--max-mismatch",
        type=int,
        default=5,
        help="Maximum genotype distance allowed before a pair is discarded. Default: %(default)s.",
    )
    parser.add_argument(
        "--mismatch-penalty",
        type=int,
        default=-10000,
        help="Penalty for mismatching genotypes between two variants. Default: %(default)s.",
    )
    parser.add_argument(
        "--phased",
        action="store_true",
        help="Calculate scores on phased haplotypes instead of genotype dosages.",
    )
    parser.set_defaults(runner=_run_score)

    # quantile
    parser = subparsers.add_parser("quantile")
    parser.add_argument(
        "--phased",
        action="store_true",
        help="Run `sstar score` on phased haplotypes during simulation scoring.",
    )
    parser.add_argument(
        "--model",
        type=str,
        required=True,
        help="Path to the demographic model file used for simulation.",
    )
    parser.add_argument(
        "--ms-dir",
        type=str,
        dest="ms_dir",
        required=True,
        help="Path to the directory containing the `ms` executable.",
    )
    parser.add_argument(
        "--N0",
        type=int,
        required=True,
        help="Reference effective population size used for ms parameter scaling.",
    )
    parser.add_argument(
        "--nsamp",
        type=int,
        required=True,
        help="Haploid sample size used in ms simulation.",
    )
    parser.add_argument(
        "--nreps",
        type=int,
        required=True,
        help="Number of simulation replicates.",
    )
    parser.add_argument(
        "--seeds",
        type=int,
        nargs=3,
        default=None,
        help="Three random seed numbers passed to ms with `-seeds`.",
    )
    parser.add_argument(
        "--ref-pop",
        type=str,
        dest="ref_pop",
        required=True,
        help="Name of the reference population in the demographic model.",
    )
    parser.add_argument(
        "--ref-size",
        type=int,
        dest="ref_size",
        required=True,
        help="Haploid sample size of the reference population.",
    )
    parser.add_argument(
        "--tgt-pop",
        type=str,
        dest="tgt_pop",
        required=True,
        help="Name of the target population in the demographic model.",
    )
    parser.add_argument(
        "--tgt-size",
        type=int,
        dest="tgt_size",
        required=True,
        help="Haploid sample size of the target population.",
    )
    parser.add_argument(
        "--mut-rate",
        type=float,
        dest="mut_rate",
        required=True,
        help="Mutation rate per site per generation.",
    )
    parser.add_argument(
        "--rec-rate",
        type=float,
        dest="rec_rate",
        required=True,
        help="Recombination rate per site per generation.",
    )
    parser.add_argument(
        "--seq-len",
        type=int,
        dest="seq_len",
        required=True,
        help="Length of the simulated sequence.",
    )
    parser.add_argument(
        "--snp-num-range",
        type=int,
        nargs=3,
        dest="snp_num_range",
        required=True,
        help="Minimum SNP count, maximum SNP count, and step size used for ms simulations.",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        dest="output_dir",
        required=True,
        help="Path to the output directory.",
    )
    parser.add_argument(
        "--quantile-step",
        type=float,
        dest="quantile_step",
        default=0.005,
        help=(
            "Step size between quantiles from 0.5 to less than 1. "
            "Default: %(default)s."
        ),
    )
    parser.add_argument(
        "--thread",
        type=int,
        default=1,
        help="Number of threads for multiprocessing. Default: %(default)s.",
    )
    parser.add_argument(
        "--keep-simulated-data",
        action="store_true",
        help="Keep intermediate simulation directories instead of deleting them after summarizing results.",
    )
    parser.set_defaults(runner=_run_quantile)

    # threshold
    parser = subparsers.add_parser("threshold")
    parser.add_argument(
        "--phased",
        action="store_true",
        help="Expect phased score rows where sample names carry haplotype suffixes (for example sample_1).",
    )
    _add_score_args(parser)
    parser.add_argument(
        "--sim-data",
        type=str,
        dest="sim_data",
        required=True,
        help="Path to the simulated quantile summary file.",
    )
    parser.add_argument(
        "--recomb-rate",
        type=float,
        dest="recomb_rate",
        default=1e-8,
        help="Uniform recombination rate used when no recombination map is provided. Default: %(default)s.",
    )
    parser.add_argument(
        "--recomb-map",
        type=str,
        dest="recomb_map",
        default=None,
        help="Path to the recombination-map file. If omitted, `--recomb-rate` is used.",
    )
    parser.add_argument(
        "--quantile",
        type=float,
        required=True,
        help="Quantile used to estimate significant S* scores.",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Path to the output threshold file.",
    )
    parser.add_argument(
        "--k",
        type=int,
        default=8,
        help="Basis dimension passed to the MGCV smooth term. Default: %(default)s.",
    )
    parser.set_defaults(runner=_run_threshold)

    # matchrate
    parser = subparsers.add_parser("matchrate")
    parser.add_argument(
        "--phased",
        action="store_true",
        help="Use haplotype-based source matching instead of dosage-based matching.",
    )
    parser.add_argument(
        "--vcf",
        type=str,
        required=True,
        help="Path to the VCF file containing genotype data.",
    )
    _add_ref_ind_args(parser)
    _add_tgt_ind_args(parser)
    _add_src_ind_args(parser)
    _add_common_args(parser)
    _add_mapped_args(parser)
    _add_score_args(parser)
    parser.set_defaults(runner=_run_match_pct)

    # tract
    parser = subparsers.add_parser("tract")
    parser.add_argument(
        "--threshold",
        type=str,
        required=True,
        help="Path to the threshold file generated by `sstar threshold`.",
    )
    parser.add_argument(
        "--match-rate",
        type=str,
        nargs="+",
        dest="match_pct",
        action=required_length(1, 2),
        help="One or two source match-rate files generated by `sstar matchrate`.",
    )
    parser.add_argument(
        "--output-prefix",
        type=str,
        dest="output",
        required=True,
        help="Prefix for output BED files.",
    )
    parser.add_argument(
        "--diff",
        type=float,
        default=0,
        help=(
            "Difference cutoff used to assign tracts when two source match-rate "
            "files are provided. Default: %(default)s."
        ),
    )
    parser.set_defaults(runner=_run_tract)

    return top_parser


def main(arg_list: Optional[Sequence[str]] = None) -> None:
    """
    Run the S* command-line interface.

    Parameters
    ----------
    arg_list : sequence of str, optional
        Command-line arguments to parse instead of `sys.argv`. Default: `None`.
    """
    _set_sigpipe_handler()
    parser = _s_star_cli_parser()
    args = parser.parse_args(arg_list)
    args.runner(args)
