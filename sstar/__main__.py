# Apache License Version 2.0
# Copyright 2022 Xin Huang
#
# Licensed under the Apache License, Version 2.0 (the "License");
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


def _set_sigpipe_handler():
    if os.name == "posix":
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)


def _run_score(args):
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


def _run_quantile(args):
    from sstar.get_quantile import get_quantile

    get_quantile(
        model=args.model,
        ms_dir=args.ms_dir,
        N0=args.N0,
        nsamp=args.nsamp,
        nreps=args.nreps,
        ref_index=args.ref_index,
        ref_size=args.ref_size,
        tgt_index=args.tgt_index,
        tgt_size=args.tgt_size,
        mut_rate=args.mut_rate,
        rec_rate=args.rec_rate,
        seq_len=args.seq_len,
        snp_num_range=args.snp_num_range,
        output_dir=args.output_dir,
        thread=args.thread,
        seeds=args.seeds,
        is_phased=args.phased,
    )


def _run_threshold(args):
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


def _run_match_pct(args):
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
    )


def _run_tract(args):
    from sstar.get_tract import get_tract

    get_tract(
        threshold_file=args.threshold,
        match_pct_files=args.match_pct,
        output_prefix=args.output,
        diff=args.diff,
    )


def _add_common_args(parser):
    parser.add_argument(
        "--anc-allele",
        type=str,
        dest="anc_allele",
        default=None,
        help="name of the BED format file containing ancestral allele information",
    )
    parser.add_argument(
        "--output", type=str, required=True, help="name of the output file"
    )
    parser.add_argument(
        "--thread",
        type=int,
        default=1,
        help="number of threads for multiprocessing; default: 1",
    )


def _add_mapped_args(parser):
    parser.add_argument(
        "--mapped-region",
        type=str,
        dest="mapped_region_file",
        help="name of the BED file containing mapped regions",
    )


def _add_window_args(parser):
    parser.add_argument("--win-len", type=int, default=50000)
    parser.add_argument("--win-step", type=int, default=10000)


def _add_ref_ind_args(parser):
    parser.add_argument("--ref", type=str, dest="ref_ind", required=True)


def _add_tgt_ind_args(parser):
    parser.add_argument("--tgt", type=str, dest="tgt_ind", required=True)


def _add_src_ind_args(parser):
    parser.add_argument("--src", type=str, dest="src_ind", required=True)


def _add_score_args(parser):
    parser.add_argument("--score", type=str, dest="score", required=True)


def required_length(nmin, nmax):
    class RequiredLength(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if not nmin <= len(values) <= nmax:
                raise argparse.ArgumentTypeError(
                    f"argument '{self.dest}' requires between {nmin} and {nmax} arguments"
                )
            setattr(args, self.dest, values)

    return RequiredLength


def _s_star_cli_parser():
    top_parser = argparse.ArgumentParser()
    subparsers = top_parser.add_subparsers(dest="subcommand")
    subparsers.required = True

    # score
    parser = subparsers.add_parser("score")
    parser.add_argument("--vcf", type=str, required=True)
    _add_ref_ind_args(parser)
    _add_tgt_ind_args(parser)
    _add_common_args(parser)
    _add_window_args(parser)
    parser.add_argument("--match-bonus", type=int, default=5000)
    parser.add_argument("--max-mismatch", type=int, default=5)
    parser.add_argument("--mismatch-penalty", type=int, default=-10000)
    parser.add_argument("--phased", action="store_true")
    parser.set_defaults(runner=_run_score)

    # quantile  ✅ FIXED
    parser = subparsers.add_parser("quantile")
    parser.add_argument("--phased", action="store_true")
    parser.add_argument("--model", type=str, required=True)
    parser.add_argument("--ms-dir", type=str, dest="ms_dir", required=True)
    parser.add_argument("--N0", type=int, required=True)
    parser.add_argument("--nsamp", type=int, required=True)
    parser.add_argument("--nreps", type=int, required=True)
    parser.add_argument("--seeds", type=int, nargs=3, default=None)
    parser.add_argument("--ref-index", type=int, dest="ref_index", required=True)
    parser.add_argument("--ref-size", type=int, dest="ref_size", required=True)
    parser.add_argument("--tgt-index", type=int, dest="tgt_index", required=True)
    parser.add_argument("--tgt-size", type=int, dest="tgt_size", required=True)
    parser.add_argument("--mut-rate", type=float, dest="mut_rate", required=True)
    parser.add_argument("--rec-rate", type=float, dest="rec_rate", required=True)
    parser.add_argument("--seq-len", type=int, dest="seq_len", required=True)
    parser.add_argument(
        "--snp-num-range", type=int, nargs=3, dest="snp_num_range", required=True
    )
    parser.add_argument("--output-dir", type=str, dest="output_dir", required=True)
    parser.add_argument("--thread", type=int, default=1)
    parser.set_defaults(runner=_run_quantile)

    # threshold
    parser = subparsers.add_parser("threshold")
    parser.add_argument("--phased", action="store_true")
    _add_score_args(parser)
    parser.add_argument("--sim-data", type=str, dest="sim_data", required=True)
    parser.add_argument("--recomb-rate", type=float, dest="recomb_rate", default=1e-8)
    parser.add_argument("--recomb-map", type=str, dest="recomb_map", default=None)
    parser.add_argument("--quantile", type=float, required=True)
    parser.add_argument("--output", type=str, required=True)
    parser.add_argument("--k", type=int, default=8)
    parser.set_defaults(runner=_run_threshold)

    # matchrate
    parser = subparsers.add_parser("matchrate")
    parser.add_argument("--vcf", type=str, required=True)
    _add_ref_ind_args(parser)
    _add_tgt_ind_args(parser)
    _add_src_ind_args(parser)
    _add_common_args(parser)
    _add_mapped_args(parser)
    _add_score_args(parser)
    parser.set_defaults(runner=_run_match_pct)

    # tract
    parser = subparsers.add_parser("tract")
    parser.add_argument("--threshold", type=str, required=True)
    parser.add_argument(
        "--match-rate",
        type=str,
        nargs="+",
        dest="match_pct",
        action=required_length(1, 2),
    )
    parser.add_argument("--output-prefix", type=str, dest="output", required=True)
    parser.add_argument("--diff", type=float, default=0)
    parser.set_defaults(runner=_run_tract)

    return top_parser


def main(arg_list=None):
    _set_sigpipe_handler()
    parser = _s_star_cli_parser()
    args = parser.parse_args(arg_list)
    args.runner(args)
