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

import argparse, os, sys, signal


def _set_sigpipe_handler():
    if os.name == "posix":
        # Set signal handler for SIGPIPE to quietly kill the program.
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)


def _run_score(args):
    from sstar.cal_s_star import cal_s_star
    cal_s_star(vcf=args.vcf, ref_ind_file=args.ref_ind, tgt_ind_file=args.tgt_ind, anc_allele_file=args.anc_allele, win_len=args.win_len, win_step=args.win_step, output=args.output, thread=args.thread, match_bonus=args.match_bonus, max_mismatch=args.max_mismatch, mismatch_penalty=args.mismatch_penalty)


def _run_quantile(args):
    from sstar.get_quantile import get_quantile
    get_quantile(model=args.model, ms_dir=args.ms_dir, N0=args.N0, nsamp=args.nsamp, nreps=args.nreps, ref_index=args.ref_index, ref_size=args.ref_size, tgt_index=args.tgt_index, tgt_size=args.tgt_size, mut_rate=args.mut_rate, rec_rate=args.rec_rate, seq_len=args.seq_len, snp_num_range=args.snp_num_range, output_dir=args.output_dir, thread=args.thread, seeds=args.seeds)


def _run_threshold(args):
    from sstar.cal_threshold import cal_threshold
    cal_threshold(simulated_data=args.sim_data, score_file=args.score, recomb_rate=args.recomb_rate, recomb_map=args.recomb_map, quantile=args.quantile, output=args.output, k=args.k)


def _run_match_pct(args):
    from sstar.cal_match_rate import cal_match_pct
    cal_match_pct(vcf=args.vcf, ref_ind_file=args.ref_ind, tgt_ind_file=args.tgt_ind, src_ind_file=args.src_ind, anc_allele_file=args.anc_allele, output=args.output, thread=args.thread, score_file=args.score, mapped_region_file=args.mapped_region_file)


def _run_tract(args):
    from sstar.get_tract import get_tract
    get_tract(threshold_file=args.threshold, match_pct_files=args.match_pct, output_prefix=args.output, diff=args.diff)


def _add_common_args(parser):
    parser.add_argument('--anc-allele', type=str, dest='anc_allele', default=None, help='name of the BED format file containing ancestral allele information, otherwise assuming the REF allele is the ancestral allele and the ALT allele is the derived allele; default: None')
    parser.add_argument('--output', type=str, required=True, help='name of the output file')
    parser.add_argument('--thread', type=int, default=1, help='number of threads for multiprocessing; default: 1')


def _add_mapped_args(parser):
    parser.add_argument('--mapped-region', type=str, dest='mapped_region_file', help='name of the BED file containing mapped regions')


def _add_window_args(parser):
    parser.add_argument('--win-len', type=int, dest='win_len', default=50000, help='length of the window to calculate S* scores; default: 50000')
    parser.add_argument('--win-step', type=int, dest='win_step', default=10000, help='step size for moving windows along genomes; default: 10000')


def _add_ref_ind_args(parser):
    parser.add_argument('--ref', type=str, dest='ref_ind', required=True, help='name of the file containing information for samples without introgression')


def _add_tgt_ind_args(parser):
    parser.add_argument('--tgt', type=str, dest='tgt_ind', required=True, help='name of the file containing information for samples for detecting introgression')


def _add_src_ind_args(parser):
    parser.add_argument('--src', type=str, dest='src_ind', required=True, help='name of the file containing information for samples from source populations')


def _add_score_args(parser):
    parser.add_argument('--score', type=str, dest='score', required=True, default=None, help='name of the file containing S* scores calculated by `sstar score`')


def required_length(nmin,nmax):
    class RequiredLength(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if not nmin<=len(values)<=nmax:
                msg='argument "{f}" requires between {nmin} and {nmax} arguments'.format(
                    f=self.dest,nmin=nmin,nmax=nmax)
                raise argparse.ArgumentTypeError(msg)
            setattr(args, self.dest, values)
    return RequiredLength


def _s_star_cli_parser():
    top_parser = argparse.ArgumentParser()
    subparsers = top_parser.add_subparsers(dest="subcommand")
    subparsers.required = True

    # Arguments for score subcommand
    parser = subparsers.add_parser('score', help='calculate S* scores from VCF files')
    parser.add_argument('--vcf', type=str, dest='vcf', required=True, help='name of the VCF file containing genotypes from samples')
    _add_ref_ind_args(parser)
    _add_tgt_ind_args(parser)
    _add_common_args(parser)
    _add_window_args(parser)
    parser.add_argument('--match-bonus', type=int, dest='match_bonus', default=5000, help='bonus for matching genotypes of two different variants; default: 5000')
    parser.add_argument('--max-mismatch', type=int, dest='max_mismatch', default=5, help='maximum genotype distance allowed; default: 5')
    parser.add_argument('--mismatch-penalty', type=int, dest='mismatch_penalty', default=-10000, help='penalty for mismatching genotypes of two different variants; default: -10000')
    parser.set_defaults(runner=_run_score)

    # Arguments for quantile subcommand
    parser = subparsers.add_parser('quantile', help='calculate quantiles for S* scores from simulated data without introgression')
    parser.add_argument('--model', type=str, required=True, help='demographic model without introgression for simulation in Demes YAML format')
    parser.add_argument('--ms-dir', type=str, dest='ms_dir', required=True, help='directory for the ms program for simulation')
    parser.add_argument('--N0', type=int, required=True, help='N0 used in ms simulation')
    parser.add_argument('--nsamp', type=int, required=True, help='sample size (haploid) used in ms simulation')
    parser.add_argument('--nreps', type=int, required=True, help='number of replicates used in ms simulation')
    parser.add_argument('--seeds', type=int, nargs=3, default=None, help='three random seed numbers used in ms simulation; default: None')
    parser.add_argument('--ref-index', type=int, dest='ref_index', required=True, help='index of the reference population in the demographic model (start from 1)')
    parser.add_argument('--ref-size', type=int, dest='ref_size', required=True, help='sample size (haploid) of the reference population')
    parser.add_argument('--tgt-index', type=int, dest='tgt_index', required=True, help='index of the target population in the demographic model (start from 1)')
    parser.add_argument('--tgt-size', type=int, dest='tgt_size', required=True, help='sample size (haploid) of the target population')
    parser.add_argument('--mut-rate', type=float, dest='mut_rate', required=True, help='mutation rate per generation per base')
    parser.add_argument('--rec-rate', type=float, dest='rec_rate', required=True, help='recombination rate per generation per base')
    parser.add_argument('--seq-len', type=int, dest='seq_len', required=True, help='length of simulated sequence')
    parser.add_argument('--snp-num-range', type=int, dest='snp_num_range', nargs=3, required=True, help='range of SNP numbers in ms simulation; the first parameter is the minimum SNP number, the second parameter is the maximum SNP number, the third parameter is the step size')
    parser.add_argument('--output-dir', type=str, dest='output_dir', required=True, help='directory for the output files')
    parser.add_argument('--thread', type=int, default=1, help='number of thread')
    parser.set_defaults(runner=_run_quantile)
    
    # Arguments for threshold subcommand
    parser = subparsers.add_parser('threshold', help='calculate S* thresholds from simulated data')
    _add_score_args(parser)
    parser.add_argument('--sim-data', type=str, dest='sim_data', required=True, help='name of the file containing simulated data for building a generalized additive model')
    parser.add_argument('--recomb-rate', type=float, dest='recomb_rate', default=1e-8, help='a uniform recombination rate used across the genome; default: 1e-8')
    parser.add_argument('--recomb-map', type=str, dest='recomb_map', default=None, help='a recombination map used across the genome; default: None')
    parser.add_argument('--quantile', type=float, required=True, default=0.99, help='a quantile for determining significant S* scores; default: 0.99')
    parser.add_argument('--output', type=str, required=True, help='name of the output file')
    parser.add_argument('--k', type=int, default=8, help='dimension(s) of the bases used to represent the smooth term')
    parser.set_defaults(runner=_run_threshold)

    # Arguments for matchpct subcommand
    parser = subparsers.add_parser('matchrate', help='calculate source match rates in target populations with genomes from source populations')
    parser.add_argument('--vcf', type=str, dest='vcf', required=True, help='name of the VCF file containing genotypes from samples')
    _add_ref_ind_args(parser)
    _add_tgt_ind_args(parser)
    _add_src_ind_args(parser)
    _add_common_args(parser)
    _add_mapped_args(parser)
    _add_score_args(parser)
    parser.set_defaults(runner=_run_match_pct)

    # Arguments for tract subcommand
    parser = subparsers.add_parser('tract', help='get candidate introgressed fragments')
    parser.add_argument('--threshold', type=str, required=True, help='threshold file from `sstar threshold`')
    parser.add_argument('--match-rate', type=str, default=None, dest='match_pct', nargs='+', action=required_length(1,2), help='match rate files from `sstar matchrate`; the first file contains match percents using genomes from the source population 1 (src1), the second file contains match percents using genomes from the source population 2 (src2)')
    parser.add_argument('--output-prefix', type=str, dest='output', required=True, help='prefix for output files')
    parser.add_argument('--diff', type=float, default=0, help='difference between src1 match rates (src1_match_rate) and src2 match rates (src2_match_rate); if src1_match_rate - src2_match_rate > diff, then this fragment is assigned to src1, if src1_match_rate - src2_match_rate < diff, then this fragment is assigned to src2; default: 0')
    parser.set_defaults(runner=_run_tract)

    return top_parser


def main(arg_list=None):
    """
    Description:
        Main entry fo sstar

    Arguments:
        arg_list list: List containing arguments for sstar
    """
    _set_sigpipe_handler()
    parser = _s_star_cli_parser()
    args = parser.parse_args(arg_list)
    args.runner(args)
