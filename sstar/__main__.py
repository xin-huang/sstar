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

def _run_likelihood(args):
    from sstar.cal_likelihood import cal_likelihood
    cal_likelihood(sim_null=args.sim_null, sim_alt1=args.sim_alt1, sim_alt2=args.sim_alt2, real_data=args.real_data, output=args.output, grid_size=args.grid_size, bdw_null=args.bdw_null, bdw_alt1=args.bdw_alt1, bdw_alt2=args.bdw_alt2, step_size=args.step_size, threshold=args.threshold)

def _run_threshold(args):
    from sstar.cal_threshold import cal_threshold
    cal_threshold(simulated_data=args.sim_data, score_file=args.score, recomb_rate=args.recomb_rate, recomb_map=args.recomb_map, quantile=args.quantile, output=args.output, k=args.k)

def _run_rmatch(args):
    from sstar.cal_ref_match_pct import cal_ref_match_pct
    cal_ref_match_pct(vcf=args.vcf, ref_ind_file=args.ref_ind, src_ind_file=args.src_ind, anc_allele_file=args.anc_allele, output=args.output, thread=args.thread, win_len=args.win_len, win_step=args.win_step, mapped_region_file=args.mapped_region_file)

def _run_pvalue(args):
    from sstar.cal_pvalue import cal_pvalue
    cal_pvalue(vcf=args.vcf, ref_ind_file=args.ref_ind, tgt_ind_file=args.tgt_ind, src_ind_file=args.src_ind, anc_allele_file=args.anc_allele, output=args.output, thread=args.thread, score_file=args.score, ref_match_pct_file=args.ref_match_pct, mapped_region_file=args.mapped_region_file, low_memory=args.low_memory, mapped_len_esp=args.mapped_len_esp, len_esp=args.len_esp, var_esp=args.var_esp, sfs_esp=args.sfs_esp)

def _add_common_args(parser):
    parser.add_argument('--anc-allele', type=str, dest='anc_allele', default=None, help='name of the BED format file containing ancestral allele information, otherwise assuming the REF allele is the ancestral allele and the ALT allele is the derived allele; default: None')
    parser.add_argument('--output', type=str, required=True, help='name of the output file')
    parser.add_argument('--thread', type=int, default=None, help='number of threads for multiprocessing; default: None')

def _add_mapped_args(parser):
    parser.add_argument('--mapped-region', type=str, dest='mapped_region_file', help='name of the BED file containing mapped regions')

def _add_window_args(parser):
    parser.add_argument('--win-len', type=int, dest='win_len', default=50000, help='length of the window to calculate S* scores; default: 50000')
    parser.add_argument('--win-step', type=int, dest='win_step', default=10000, help='step size for moving windows along genomes; default: 10000')

def _add_ref_ind_args(parser):
    parser.add_argument('--ref-ind', type=str, dest='ref_ind', required=True, help='name of the file containing information for samples without introgression')

def _add_tgt_ind_args(parser):
    parser.add_argument('--tgt-ind', type=str, dest='tgt_ind', required=True, help='name of the file containing information for samples for detecting introgression')

def _add_src_ind_args(parser):
    parser.add_argument('--src-ind', type=str, dest='src_ind', required=True, help='name of the file containing information for samples from source populations')

def _add_score_args(parser):
    parser.add_argument('--score', type=str, dest='score', required=True, default=None, help='name of the file containing S* scores calculated by `sstar score`')

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
    
    # Arguments for assign threshold subcommand
    parser = subparsers.add_parser('threshold', help='calculate S* thresholds from simulated data')
    _add_score_args(parser)
    parser.add_argument('--sim-data', type=str, dest='sim_data', required=True, help='name of the file containing simulated data for building a generalized additive model')
    parser.add_argument('--recomb-rate', type=float, dest='recomb_rate', default=1e-8, help='a uniform recombination rate used across the genome; default: 1e-8')
    parser.add_argument('--recomb-map', type=str, dest='recomb_map', default=None, help='a recombination map used across the genome; default: None')
    parser.add_argument('--quantile', type=float, required=True, default=0.9, help='a quantile for determining significant S* scores; default: 0.9')
    parser.add_argument('--output', type=str, required=True, help='name of the output file')
    parser.add_argument('--k', type=int, default=8, help='dimension(s) of the bases used to represent the smooth term')
    parser.set_defaults(runner=_run_threshold)

    # Arguments for rmatch subcommand
    parser = subparsers.add_parser('rmatch', help='calculate match percents in reference populations with genomes from source populations')
    parser.add_argument('--vcf', type=str, dest='vcf', required=True, help='name of the VCF file containing genotypes from samples')
    _add_ref_ind_args(parser)
    _add_src_ind_args(parser)
    _add_common_args(parser)
    _add_window_args(parser)
    _add_mapped_args(parser)
    parser.set_defaults(runner=_run_rmatch)

    # Arguments for tmatch subcommand
    parser = subparsers.add_parser('pvalue', help='calculate p-values for S* haplotypes in target populations with genomes from source populations')
    parser.add_argument('--vcf', type=str, dest='vcf', required=True, help='name of the VCF file containing genotypes from samples')
    _add_ref_ind_args(parser)
    _add_tgt_ind_args(parser)
    _add_src_ind_args(parser)
    _add_common_args(parser)
    _add_mapped_args(parser)
    _add_score_args(parser)
    parser.add_argument('--ref-match-pct', type=str, dest='ref_match_pct', required=True, help='name of the file containing match percents in reference populations calculated by `sstar rmatch`')
    parser.add_argument('--mapped-len-esp', type=int, default=1000, dest='mapped_len_esp', help='increment of the length of the mapped region')
    parser.add_argument('--len-esp', type=int, default=1000, dest='len_esp', help='increment of the length of the haplotype')
    parser.add_argument('--var-num-esp', type=int, default=1, dest='var_esp', help='increment of the number of derived alleles on the haplotype')
    parser.add_argument('--site-var-num-esp', type=float, default=0.02, dest='sfs_esp', help='increment of the average number of dervied alleles per haplotype per site')
    parser.add_argument('--low-memory', dest='low_memory', action='store_true', help='Determine whether using low memory to query reference match percentages or not')
    parser.set_defaults(runner=_run_pvalue)

    # Arguments for likelihood subcommand
    parser = subparsers.add_parser('prob', help='calculate posterior probabilities')
    parser.add_argument('--sim-null', type=str, nargs=2, dest='sim_null', required=True, help='two files containing p-values for simulated data without introgression (null model)')
    parser.add_argument('--sim-alt1', type=str, nargs=2, dest='sim_alt1', required=True, help='two files containing p-values for simulated data from alternative model 1')
    parser.add_argument('--sim-alt2', type=str, nargs=2, dest='sim_alt2', required=True, help='two files containing p-values for simulated data from alternative model 2')
    parser.add_argument('--real-data', type=str, nargs=2, dest='real_data', required=True, help='two files containing p-values for real data')
    parser.add_argument('--threshold', type=str, required=True, help='name of the file containing significant S* scores from `sstar score`')
    parser.add_argument('--grid-size', type=int, dest='grid_size', default=200, help='a grid size for kernel density estimatin; default: 200')
    parser.add_argument('--null-bandwidth', type=float, dest='bdw_null', default=1.5, help='a bandwidth for kernel density estimation under the null model; default: 1.5')
    parser.add_argument('--alt1-bandwidth', type=float, dest='bdw_alt1', default=6, help='a bandwidth for kernel density estimation under the alternative model 1; default: 6')
    parser.add_argument('--alt2-bandwidth', type=float, dest='bdw_alt2', default=6, help='a bandwidth for kernel density estimation under the alternative model 2; default: 6')
    parser.add_argument('--step-size', type=float, dest='step_size', default=0.005, help='a step size for dividing prior probabilities; default: 0.005')
    parser.add_argument('--output', type=str, required=True, help='name of the output file')
    parser.set_defaults(runner=_run_likelihood)

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
