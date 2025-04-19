from sstar.parsers.common_args import add_common_args, add_window_args, add_ref_ind_args, add_tgt_ind_args

def run_score(args):
    from sstar.calc_s_star import calc_s_star

    calc_s_star(
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
    )

def add_score_parser(subparsers):
    parser = subparsers.add_parser("score", help="calculate S* scores from VCF files")
    parser.add_argument(
        "--vcf",
        type=str,
        dest="vcf",
        required=True,
        help="name of the VCF file containing genotypes from samples",
    )
    add_ref_ind_args(parser)
    add_tgt_ind_args(parser)
    add_common_args(parser)
    add_window_args(parser)
    parser.add_argument(
        "--match-bonus",
        type=int,
        dest="match_bonus",
        default=5000,
        help="bonus for matching genotypes of two different variants; default: 5000",
    )
    parser.add_argument(
        "--max-mismatch",
        type=int,
        dest="max_mismatch",
        default=5,
        help="maximum genotype distance allowed; default: 5",
    )
    parser.add_argument(
        "--mismatch-penalty",
        type=int,
        dest="mismatch_penalty",
        default=-10000,
        help="penalty for mismatching genotypes of two different variants; default: -10000",
    )
    parser.set_defaults(runner=run_score)
