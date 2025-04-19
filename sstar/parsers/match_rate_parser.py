from sstar.parsers.common_args import (
    add_common_args, 
    add_ref_ind_args, 
    add_tgt_ind_args, 
    add_src_ind_args,
    add_mapped_args,
    add_score_args
)

def run_match_rate(args):
    from sstar.calc_match_rate import calc_match_pct

    calc_match_pct(
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

def add_match_rate_parser(subparsers):
    parser = subparsers.add_parser(
        "matchrate",
        help="calculate source match rates in target populations with genomes from source populations",
    )
    parser.add_argument(
        "--vcf",
        type=str,
        dest="vcf",
        required=True,
        help="name of the VCF file containing genotypes from samples",
    )
    add_ref_ind_args(parser)
    add_tgt_ind_args(parser)
    add_src_ind_args(parser)
    add_common_args(parser)
    add_mapped_args(parser)
    add_score_args(parser)
    parser.set_defaults(runner=run_match_rate)

