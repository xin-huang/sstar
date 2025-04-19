from sstar.parsers.argument_validation import required_length

def run_tract(args):
    from sstar.get_tract import get_tract

    get_tract(
        threshold_file=args.threshold,
        match_pct_files=args.match_pct,
        output_prefix=args.output,
        diff=args.diff,
    )

def add_tract_parser(subparsers):
    parser = subparsers.add_parser("tract", help="get candidate introgressed fragments")
    parser.add_argument(
        "--threshold",
        type=str,
        required=True,
        help="threshold file from `sstar threshold`",
    )
    parser.add_argument(
        "--match-rate",
        type=str,
        default=None,
        dest="match_pct",
        nargs="+",
        action=required_length(1, 2),
        help="match rate files from `sstar matchrate`; the first file contains match percents using genomes from the source population 1 (src1), the second file contains match percents using genomes from the source population 2 (src2)",
    )
    parser.add_argument(
        "--output-prefix",
        type=str,
        dest="output",
        required=True,
        help="prefix for output files",
    )
    parser.add_argument(
        "--diff",
        type=float,
        default=0,
        help="difference between src1 match rates (src1_match_rate) and src2 match rates (src2_match_rate); if src1_match_rate - src2_match_rate > diff, then this fragment is assigned to src1, if src1_match_rate - src2_match_rate < diff, then this fragment is assigned to src2; default: 0",
    )
    parser.set_defaults(runner=run_tract)


