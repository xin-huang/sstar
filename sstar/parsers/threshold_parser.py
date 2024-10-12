from sstar.parsers.common_args import add_score_args

def run_threshold(args):
    from sstar.calc_threshold import calc_threshold

    calc_threshold(
        simulated_data=args.sim_data,
        score_file=args.score,
        recomb_rate=args.recomb_rate,
        recomb_map=args.recomb_map,
        quantile=args.quantile,
        output=args.output,
        k=args.k,
    )

def add_threshold_parser(subparsers):
    parser = subparsers.add_parser(
        "threshold", help="calculate S* thresholds from simulated data"
    )
    add_score_args(parser)
    parser.add_argument(
        "--sim-data",
        type=str,
        dest="sim_data",
        required=True,
        help="name of the file containing simulated data for building a generalized additive model",
    )
    parser.add_argument(
        "--recomb-rate",
        type=float,
        dest="recomb_rate",
        default=1e-8,
        help="a uniform recombination rate used across the genome; default: 1e-8",
    )
    parser.add_argument(
        "--recomb-map",
        type=str,
        dest="recomb_map",
        default=None,
        help="a recombination map used across the genome; default: None",
    )
    parser.add_argument(
        "--quantile",
        type=float,
        required=True,
        default=0.99,
        help="a quantile for determining significant S* scores; default: 0.99",
    )
    parser.add_argument(
        "--output", type=str, required=True, help="name of the output file"
    )
    parser.add_argument(
        "--k",
        type=int,
        default=8,
        help="dimension(s) of the bases used to represent the smooth term",
    )
    parser.set_defaults(runner=run_threshold)

