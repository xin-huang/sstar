def sstar.parsers.argument_validation import required_length

def run_quantile(args):
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
    )


def add_quantile_parser(subparsers):
    parser = subparsers.add_parser(
        "quantile",
        help="calculate quantiles for S* scores from simulated data without introgression",
    )
    parser.add_argument(
        "--model",
        type=str,
        required=True,
        help="demographic model without introgression for simulation in Demes YAML format",
    )
    parser.add_argument(
        "--ms-dir",
        type=str,
        dest="ms_dir",
        required=True,
        help="directory for the ms program for simulation",
    )
    parser.add_argument(
        "--N0", type=int, required=True, help="N0 used in ms simulation"
    )
    parser.add_argument(
        "--nsamp",
        type=int,
        required=True,
        help="sample size (haploid) used in ms simulation",
    )
    parser.add_argument(
        "--nreps",
        type=int,
        required=True,
        help="number of replicates used in ms simulation",
    )
    parser.add_argument(
        "--seeds",
        type=int,
        nargs=3,
        default=None,
        help="three random seed numbers used in ms simulation; default: None",
    )
    parser.add_argument(
        "--ref-index",
        type=int,
        dest="ref_index",
        required=True,
        help="index of the reference population in the demographic model (start from 1)",
    )
    parser.add_argument(
        "--ref-size",
        type=int,
        dest="ref_size",
        required=True,
        help="sample size (haploid) of the reference population",
    )
    parser.add_argument(
        "--tgt-index",
        type=int,
        dest="tgt_index",
        required=True,
        help="index of the target population in the demographic model (start from 1)",
    )
    parser.add_argument(
        "--tgt-size",
        type=int,
        dest="tgt_size",
        required=True,
        help="sample size (haploid) of the target population",
    )
    parser.add_argument(
        "--mut-rate",
        type=float,
        dest="mut_rate",
        required=True,
        help="mutation rate per generation per base",
    )
    parser.add_argument(
        "--rec-rate",
        type=float,
        dest="rec_rate",
        required=True,
        help="recombination rate per generation per base",
    )
    parser.add_argument(
        "--seq-len",
        type=int,
        dest="seq_len",
        required=True,
        help="length of simulated sequence",
    )
    parser.add_argument(
        "--snp-num-range",
        type=int,
        dest="snp_num_range",
        nargs=3,
        required=True,
        help="range of SNP numbers in ms simulation; the first parameter is the minimum SNP number, the second parameter is the maximum SNP number, the third parameter is the step size",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        dest="output_dir",
        required=True,
        help="directory for the output files",
    )
    parser.add_argument("--thread", type=int, default=1, help="number of thread")
    parser.set_defaults(runner=run_quantile)
