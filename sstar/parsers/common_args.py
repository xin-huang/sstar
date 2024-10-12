def add_common_args(parser):
    parser.add_argument(
        "--anc-allele",
        type=str,
        dest="anc_allele",
        default=None,
        help="name of the BED format file containing ancestral allele information, otherwise assuming the REF allele is the ancestral allele and the ALT allele is the derived allele; default: None",
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

def add_mapped_args(parser):
    parser.add_argument(
        "--mapped-region",
        type=str,
        dest="mapped_region_file",
        help="name of the BED file containing mapped regions",
    )


def add_window_args(parser):
    parser.add_argument(
        "--win-len",
        type=int,
        dest="win_len",
        default=50000,
        help="length of the window to calculate S* scores; default: 50000",
    )
    parser.add_argument(
        "--win-step",
        type=int,
        dest="win_step",
        default=10000,
        help="step size for moving windows along genomes; default: 10000",
    )

def add_ref_ind_args(parser):
    parser.add_argument(
        "--ref",
        type=str,
        dest="ref_ind",
        required=True,
        help="name of the file containing information for samples without introgression",
    )


def add_tgt_ind_args(parser):
    parser.add_argument(
        "--tgt",
        type=str,
        dest="tgt_ind",
        required=True,
        help="name of the file containing information for samples for detecting introgression",
    )


def add_src_ind_args(parser):
    parser.add_argument(
        "--src",
        type=str,
        dest="src_ind",
        required=True,
        help="name of the file containing information for samples from source populations",
    )


def add_score_args(parser):
    parser.add_argument(
        "--score",
        type=str,
        dest="score",
        required=True,
        default=None,
        help="name of the file containing S* scores calculated by `sstar score`",
    )
