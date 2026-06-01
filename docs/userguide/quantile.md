# Calculate quantiles of expected S*

`sstar` provides a `quantile` subcommand to estimate quantiles of expected S\* scores from simulated data **without** introgression.

### Input

Users need to install the `ms` program from [Hudson Lab](https://home.uchicago.edu/~rhudson1/source/mksamples.html) and provide a demographic model (e.g., [BonoboGhost_4K19_no_introgression.yaml](https://github.com/xin-huang/sstar/blob/main/examples/models/BonoboGhost_4K19_no_introgression.yaml)) in [Demes](https://popsim-consortium.github.io/demes-spec-docs/main/introduction.html) YAML format. The reference sample size should match the number of chromosomes in the real data to be analyzed, whereas the target sample size should be `2`, because `sstar` predicts introgressed fragments for each target sample separately. Then users can use the following command:

	sstar quantile --model examples/models/BonoboGhost_4K19_no_introgression.yaml \
                   --ms-dir ext/msdir/ \
                   --N0 1000 \
                   --nsamp 22 \
                   --nreps 20000 \
                   --ref-pop Western \
                   --ref-size 20 \
                   --tgt-pop Central \
                   --tgt-size 2 \
                   --mut-rate 1.2e-8 \
                   --rec-rate 0.7e-8 \
                   --seq-len 40000 \
                   --snp-num-range 25 30 5 \
                   --output-dir quantiles \
                   --thread 2 \
                   --seeds 1 2 3

The expected result above can be found in [test.quantile.exp.summary](https://github.com/xin-huang/sstar/blob/main/tests/results/test.quantile.exp.summary).

If users run `sstar` from [the Apptainer image](https://github.com/xin-huang/sstar/pkgs/container/sstar), `ms` is available at `/opt/ms`.

### Output

An example of the output is below:

| S*_score | SNP_num | quantile | log(local_recomb_rate) |
| -        | -       | -        | -                      |
| 62284.5  | 25      | 0.5      | -8.154901959985743     |

The meaning of each column:

- The `S*_score` column is the S\* score.
- The `SNP_num` column is the number of SNPs.
- The `quantile` column is the quantile of the expected S\* score in the simulated dataset.
- The `log(local_recomb_rate)` column is the base-10 logarithm of the local recombination rate. This column is not present in the output of the intermediate simulation.

### Settings

| Argument | Description |
| - | - |
| `--model` | Path to the demographic model file used for simulation. |
| `--ms-dir` | Path to the directory containing the `ms` executable. |
| `--N0` | Reference effective population size used for `ms` parameter scaling. |
| `--nsamp` | Haploid sample size used in `ms` simulation. Should be equal to the sum of `--ref-size` and `--tgt-size`. |
| `--nreps` | Number of simulation replicates. |
| `--seeds` | Three random seed numbers passed to ms with `-seeds`. |
| `--ref-pop` | Name of the reference population in the demographic model. |
| `--ref-size` | Haploid sample size of the reference population. |
| `--tgt-pop` | Name of the target population in the demographic model. |
| `--tgt-size` | Haploid sample size of the target population. |
| `--mut-rate` | Mutation rate per site per generation. |
| `--rec-rate` | Recombination rate per site per generation. |
| `--seq-len` | Length of the simulated sequence. |
| `--snp-num-range` | Minimum SNP count, maximum SNP count, and step size used for `ms` simulations. |
| `--output-dir` | Path to the output directory. |
| `--quantile-step` | Step size between quantiles from 0.5 to less than 1. Default: `0.005`. |
| `--thread` | Number of threads for multiprocessing. Default: `1`. |
| `--keep-simulated-data` | Keep intermediate simulation directories instead of deleting them after summarizing results. |
| `--phased` | Run `sstar score` on phased haplotypes during simulation scoring. |

To simulate a range of recombination rates, run `sstar quantile` separately for each recombination-rate step and merge the resulting summary tables before using them with `sstar threshold`.
