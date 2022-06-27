# Calculate quanitles of expected S*

`sstar` provides a function `quantile` to estimate different quantiles of expected S\* scores from simualted data without introgression. If users need to estimate quantiles of expected S\* with large scale simulation, we recommend users follow our `sstar` pipelines in [sstar-analysis](https://github.com/admixVIE/sstar-analysis), for example, [sstar.snake](https://github.com/admixVIE/sstar-analysis/blob/main/workflows/1src/sstar.snake).

### Input

Users need to install the `ms` program from [Hudson Lab](https://home.uchicago.edu/~rhudson1/source/mksamples.html) and provide a demographic model (e.g. BonoboGhost_4K19_no_introgression.yaml[](https://github.com/admixVIE/sstar-analysis/blob/main/config/simulation/models/BonoboGhost_4K19_no_introgression.yaml)) in [Demes](https://popsim-consortium.github.io/demes-spec-docs/main/introduction.html) YAML format. Then users can use the following command:

	sstar quantile --model BonoboGhost_4K19_no_introgression.yaml --ms-dir ./msdir/ --N0 1000 --nsamp 22 --nreps 20000 --ref-index 4 --ref-size 20 --tgt-index 3 --tgt-size 2 --mut-rate 1.2e-8 --rec-rate 0.7e-8 --seq-len 40000 --snp-num-range 25 30 5 --output-dir quantiles --thread 2

The meaning of each option can be found using `sstar quantile -h`.

The expected result above can be found in [test.quantile.exp.summary](https://github.com/xin-huang/sstar/blob/main/tests/results/test.quantile.exp.summary).

### Output

An example for the output is below:

| S*_score | SNP_num | quantile | log(local_recomb_rate) |
| -        | -       | -        | -                      |
| 62168.0  | 25      | 0.5      | -8.154901959985743     |

The meaning of each column:

- The `S*_score` column is the S\* score.
- The `SNP_number` column is the number of SNPs.
- The `quantile` column is the quantile of the expected S\*_score in the simulated dataset.
- The `rec_rate` column is the logarithm of the local recombination rate.
