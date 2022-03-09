# Calculating significant S* score thresholds

### Input

To determing whether an S\* score is statistically significant, users could use simulation under the demographic model without introgression. Users need to provide a simulated dataset (e.g. [gravel_asn_scale_60k.simulated.data](https://raw.githubusercontent.com/xin-huang/sstar/main/examples/data/simulated_data/gravel_asn_scale_60k.simulated.data)) for building up a generalized additive model to predict the expected S\* score when no introgression occurs.

An example for the simulated dataset is below:

| S*_score | SNP_number | log(local_recomb_rate) | quantile |
| - | - | - | - |
| 0 | 10 | -7.7500000013365 | 0.8 |

The meaning of each column:

- The `S*_score` column is the S\* score in the current window.
- The `SNP_number` column is the number of SNPs in the current window.
- The `log(local_recomb_rate)` column is the logarithm of the local recombination rate in the current window.
- The `quantile` column is the quantile of the S\*_score in the simulated dataset.

Users also need to provide the file (e.g. [test.score.exp.results](https://github.com/xin-huang/sstar/blob/main/tests/results/test.score.exp.results)) containing S\* score estimated from real data. Users can use a BED file (e.g. [hum.windows.50k.10k.recomb.map](https://raw.githubusercontent.com/xin-huang/sstar/main/examples/data/real_data/hum.windows.50k.10k.recomb.map)) to specify local recombination rates across the genome.

Users can estimate significant S* score thresholds giving a quantile 0.99 with the following command:

	sstar threshold --score test.score.exp.results --sim-data gravel_asn_scale_60k.simulated.data --recomb-map hum.windows.50k.10k.recomb.map --quantile 0.99 --output test.threshold.results

The expected result above can be found in [test.threshold.exp.results](https://github.com/xin-huang/sstar/blob/main/tests/results/test.threshold.exp.results)

### Output

An example for the output is below:

| chrom | start | end | sample | S*_score | expected_S*_score | local_recomb_rate | quantile | significant |
| -     | -     | -   | -      | -        | -                 | -                 | -        | -           |
| 21    |  0    | 50000 | ind1 | 51470 | -8539.774316 | 1.29162 | 0.99 | True |

The meanings of the first to fifth columns are the same as those in the [output](https://sstar.readthedocs.io/en/latest/userguide/score/#output) from `sstar score`. The meanings of the remaining columns:

- The `expected_S*_score` column is the expected S\* score calculated under the null model (i.e. without introgression).
- The `local_recomb_rate` column is the local recombination rate in the current window.
- The `quantile` column is the quantile used to determing whether the S\* score is significant.
- The `significant` column is an indicator indicating whether the S\* score is significant (i.e. `S*_score` > `expected_S*_score`).

### Settings

Users also could use a uniform recombination rate with the argument `--recomb-rate` instead of a recombination map to estimate thresholds.
