# Calculate significant S* score thresholds

### Input

To determine the significance of S\* scores, users should use data from simulations. Users need to provide data from a simulated dataset (e.g. [gravel_asn_scale_60k.simulated.data](https://raw.githubusercontent.com/xin-huang/sstar/main/examples/data/simulated_data/gravel_asn_scale_60k.simulated.data)) for building a generalized additive model to predict the expected S\* score under a demographic model without introgression. This file can be created with the [command](https://sstar.readthedocs.io/en/latest/userguide/quantile/) `sstar quantile` or the `sstar` pipeline [sstar.snake](https://github.com/admixVIE/sstar-analysis/blob/main/workflows/1src/sstar.snake) in [sstar-analysis](https://github.com/admixVIE/sstar-analysis). In the simulation, S\* is calculated for a given number of SNPs per window and a given recombination rate. For each value of number of SNPs (or a grid of values in case both number of SNPs and recombination rate are used), the simulation should be of succifient scale (e.g. simulating 20,000 windows), of the same window size used for calculating S\* in real data (e.g. 50kbp), and with the same number of individuals used in real data.

An example for the simulated dataset (using the software ms) is below:

| S*_score | SNP_number | quantile | log(local_recomb_rate) |
| - | - | - | - |
| 0 | 10 | 0.8 | -7.7500000013365 |

The meaning of each column:

- The S\*_score column is the S\* score at a given quantile (column 3), in a simulated dataset for a given number of SNPs (column 2) and recombination rate (column 4).
- The SNP_number column is the number of SNPs defined to occur within all windows of a simulated dataset for which quantiles and S\* scores are obtained.
- The quantile column is the quantile at which the S\* score in column 1 is observed in the simulated dataset 
- The log(local_recomb_rate) column is the logarithm of the local recombination rate of a simulated dataset.

Users also need to provide the file containing S\* scores estimated from real data (e.g. [test.score.exp.results](https://github.com/xin-huang/sstar/blob/main/tests/results/test.score.exp.results)). Users can use a BED file (e.g. [hum.windows.50k.10k.recomb.map](https://raw.githubusercontent.com/xin-huang/sstar/main/examples/data/real_data/hum.windows.50k.10k.recomb.map)) to specify local recombination rates across the genome. If a window is not found in the recombination map, while it is in the output file from `sstar score`, then this window will be ignored when calculating significant thresholds.

Users can estimate significant S\* score thresholds giving a quantile 0.99 with the following command:

	sstar threshold --score test.score.exp.results --sim-data gravel_asn_scale_60k.simulated.data --recomb-map hum.windows.50k.10k.recomb.map --quantile 0.99 --output test.threshold.results

The expected result above can be found in [test.threshold.exp.results](https://github.com/xin-huang/sstar/blob/main/tests/results/test.threshold.exp.results)

### Output

An example for the output is below:

| chrom | start | end | sample | S*_score | expected_S*_score | local_recomb_rate | quantile | significant |
| -     | -     | -   | -      | -        | -                 | -                 | -        | -           |
| 21    | 9480000 | 9530000 | NA06984 | 48493 | 53263.269232  | 1.29162 | 0.99 | False |
| 21    | 9650000 | 9700000 | NA06984 | 55639 | 53754.835563  | 1.29162 | 0.99 | True  |

The meanings of the first to fifth columns are the same as those in the [output](https://sstar.readthedocs.io/en/latest/userguide/score/#output) from `sstar score`. The meanings of the remaining columns:

- The `expected_S*_score` column is the expected S\* score calculated under the null model (i.e. without introgression).
- The `local_recomb_rate` column is the local recombination rate in the current window.
- The `quantile` column is the quantile used to determing whether the S\* score is significant.
- The `significant` column is an indicator indicating whether the S\* score is significant (i.e. `S*_score` > `expected_S*_score`).

### Settings

In the absence of a recombination map for the tested population, a uniform recombination rate can be used. Users may add the argument `--recomb-rate` instead of a recombination map to estimate thresholds.
