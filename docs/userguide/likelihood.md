# Calculating posterior probabilities

If genomes from two different source populations are available, `sstar` can estimate posterior probabilities to determine whether an S* haplotype is more likely from the source population 1 or the source population 2. The details can be found in Vernot et al. (2016).

### Input

To estimate posterior probabilities, users need to simulate data under the demographic model without introgression (i.e. the null model), the demographic model with introgression from the source population 1 (i.e. the alternative model 1), and the demographic model with introgression from the source population 2 (i.e. the alternative model 2). For each model, users also need to calculate p-values of archaic match percentages with the simulated genomes from the source population 1 and the source population 2, respectively. The input simulated data (e.g. [sstar.example.simulated.alt1.src1.pval](https://raw.githubusercontent.com/xin-huang/sstar/main/examples/data/simulated_data/sstar.example.simulated.alt1.src1.pval.txt)) should be formatted as the [output](https://sstar.readthedocs.io/en/latest/userguide/pvalue/#output) from `sstar pvalue`. These input file can be compressed by `Gzip`.

For real data, users should provide p-values (e.g. [sstar.example.match.nean.filtered.pval](https://github.com/xin-huang/sstar/blob/main/examples/results/sstar.example.match.nean.filtered.pval.txt) and [sstar.example.match.den.filtered.pval](https://github.com/xin-huang/sstar/blob/main/examples/results/sstar.example.match.den.filtered.pval.txt)) estimated by `sstar pvalue`. A file (e.g. [star.example.threshold.txt](https://raw.githubusercontent.com/xin-huang/sstar/main/examples/results/sstar.example.threshold.txt)) from `sstar threshold` is also needed, which contains significant S\* scores.

Users can estimate the posterior probabilities under different models with the following command:

	sstar prob --sim-null sstar.example.simulated.null.src1.pval.txt sstar.example.simulated.null.src2.pval.txt --sim-alt1 sstar.example.simulated.alt1.src1.pval.txt sstar.example.simulated.alt1.src2.pval.txt --sim-alt2 sstar.example.simulated.alt2.src1.pval.txt sstar.example.simulated.alt2.src2.pval.txt --real-data sstar.example.match.nean.filtered.pval.txt sstar.example.match.den.filtered.pval.txt --threshold sstar.example.threshold.txt --output test.posterior.prob.results --grid-size 10

The expected result above can be found in [test.posterior.prob.exp.results](https://github.com/xin-huang/sstar/blob/main/tests/results/test.posterior.prob.exp.results).

### Output

An example for the output is below:

| chrom | start | end | sample | hap_index | S*_start | S*_end | post_null | post_src1 | post_src2 |
| -     | -     | -   | -      | -         | -        | -      | -         | -         | -         |
| 21 | 15270000 | 15320000 | NA07048 | 1 | 15307909 | 15319809 | 0.99999924 | 4.5e-07 | 3.1e-07 |

The meaning of each column:

- The `chrom` column is the name of the chromosome.
- The `start` column is  the start position of the current window.
- The `end` column is the end position of the current window.
- The `sample` column is the name of the individual.
- The `hap_index` column is the index of the haplotype.
- The `S*_start` column is the start position of the current S* haplotype.
- The `S*_end` column is the end position of the current S* haplotype.
- The `post_null` column is the posterior probability under the null modle (i.e. without introgression).
- The `post_src1` column is the posterior probability under the alternative model 1 (i.e. introgression from the source population 1).
- The `post_src2` column is the posterior probability under the alternative model 2 (i.e. introgression fromt the source population 2).

### Settings

Users can use the argument `--grid-size` to change the grid size for the kernel density estimation. Also, users can change the bandwidth for the kernel density estimation under the null model with the argument `--null-bandwidth`, under the alternative model 1 with the argument `--alt1-bandwidth`, and under the alternative model 2 with the argument `--alt2-bandwidth`. Finally, users can specify a step size with the argument `--step-size` to split prior distributions (uniform distribution `(0, 1)`) into different proportions.
