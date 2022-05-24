# Calculate S* scores

### Input

To calculate S* scores, users should provide a VCF file containing genotypes from the reference and target populations (e.g. [test.data.vcf](https://github.com/xin-huang/sstar/blob/main/tests/data/test.data.vcf)). Users also need to provide two files containing names of individuals from the reference and target populations (e.g. [test.ref.ind.list](https://github.com/xin-huang/sstar/blob/main/tests/data/test.ref.ind.list) and [test.tgt.ind.list](https://github.com/xin-huang/sstar/blob/main/tests/data/test.tgt.ind.list)) for analysis.

Users can calculate S* scores with the following command:

        sstar score --vcf test.data.vcf --ref test.ref.ind.list --tgt test.tgt.ind.list --output test.score.results

The expected result above can be found in [test.score.exp.results](https://github.com/xin-huang/sstar/blob/main/tests/results/test.score.exp.results).

### Output

An example for the output is below:

| chrom | start | end | sample | S*_score | region_ind_SNP_number | S*_SNP_number | S*_SNPs |
| -     | -     | -   | -      | -        | -                     | -             | -       |
| 21    | 0     | 50000 | ind1 | 51470    | 11                    | 6             | 2309,25354,26654,29724,40809,45079 |

The meaning of each column:

- The `chrom` column is the name of the chromosome.
- The `start` column is the start position of the current window for calculating S* scores.
- The `end` column is the end position of the current window for calculating S* scores.
- The `sample` column is the name of the individual.
- The `S*_score` column is the estimated S* score.
- The `region_ind_SNP_number` column is the number of shared derived variants in the current window between all the individuals from the reference populations and the current individual from the target population.
- The `S*_SNP_number` column is the number of S* SNPs found in the current individual.
- The `S*_SNPs` column is the positions for S* SNPs found in the current individual.

### Settings

By default, `sstar` assumes the reference allele is the ancestral allele and the alternative allele is the dervied allele. Users can use the argument `--anc-allele` with a BED format file (e.g. [test.anc.allele.bed](https://github.com/xin-huang/sstar/blob/main/tests/data/test.anc.allele.bed)) to define the ancestral allele for each variant. If `--anc-allele` is used, then variants without ancestral allele information will be removed. Besides, `sstar` uses a window size with 50,000 bp and step size with 10,000 bp for calculating S\* scores. Users can change these settings with the arguments `--win-len` and `--win-step`. Finally, users can use `--thread` to specifiy numbers of CPUs in order to speed up the calculation.
