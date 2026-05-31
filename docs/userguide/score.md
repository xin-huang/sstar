# Calculate S* scores

### Input

To calculate S* scores, users should provide a VCF file containing genotypes from the reference and target populations (e.g. [test.score.data.vcf](https://github.com/xin-huang/sstar/blob/main/tests/data/test.score.data.vcf)). The genotypes in the VCF file should be segregating in the combined reference and target populations used for scoring, and should not contain non-variable sites across these individuals. If additional individuals are included in the VCF file, such as potential source individuals used in later analyses, sites that are variable only in those additional individuals should not be included. Users also need to provide two files containing names of individuals from the reference and target populations (e.g. [test.ref.ind.list](https://github.com/xin-huang/sstar/blob/main/tests/data/test.ref.ind.list) and [test.tgt.ind.list](https://github.com/xin-huang/sstar/blob/main/tests/data/test.tgt.ind.list)) for analysis.

Users can calculate S* scores for unphased data with the following command:

	sstar score --vcf tests/data/test.score.data.vcf \
                --ref tests/data/test.ref.ind.list \
                --tgt tests/data/test.tgt.ind.list \
                --output test.score.unphased.results.tsv

The expected result above can be found in [test.score.unphased.exp.results.tsv](https://github.com/xin-huang/sstar/blob/main/tests/results/test.score.unphased.exp.results.tsv).

To calculate S* scores for phased data, add the --phased argument to the command:

	sstar score --vcf tests/data/test.score.data.vcf \
                --ref tests/data/test.ref.ind.list \
                --tgt tests/data/test.tgt.ind.list \
                --output test.score.phased.results.tsv \
                --phased

The expected result can be found in [test.score.phased.exp.results.tsv](https://github.com/xin-huang/sstar/blob/main/tests/results/test.score.phased.exp.results.tsv).

### Output

An example of the output is below:

| chrom | start | end | sample | S*_score | region_ind_SNP_number | S*_SNP_number | S*_SNPs |
| -     | -     | -   | -      | -        | -                     | -             | -       |
| 21    | 0     | 50000 | ind1 | 51470 | 11 | 6 | 2309,25354,26654,29724,40809,45079 |

The meaning of each column:

- The `chrom` column is the name of the chromosome.
- The `start` column is the start position of the current window for calculating S* scores.
- The `end` column is the end position of the current window for calculating S* scores.
- The `sample` column is the name of the target individual. When S* scores are calculated with the `--phased` argument, suffixes such as `_1` and `_2` are added to distinguish the two haplotypes.
- The `S*_score` column is the estimated S* score.
- The `region_ind_SNP_number` column is the number of shared derived variants in the current window between all the individuals from the reference population and the current target individual or haplotype.
- The `S*_SNP_number` column is the number of S* SNPs found in the current target individual or haplotype.
- The `S*_SNPs` column is the positions for S* SNPs found in the current target individual or haplotype.

### Settings

| Argument | Description |
| - | - |
| `--vcf` | Path to the VCF file containing genotype data. |
| `--ref` | Path to the file containing reference individual IDs. |
| `--tgt` | Path to the file containing target individual IDs. |
| `--output` | Path to the output score file. |
| `--anc-allele` | Path to a BED format file containing ancestral allele information. Variants without ancestral allele information are removed. Default: `None`. |
| `--win-len` | Length of sliding windows in base pairs. Default: `50000`. |
| `--win-step` | Step size for sliding windows in base pairs. Default: `10000`. |
| `--match-bonus` | Bonus for matching genotypes between two variants. Default: `5000`. |
| `--max-mismatch` | Maximum genotype distance allowed before a pair is discarded. Default: `5`. |
| `--mismatch-penalty` | Penalty for mismatching genotypes between two variants. Default: `-10000`. |
| `--thread` | Number of CPUs used for multiprocessing. Default: `1`. |
| `--phased` | Calculate scores on phased haplotypes instead of genotype dosages. |

By default, `sstar` assumes the reference allele is the ancestral allele and the alternative allele is the derived allele. Users can use the argument `--anc-allele` with a BED format file (e.g. [test.anc.allele.bed](https://github.com/xin-huang/sstar/blob/main/tests/data/test.anc.allele.bed)) to define the ancestral allele for each variant. If `--anc-allele` is used, variants without ancestral allele information will be removed.
