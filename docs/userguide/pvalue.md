# Calculating p-values for source match percentages

`sstar` calculates source match percentages as defined in Vernot et al. (2016).

### Input

To calculate p-values for source match percentages, users should provide a VCF file containing genotypes from the reference, target, and source populations (e.g. [test.pvalue.data.vcf](https://github.com/xin-huang/sstar/blob/main/tests/data/test.cal.pvalue.ref.data.vcf)). Users also need to provide three files containing names of individuals from the reference, target and source populations (e.g. [ref.ind.list](https://github.com/xin-huang/sstar/blob/main/examples/data/ind_list/ref.ind.list), [tgt.ind.list](https://github.com/xin-huang/sstar/blob/main/examples/data/ind_list/tgt.ind.list) and [nean.ind.list](https://github.com/xin-huang/sstar/blob/main/examples/data/ind_list/nean.ind.list)) for analysis.

Another two files are required. One is the file (e.g. [test.pvalue.score.exp.results](https://github.com/xin-huang/sstar/blob/main/tests/results/test.pvalue.score.exp.results)) containing S* scores from `sstar score`. The other is the file (e.g. [test.pvalue.ref.match.pct](https://raw.githubusercontent.com/xin-huang/sstar/main/tests/data/test.pvalue.ref.match.pct)) containing match percentages in the reference population from `sstar rmatch`.

Users can calculate p-values with the following command:

	sstar pvalue --vcf test.pvalue.data.vcf --ref-ind ref.ind.list --tgt-ind tgt.ind.list --src-ind nean.ind.list --score test.pvalue.score.exp.results --ref-match-pct test.pvalue.ref.match.pct --output test.pvalue.results

The expected result above can be found in [test.pvalue.exp.results](https://github.com/xin-huang/sstar/blob/main/tests/results/test.pvalue.exp.results).

### Output

An example for the output is below:

| chrom | start | end | sample | p-value | src_sample | hap_index | S*_start | S*_end | S*_SNP_num | hap_dSNV_num | hap_len | hap_mapped_len | hap_match_num | hap_tot_num | hap_dSNP_per_site_num | hap_S*_match(%) | hap_num_match_ref |
| -  | -       | -       | -       | -  | -    | - | -       | -       | - | - | -     | -     | -   | - | -        | -   | - |
| 21 | 9400000 | 9450000 | NA06986 | NA | Nean | 2 | 9432456 | 9442545 | 2 | 7 | 10089 | 10000 | 0.0 | 7 | 0.428571 | 0.0 | 0 |

The meanings of the first to fourth columns are the same as those in the [output](https://sstar.readthedocs.io/en/latest/userguide/score/#output) from `sstar score`. The meanings of the remaining columns:

- The `p-value` column is the p-value of the source match percentage on the current S\* haplotype.
- The `src_sample` column is the name of the individual from the source population for calculating the source match percentage.
- The `hap_index` column is the index of the haplotype, i.e., the first or second haplotype of the current individual, because `sstar` assumes an individual is diploid.
- The `S*_start` column is the start position of the current S\* haplotype. We define an S\* haplotype using S\* SNPs found in a window. We further require an S\* haplotype should contain three SNPs at least and it also should both start and end with derived alleles. For example, the example above has three S\* SNPs at the positions 9420130, 9432456, and 9442545. We could define an S\* haplotype from the position 9420130 to the position 9442545. However, there are no derived alleles at these two positions and only one SNP is left on the haplotype. Therefore, we find no S\* haplotype on the current window. 
- The `S*_end` column is the end position of the current S\* haplotype.
- The `S*_SNP_num` column is the number of S\* SNPs on the current S\* haplotype.
- The `hap_dSNV_num` column is the number of SNPs with derived alleles on the current S\* haplotype.
- The `hap_len` column is the length of the current S\* haplotype.
- The `hap_mapped_len` column is the length of the mapped region on the current S\* haplotype.
- The `hap_match_num` column is the number of SNPs with derived alleles both on the current S\* haplotype and the source genomes.
- The `hap_tot_num` column is the number of SNPs with derived alleles either on the current S\* haplotype or the source genomes.
- The `hap_dSNP_per_site_num` column is the average number of SNPs with derived alleles per site per haplotype within the range of the current S\* haplotype, which is estimated using all the individuals from the target population.
- The `hap_S*_match(%)` column is the source match percentage of the current S\* haplotype.
- The `hap_num_match_ref` column is the number of haplotypes found in the reference population, which match up haplotype statistics of the current S\* haplotype.

### Settings

By default, `sstar` assumes the reference allele is the ancestral allele and the alternative allele is the dervied allele. Users can use the argument `--anc-allele` with a BED format file (e.g. [test.anc.allele.bed](https://github.com/xin-huang/sstar/blob/main/tests/data/test.anc.allele.bed)) to define the ancestral allele for each variant. If `--anc-allele` is used, then variants without ancestral allele information will be removed. Users also can provide a BED file (e.g. [test.mapped.region.bed](https://github.com/xin-huang/sstar/blob/main/tests/data/test.mapped.region.bed)) defining non-overlapping mapped regions with the argument `--mapped-region`. Finally, users can use `--thread` to specifiy numbers of CPUs in order to speed up the calculation. If users do not specify numbers of CPUs, `sstar` will use all CPUs for calculating p-values.
