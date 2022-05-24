# Calculate source match rates

### Input

To calculate source match rates, users should provide a VCF file containing genotypes from the reference, target, and source populations (e.g. [test.match.rate.data.vcf](https://github.com/xin-huang/sstar/blob/main/tests/data/test.match.rate.data.vcf)). Users also need to provide three files containing names of individuals from the reference, target and source populations (e.g. [ref.ind.list](https://github.com/xin-huang/sstar/blob/main/examples/data/ind_list/ref.ind.list), [tgt.ind.list](https://github.com/xin-huang/sstar/blob/main/examples/data/ind_list/tgt.ind.list) and [nean.ind.list](https://github.com/xin-huang/sstar/blob/main/examples/data/ind_list/nean.ind.list)) for analysis. The file (e.g. [test.match.rate.score.exp.results](https://github.com/xin-huang/sstar/blob/main/tests/results/test.match.rate.score.exp.results)) containing S* scores from `sstar score` is also required.

Users can calculate source match rates with the following command:

	sstar matchrate --vcf test.match.rate.data.vcf --ref ref.ind.list --tgt tgt.ind.list --src nean.ind.list --score test.match.rate.score.exp.results --output test.match.rate.results

The expected result above can be found in [test.match.rate.exp.results](https://github.com/xin-huang/sstar/blob/main/tests/results/test.match.rate.exp.results).

### Output

An example for the output is below:

| chrom | start | end | sample | match_rate | src_sample |
|    -  | -     | -   | -      | -          | -          | 
| 21 | 9400000 | 9450000 | NA06986 | 0.0454545 | Nean | 

The meanings of the first to fourth columns are the same as those in the [output](https://sstar.readthedocs.io/en/latest/userguide/score/#output) from `sstar score`. The meanings of the remaining columns:

- The `match_rate` column is the source match rate on the current region.
- The `src_sample` column is the name of the individual from the source population for calculating the source match percentage.

### Settings

By default, `sstar` assumes the reference allele is the ancestral allele and the alternative allele is the dervied allele. Users can use the argument `--anc-allele` with a BED format file (e.g. [test.anc.allele.bed](https://github.com/xin-huang/sstar/blob/main/tests/data/test.anc.allele.bed)) to define the ancestral allele for each variant. If `--anc-allele` is used, then variants without ancestral allele information will be removed. Users also can provide a BED file (e.g. [test.mapped.region.bed](https://github.com/xin-huang/sstar/blob/main/tests/data/test.mapped.region.bed)) defining non-overlapping mapped regions with the argument `--mapped-region`. Finally, users can use `--thread` to specifiy numbers of CPUs in order to speed up the calculation. If users do not specify numbers of CPUs.
