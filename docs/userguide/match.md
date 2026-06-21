# match

The `match` command calculates source match rates for inferred tracts.

For each tract, it compares the target sample with all source samples in the input VCF and reports the mean match rate across source individuals. The tract file must contain at least four columns: chromosome, start, end, and sample. Additional columns are ignored. Target and source sample lists must not overlap. For phased target labels, names ending in `_1` or `_2` are interpreted as haplotypes of the corresponding base sample.

Tract coordinates are interpreted as BED-style intervals. A VCF position is included in a tract when:

```
start < POS <= end
```

If a tract contains no variants from the VCF, `match_rate` is written as `NA`.

### Example

Users can calculate the match rate between each inferred tract and a Neanderthal sample with the following command:

```
sstar2 match --vcf examples/data/sstar2.example.biallelic.snps.vcf.gz \
             --tgt examples/data/tgt.samples.list \
             --src examples/data/nean.samples.list \
             --tract-file examples/results/sstar2.example.inferred.tracts.bed \
             --output examples/results/sstar2.example.nean.match.rate.bed
```

The output can be found [here](https://github.com/xin-huang/sstar/blob/main/examples/results/sstar2.example.nean.match.rate.bed).

### Outputs

The output is a tab-separated BED file without a header:

```
chrom  start  end  sample  match_rate
```

### Settings

| Argument | Description |
| - | - |
| `--vcf` | Path to the VCF file. |
| `--tgt` | Path to the target individual list. |
| `--src` | Path to the source individual list. |
| `--tract-file` | Path to the inferred tract BED file from `sstar2 infer`. |
| `--output` | Path to the output match-rate BED file. |
| `--ploidy` | Ploidy used to normalize dosage differences. Default: 2. |
