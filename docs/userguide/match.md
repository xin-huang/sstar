# match

The `match` command calculates source match rates for inferred tracts.

For each tract, it compares the target sample with all source samples in the input VCF and reports the mean match rate across source individuals. The tract file must contain at least four columns: chromosome, start, end, and sample. Additional columns are ignored.

Tract coordinates are interpreted as BED-style intervals. A VCF position is included in a tract when:

```text
start < POS <= end
```

```
sstar2 match --vcf examples/sstar2.example.biallelic.snps.vcf.gz \
             --tgt examples/tgt.samples.list \
             --src examples/nean.samples.list \
             --tract-file examples/sstar2.example.inferred.tracts.bed \
             --output examples/sstar2.example.nean.match.rate.bed
```

### Outputs

The output is a tab-separated file without a header:

```text
chrom  start  end  sample  match_rate
```

If a tract contains no variants from the VCF, `match_rate` is written as `NA`.

Target and source sample lists must not overlap. For phased target labels, names ending in `_1` or `_2` are interpreted as haplotypes of the corresponding base sample.

### Settings

| Argument | Description |
| - | - |
| `--vcf` | Path to the VCF file. |
| `--tgt` | Path to the target individual list. |
| `--src` | Path to the source individual list. |
| `--tract-file` | Path to the inferred tract BED file from `sstar2 infer`. |
| `--output` | Path to the output match-rate BED file. |
| `--ploidy` | Ploidy used to normalize dosage differences. Default: 2. |
