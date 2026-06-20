# infer

The `infer` command identifies candidate introgressed tracts from genotype data.

It reads preprocessing settings from the `sstar2` configuration file, computes S\* features for genomic windows, applies a trained ONNX model, and writes windows whose observed `S*_score` is greater than the predicted score to a BED file.

```
sstar2 infer --model examples/data/trained.model.onnx \
             --config examples/data/sstar2.config.yaml \
             --feat-file test.score.tsv \
             --pred-file test.pred.tsv \
             --tract-file test.inferred.tracts.bed
```

### Outputs

- `--feat-file`: feature TSV generated from the input genotype data.
- `--pred-file`: prediction TSV containing the observed and predicted S\* scores.
- `--tract-file`: BED file of inferred tracts.

The BED file is tab-separated and has no header:

```text
chrom  start  end  sample
```

The start coordinate is written in BED-style 0-based format.

### Settings

| Argument | Description |
| - | - |
| `--model` | Path to the trained model file. |
| `--config` | Path to the `sstar2` configuration YAML file. |
| `--feat-file` | Path to the feature TSV file. |
| `--pred-file` | Path to the prediction TSV file. |
| `--tract-file` | Path to the output BED file. | 
| `--match-bonus` | Bonus for matching genotypes between two variants. Default: `5000`. |
| `--max-mismatch` | Maximum genotype distance allowed before a pair is discarded. Default: `5`. |
| `--mismatch-penalty` | Penalty for mismatching genotypes between two variants. Default: `-10000`. | 
