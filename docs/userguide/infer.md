# infer

```
sstar2 infer --model examples/data/trained.model.onnx --config examples/data/sstar2.config.yaml --feat-file test.score.tsv --pred-file test.pred.tsv --tract-file test.inferred.tracts.bed
```

### Settings

| Argument | Description |
| - | - |
| `--model` | Path to the model file. |
| `--config` | Path to the config file. |
| `--feat-file` | Path to the feature TSV file. |
| `--pred-file` | Path to the prediction TSV file. |
| `--tract-file` | Path to the output BED file. | 
| `--match-bonus` | Bonus for matching genotypes between two variants. Default: `5000`. |
| `--max-mismatch` | Maximum genotype distance allowed before a pair is discarded. Default: `5`. |
| `--mismatch-penalty` | Penalty for mismatching genotypes between two variants. Default: `-10000`. | 
