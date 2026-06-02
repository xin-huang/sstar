# train


```
sstar2 train --demes examples/data/BonoboGhost_4K19_wo_introgression.yaml --config examples/data/sstar2.config.yaml --output test.model.onnx
```

### Settings

| Argument | Description |
| - | - |
| `--demes` | Path to the demes file. |
| `--config` | Path to the config file. |
| `--output` | Path to the output file. |
| `--match-bonus` | Bonus for matching genotypes between two variants. Default: `5000`. |
| `--max-mismatch` | Maximum genotype distance allowed before a pair is discarded. Default: `5`. |
| `--mismatch-penalty` | Penalty for mismatching genotypes between two variants. Default: `-10000`. |
