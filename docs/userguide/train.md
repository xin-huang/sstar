# train

The `train` command trains an ONNX model for S\*-based archaic introgression detection.

It reads a DEMES demographic model and an `sstar2` configuration file. If the training feature table is not already present, `sstar2 train` first simulates training data and writes the feature table. It then trains a `GradientBoostingRegressor` model to predict `S*_score` from `Region_ind_SNP_number`.

The training feature table is derived from the model output path:

```text
<output-prefix>.training.features.tsv
```

If this file already exists, `sstar2 train` reuses it instead of running simulation again.

```
sstar2 train --demes examples/data/BonoboGhost_4K19_wo_introgression.yaml \
             --config examples/data/sstar2.config.yaml \
             --output test.model.onnx
```

### Outputs

- `--output`: trained ONNX model file.
- `<output-prefix>.training.features.tsv`: simulated training feature table.

### Settings

| Argument | Description |
| - | - |
| `--demes` | Path to the DEMES demographic model file. |
| `--config` | Path to the `sstar2` configuration YAML file. |
| `--output` | Path to the trained model output file. |
| `--match-bonus` | Bonus for matching genotypes between two variants. Default: `5000`. |
| `--max-mismatch` | Maximum genotype distance allowed before a pair is discarded. Default: `5`. |
| `--mismatch-penalty` | Penalty for mismatching genotypes between two variants. Default: `-10000`. |
