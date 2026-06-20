# train

The `train` command trains an [ONNX](https://onnx.ai/) model for *S*\*-based archaic introgression detection.

It reads a [DEMES](https://popsim-consortium.github.io/demes-spec-docs/main/introduction.html) demographic model **without introgression** and an `sstar2` [configuration](https://xin-huang.github.io/sstar/2.0/userguide/configuration/) file. If the training feature table is not already present, `sstar2 train` first simulates training data and writes the feature table. It then trains a `GradientBoostingRegressor` model to predict `S*_score` from `Region_ind_SNP_number`.

The training feature table is derived from the model output path:

```
<output-prefix>.training.features.tsv
```

If this file already exists, `sstar2 train` reuses it instead of running simulation again.

### Example

Using the example [demographic model](https://github.com/xin-huang/sstar/blob/main/examples/data/HumanNeanderthal_4G21_wo_introgression.yaml) and example [configuration file](https://github.com/xin-huang/sstar/blob/main/examples/data/sstar2.example.config.yaml), users can train a quantile regression model with the following command:

```
sstar2 train --demes examples/data/HumanNeanderthal_4G21_wo_introgression.yaml \
             --config examples/data/sstar2.example.config.yaml \
             --output examples/results/sstar2.example.trained.model.onnx
```

The trained model can be found [here](https://github.com/xin-huang/sstar/blob/main/examples/results/sstar2.example.trained.model.onnx) and the features used to train it are available [here](https://github.com/xin-huang/sstar/blob/main/examples/results/sstar2.example.trained.model.training.features.tsv).

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
