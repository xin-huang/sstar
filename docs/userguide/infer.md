# infer

The `infer` command identifies candidate introgressed tracts from genotype data.

It reads preprocessing settings from the `sstar2` [configuration file](https://xin-huang.github.io/sstar/latest/userguide/configuration/), computes features for genomic windows, applies a trained ONNX model, and writes windows whose observed `S*_score` is greater than the predicted score to a BED file.

### Example

Using the [trained model](https://github.com/xin-huang/sstar/blob/main/examples/results/sstar2.example.trained.model.onnx) from the [train](https://xin-huang.github.io/sstar/latest/userguide/train/) command and the [configuration file](https://github.com/xin-huang/sstar/blob/main/examples/data/sstar2.example.config.yaml), users can infer the introgressed fragments with the following command:

```
sstar2 infer --model examples/results/sstar2.example.trained.model.onnx \
             --config examples/data/sstar2.example.config.yaml \
             --feat-file examples/results/sstar2.example.inference.features.tsv \
             --pred-file examples/results/sstar2.example.pred.tsv \
             --tract-file examples/results/sstar2.example.inferred.tracts.bed
```

The inferred tracts can be found [here](https://github.com/xin-huang/sstar/blob/main/examples/results/sstar2.example.inferred.tracts.bed). Two additional files are generated: [one](https://github.com/xin-huang/sstar/blob/main/examples/results/sstar2.example.inference.features.tsv) records the features for prediction, and [the other](https://github.com/xin-huang/sstar/blob/main/examples/results/sstar2.example.pred.tsv) contains the observed and predicted *S*\* scores.

### Outputs

- `--feat-file`: feature TSV generated from the input genotype data specified in the configuration file.
- `--pred-file`: prediction TSV containing the observed and predicted *S*\* scores.
- `--tract-file`: BED file of inferred tracts.

The BED file is tab-separated and has no header:

```
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
