# Configuration

For the `train` and `infer` commands, `sstar2` requires a [YAML](https://en.wikipedia.org/wiki/YAML) configuration file specifying parameters for simulation, preprocessing, and the machine learning model ([GradientBoostingRegressor](https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.GradientBoostingRegressor.html)). An example configuration file can be found [here](https://github.com/xin-huang/sstar/blob/docs/examples/data/sstar2.config.template.yaml) and is shown below:

```
simulation:
  nrep: 100
  nref: 10
  ntgt: 10
  ref_id: "Western"
  tgt_id: "Bonobo"
  seq_len: 40000
  mut_rate: 1.2e-8
  rec_rate: 0.7e-8
  ploidy: 2
  is_phased: False
  nfeature: 100
  is_shuffled: False
  nprocess: 2
  seed: 4836

preprocessing:
  vcf_file: "examples/data/example.vcf"
  chr_name: "1"
  ref_ind_file: "examples/data/ref.ind.list"
  tgt_ind_file: "examples/data/tgt.ind.list"
  win_len: 40000
  win_step: 10000
  is_phased: False
  nprocess: 2
  ploidy: 2

model:
  params:
    loss: "quantile"
    alpha: 0.999
    n_estimators: 200
    max_depth: 3
    random_state: 4836
```

The configuration file has three top-level sections: `simulation`, `preprocessing`, and `model`.

In the `simulation` section:

| Parameter | Description |
|---|---|
| `nrep` | Number of simulations run in each batch. Additional batches are simulated until at least `nfeature` genomic windows are obtained. |
| `nref` | Number of reference individuals. |
| `ntgt` | Number of target individuals. |
| `ref_id` | Population label used for the reference population in the simulated data. |
| `tgt_id` | Population label used for the target population in the simulated data. |
| `seq_len` | Length of the simulated sequence. |
| `mut_rate` | Mutation rate used in simulation. |
| `rec_rate` | Recombination rate used in simulation. |
| `ploidy` | Ploidy of the simulated individuals. |
| `is_phased` | Whether the simulated data are treated as phased. |
| `nfeature` | Minimum number of genomic windows to generate for model training. |
| `is_shuffled` | Whether the generated features are shuffled before model training. |
| `nprocess` | Number of processes used for simulation. |
| `seed` | Random seed for reproducibility. |

In the `preprocessing` section:

| Parameter | Description |
|---|---|
| `vcf_file` | Path to the input VCF file for inference. |
| `chr_name` | Chromosome name to analyze. |
| `ref_ind_file` | Path to the file listing reference individuals. |
| `tgt_ind_file` | Path to the file listing target individuals. |
| `win_len` | Window length used to divide the genome. |
| `win_step` | Step size between adjacent windows. |
| `is_phased` | Whether the input VCF is treated as phased. |
| `nprocess` | Number of processes used for preprocessing. |
| `ploidy` | Ploidy of the individuals in the input VCF. |

In the `model` section, parameters under `params` are passed to `GradientBoostingRegressor`:

| Parameter | Description |
|---|---|
| `loss` | Loss function used by the gradient boosting model. Must be set to `quantile`. |
| `alpha` | Quantile level used |
| `n_estimators` | Number of boosting stages. |
| `max_depth` | Maximum depth of each regression tree. |
| `random_state` | Random seed used by the model. |

For other available model parameters, see the official scikit-learn documentation for [GradientBoostingRegressor](https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.GradientBoostingRegressor.html).
