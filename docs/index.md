# sstar

### Introduction

`sstar` is a Python package for detecting archaic admixture from population genetic data with S* scores (Plagnol and Wall 2006), which can be used for detecting introgression not only in humans, but also in other species, e.g. Kulhwilm et al. (2019). `sstar` re-implements the algorithm and pipeline in [freezing-archer](https://github.com/bvernot/freezing-archer) (Vernot and Akey 2014; Vernot et al. 2016) with [Python 3](https://www.python.org/downloads/), because `freezing-archer` was originally developed with [Python 2](https://www.python.org/doc/sunset-python-2/) that is no longer supported from 2020.

### Requirements

`sstar` was implemented and tested on UNIX/LINUX operating systems with the following:

- Python 3.8/3.9
- R 4.0.0/4.1.1
- Python packages:
	- pandas
	- scikit-allel
	- rpy2
- R packages:
	- MASS
	- mgcv
	- spatstat
	- stat

### Installation

Users should first install [R](https://cran.r-project.org/) and all the R dependencies listed above. Please ensure the path for the dynamic libraries in R is in the environment variable `$LD_LIBRARY_PATH`. For example, if the dynamic libraries in R are in the path `/usr/local/lib/R/lib/`, users can add this path with the following command:

	export LD_LIBRARY_PATH=/usr/local/lib/R/lib/:$LD_LIBRARY_PATH

Then users can install `sstar` with `pip`.

	pip install sstar

Users can also use `conda` to create a virtual environment and install `sstar` with this [conda-env.yaml](https://github.com/xin-huang/sstar/blob/main/examples/snakepipe/conda-env.yaml). To install `conda`, please follow the [instruction](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html). Then users can use the following commands:

	conda env create -f conda-env.yaml
	conda activate sstar
	export R_LIBS=$CONDA_PREFIX/lib/R/library
