# sstar

`sstar` is a Python package for detecting archaic introgression from population genetic data with S\* scores ([Plagnol and Wall 2006](https://doi.org/10.1371/journal.pgen.0020105)), which can be used for detecting introgression not only in humans, but also in other species.

### Requirements

`sstar` works on Unix/Linux operating systems and tested with the following:

- Python 3.8.19
- R 4.1
- Python packages:
	- demes=0.2.3
	- numpy=1.24.4
	- pandas=2.0.3
	- rpy2=3.5.11
	- scikit-allel=1.3.7
	- scipy=1.10.1
- R packages:
	- MASS
	- mgcv
	- stat

### Installation

Users should first install [R](https://cran.r-project.org/) and all the R dependencies listed above. Please ensure the path for the dynamic libraries in R is in the environment variable `$LD_LIBRARY_PATH`. For example, if the dynamic libraries in R are in the path `/usr/local/lib/R/lib/`, users can add this path with the following command:

	export LD_LIBRARY_PATH=/usr/local/lib/R/lib:$LD_LIBRARY_PATH

Then users can install `sstar` with `pip`.

	pip install sstar

Users can also use [mamba](https://github.com/mamba-org/mamba) to create a virtual environment and install `sstar` with this [build-env.yaml](https://github.com/xin-huang/sstar/blob/main/build-env.yaml).

	mamba env create -f build-env.yaml
	conda activate sstar
	export R_LIBS=$CONDA_PREFIX/lib/R/library

Users may verify whether using the correct versions of Python and R under the virtual environment sstar or not with the following commands:

	which python
	which R

And it should return some paths similar to

	/path/to/conda/envs/sstar/bin/python
	/path/to/conda/envs/sstar/bin/R

Users can also use [Apptainer](https://apptainer.org/) to pull a container image that includes `sstar` and `ms`:

    apptainer pull sstar.sif oras://ghcr.io/xin-huang/sstar:1.2.1

Then run `sstar` with:

    apptainer exec sstar.sif sstar

### Citations

If you find `sstar` is useful, please cite

- Plagnol V, Wall JD. 2006. Possible ancestral structure in human populations. *PLoS Genet*. **2**: e105.
- Vernot B, Akey JM. 2014. Resurrecting surviving Neandertal lineages from modern human genomes. *Science* **343**: 1017–1021.
- Vernot B, et al. 2016. Excavating Neandertal and Denisovan DNA from the genomes of Melanesian individuals. *Science* **352**: 235–239.
- Huang X, Kruisz P, Kuhlwilm M. 2022. sstar: A Python package for detecting archaic introgression from population genetic data with *S*\*. *Mol Biol Evol* **39**: msac212.
