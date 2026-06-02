# sstar2

`sstar2` is a Python package for *S*\*-based archaic introgression detection with machine learning.

### Requirements

`sstar2` works on Unix/Linux operating systems and was tested with the following:

- Python 3.12.13
- Python packages:
    - black=26.3.1
    - codecov=2.1.13
    - demes=0.2.3
    - flake8=7.3.0
    - msprime=1.4.1
    - numpy=2.4.5
    - onnx=1.19.1
    - onnxruntime=1.26.0
    - pandas=3.0.3
    - pip=26.1.1
    - pydantic=2.12.4
    - pytest=9.0.3
    - pytest-cov=7.1.0
    - scikit-allel=1.3.13
    - scikit-learn=1.8.0
    - scipy=1.17.1
    - skl2onnx=1.19.1

### Installation

Users can install `sstar2` with `pip`:

	pip install sstar==2.0.0

Users can also use [mamba](https://github.com/mamba-org/mamba) to create a virtual environment and install `sstar2` with this [build-env.yaml](https://github.com/xin-huang/sstar/blob/main/build-env.yaml):

	mamba env create -f build-env.yaml
	conda activate sstar2

Users can also use [Apptainer](https://apptainer.org/) to pull a container image that includes `sstar2`:

    apptainer pull sstar2.sif oras://ghcr.io/xin-huang/sstar:2.0.0

Then run `sstar2` with:

    apptainer exec sstar2.sif sstar2

### Citations

If you find `sstar2` is useful, please cite:

- Koça A, Stöckl A, Chen S, Kuhlwilm M, Huang X. 2026. sstar2: A Python package for *S*\*-based archaic introgression detection with machine learning. bioRxiv: 2026.05.31.729079.
- Huang X, Kruisz P, Kuhlwilm M. 2022. sstar: A Python package for detecting archaic introgression from population genetic data with *S*\*. *Molecular Biology and Evolution* **39**: msac212.
- Vernot B, et al. 2016. Excavating Neandertal and Denisovan DNA from the genomes of Melanesian individuals. *Science* **352**: 235–239.
- Vernot B, Akey JM. 2014. Resurrecting surviving Neandertal lineages from modern human genomes. *Science* **343**: 1017–1021.
- Plagnol V, Wall JD. 2006. Possible ancestral structure in human populations. *PLoS Genetics*. **2**: e105.
