# Get candidate introgressed fragments

### Input

To obtain candidate introgressed fragments from `sstar`, users should provide the output (e.g. [test.tract.threshold](https://github.com/xin-huang/sstar/blob/main/tests/data/test.tract.threshold)) from `sstar threshold`. If no source genome is available, users can use the following command:

	sstar tract --threshold test.tract.threshold --output-prefix test

If source genomes from two different source populations are available, users could provide the output (e.g. [test.tract.src1.match.rate](https://github.com/xin-huang/sstar/blob/main/tests/data/test.tract.src1.match.rate) and [test.tract.src2.match.rate](https://github.com/xin-huang/sstar/blob/main/tests/data/test.tract.src2.match.rate)) from `sstar matchrate`, and use the following command:

	sstar tract --threshold test.tract.threshold --output-prefix test --match-rate test.tract.src1.match.rate test.tract.src2.match.rate

### Output

Bed files containing candidate introgressed fragments will be generated, e.g. [test.tract.exp.bed](https://github.com/xin-huang/sstar/blob/main/tests/results/test.tract.exp.bed)
