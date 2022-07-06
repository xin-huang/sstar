# Get candidate introgressed fragments

### Input

To obtain candidate introgressed fragments from `sstar`, users should provide the output (e.g. [test.tract.threshold](https://github.com/xin-huang/sstar/blob/main/tests/data/test.tract.threshold)) from `sstar threshold`. If no source genome is available, users can use the following command:

	sstar tract --threshold test.tract.threshold --output-prefix test

If one source individual is available, users also can calculate source match rates and output these source match rates and candidate introgressed regions in a single BED file.

	sstar --threshold test.tract.threshold --output-prefix test --match-rate test.tract.src1.match.rate

If source genomes from two different source populations are available, users could provide the output (e.g. [test.tract.src1.match.rate](https://github.com/xin-huang/sstar/blob/main/tests/data/test.tract.src1.match.rate) and [test.tract.src2.match.rate](https://github.com/xin-huang/sstar/blob/main/tests/data/test.tract.src2.match.rate)) from `sstar matchrate`, and use the following command:

	sstar tract --threshold test.tract.threshold --output-prefix test --match-rate test.tract.src1.match.rate test.tract.src2.match.rate

### Output

BED files containing candidate introgressed fragments will be generated, e.g. [test.tract.exp.bed](https://github.com/xin-huang/sstar/blob/main/tests/results/test.tract.exp.bed)

An example for the output is below:

|   |          |          |        |          |
| - | -        | -        | -      | -        |
| 1 | 75030000 | 75080000 | tsk_10 | 0.280797 |

- The first column is the name of the chromosome.
- The second column is the start position.
- The third column is the end position.
- The fourth column is the name of the sample.
- The fifth column is the source match rate.
