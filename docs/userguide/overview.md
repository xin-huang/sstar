# Overview

`sstar` assumes individuals are diploid and variants are bi-allelic in the combined reference and target populations used for scoring. Source populations are used for downstream annotation and match-rate calculation, and should not determine which variants are retained for S* scoring. Please check the data carefully before using `sstar`.

There are five subcommands in `sstar`:

- `score`
- `quantile`
- `threshold`
- `matchrate`
- `tract`

To display help information for each subcommand, users can use:

	sstar subcommand -h

For example:

	sstar score -h

**Note:** In this manual, we define the population without introgressed fragments as the **reference population**, the population receiving introgressed fragments as the **target population**, and the population donating introgressed fragments as the **source population**.

The example commands assume that users have cloned the `sstar` [GitHub repository](https://github.com/xin-huang/sstar) and run the commands from the root directory of the repository.
