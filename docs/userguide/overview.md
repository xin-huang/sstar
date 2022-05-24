# Overview

There are five subcommands in `sstar`:

- score
- quantile
- threshold
- matchrate
- tract

To display help information for each subcommand, users can use

	sstar subcommand -h

For example,

	sstar score -h 

**Note:** `sstar` assumes a variant is bi-allelic among reference, target, and source populations. Please check your data carefully before using `sstar`. In this manual, we define the population without introgressed fragments as the reference population, the population receiving introgressed fragments as the target population, and the population donating introgressed fragments as the source population. Besides, `sstar` assumes an individual is diploid.
