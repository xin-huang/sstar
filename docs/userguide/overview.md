# Overview

`sstar2` assumes biallelic genotype data with derived alleles coded as `1` and supports both phased and unphased data.

There are four subcommands in `sstar2`:

- `train`
- `infer`
- `match`
- `assign`

To display help information for each subcommand, users can use:

```
sstar2 <subcommand> -h
```

For example:

```
sstar2 train -h
```

If you need further help, such as reporting a bug or suggesting a feature, please open an [issue](https://github.com/xin-huang/sstar/issues).

**Note:** In this manual, we define the population without introgressed fragments as the **reference population**, the population receiving introgressed fragments as the **target population**, and the population donating introgressed fragments as the **source population**.

The example commands assume that users have cloned the `sstar` [GitHub repository](https://github.com/xin-huang/sstar) and run the commands from the root directory of the repository:

```
git clone https://github.com/xin-huang/sstar
cd sstar
```
