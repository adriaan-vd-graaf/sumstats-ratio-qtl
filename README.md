# Sumstats Ratio QTL

A command-line tool to calculate summary statistics for the ratio of two metabolites from their individual summary statistics files.

## Installation

### With `uv`

1. Create a virtual environment:
   ```bash
   python -m venv .venv
   source .venv/bin/activate
   ```
2. Install `uv`:
   ```bash
   pip install uv
   ```
3. Install the required packages and the tool itself with `uv`:
    ```bash
    uv pip install -r requirements.txt
    uv pip install .
    ```

## Usage

You can use the `sumstats-ratio` command to calculate the ratio of two metabolite summary statistics files.

### Arguments

*   `--met1`: Path to the first metabolite summary statistics file (numerator). (Required)
*   `--met2`: Path to the second metabolite summary statistics file (denominator). (Required)
*   `--out`: Path for the output file (e.g., `ratio_sumstats.parquet`). (Required)
*   `--maf`: The minor allele frequency threshold. (Optional, default: 0.01)

### Example

```bash
sumstats-ratio --met1 /path/to/metabolite1.parquet --met2 /path/to/metabolite2.parquet --out ./ratio_output.parquet --maf 0.01
```

This command will generate a new file named `ratio_output.parquet` containing the summary statistics for the ratio of the two input metabolites.

## How it works

The tool merges the two metabolite summary statistics files based on common SNPs. It then calculates the optimal correlation between the two metabolites that results in a genomic inflation factor (lambda) of 1 for the ratio. This correlation is then used to calculate the beta, standard error, and Z-score for the ratio of the two metabolites.

## System Requirements

-   **Software Version**: This tool has been tested with `sumstats-ratio-qtl` version `0.1.0`.
-   **Operating System**: The software has been tested on macOS. It is expected to be compatible with other Unix-like systems such as Linux.
-   **Python**: Python >=3.13 is required.
-   **Hardware**: No non-standard hardware is required for use.

Installation takes approximately a minute to run. 