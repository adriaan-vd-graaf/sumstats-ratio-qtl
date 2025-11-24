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

This repository includes example data to demonstrate the usage of `sumstats-ratio`.

1.  **Run the example:**

    You can run the tool with the example data files included in the repository.

    ```bash
    sumstats-ratio --met1 sumstats_ratio/example_data/met1.parquet --met2 sumstats_ratio/example_data/met2.parquet --out ratio_output.parquet
    ```

2.  **Expected Output:**

    This command will generate a file named `ratio_output.parquet` in your root directory. The head of the output DataFrame should look like this:

    ```
       pos_name_1 chromosome  base_pair_location other_allele effect_allele    beta_1      se_1  ...      se_2       z_2 effect_allele_frequency_2       maf        se      beta         z
    0  rs55914848         15            26407429            C             T -0.014968  0.016960  ...  0.015712  1.259270                  0.523108  0.474779  0.017916 -0.034754 -1.939824
    1  rs58408126          6           139360300            T             C -0.038087  0.026572  ...  0.024557 -0.472749                  0.113862  0.112691  0.028041 -0.026478 -0.944258
    2    rs337899          7           125887158            A             G  0.022595  0.017129  ...  0.015856  0.703436                  0.402869  0.400439  0.018088  0.011442  0.632564
    3  rs60593995          5            92279193            T             C -0.018380  0.041404  ...  0.038018 -0.635823                  0.044147  0.042964  0.043571  0.005792  0.132938
    4  rs76553537         12            62426616            A             T -0.002239  0.026997  ...  0.024982  1.230782                  0.107424  0.106662  0.028505 -0.032986 -1.157182

    [5 rows x 18 columns]
    ```

3.  **Expected Runtime:**

    The command should take approximately **1-2 seconds** to run on a standard machine.

## How it works

The tool merges the two metabolite summary statistics files based on common SNPs. It then calculates the optimal correlation between the two metabolites that results in a genomic inflation factor (lambda) of 1 for the ratio. This correlation is then used to calculate the beta, standard error, and Z-score for the ratio of the two metabolites.

## System Requirements

-   **Software Version**: This tool has been tested with `sumstats-ratio-qtl` version `0.1.0`.
-   **Operating System**: The software has been tested on macOS. It is expected to be compatible with other Unix-like systems such as Linux.
-   **Python**: Python >=3.13 is required.
-   **Hardware**: No non-standard hardware is required for use.

Installation takes approximately a minute to run. 