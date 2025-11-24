import pandas as pd
import numpy as np
from scipy.optimize import root_scalar


def eaf_to_maf(eaf):
    """A function to compute minor allele frequency."""
    if eaf <= 0.5:
        return eaf
    else:
        return 1 - eaf


def objective(correlation, df):
    """Objective function for root_scalar, now accepting a dataframe."""
    beta_built_ratio = df['beta_1'] - df['beta_2']
    se_built_ratio = df['se_1']**2 + df['se_2']**2 - 2 * correlation * df['se_1'] * df['se_2']
    se_built_ratio[se_built_ratio < 0] = 0  # Handle potential floating point errors
    se_built_ratio = np.sqrt(se_built_ratio)

    z_built_ratio = np.divide(beta_built_ratio, se_built_ratio,
                              out=np.zeros_like(beta_built_ratio),
                              where=se_built_ratio != 0)

    lamb = np.median(z_built_ratio**2) / 0.454936423119572
    return lamb - 1


def calculate_ratio(met1_file, met2_file, maf_threshold=0.01):
    """
    Calculates the summary statistics of a ratio of two metabolites.

    :param met1_file: Path to the first metabolite summary statistics file (numerator).
    :param met2_file: Path to the second metabolite summary statistics file (denominator).
    :param maf_threshold: Minor allele frequency threshold.
    :return: A pandas DataFrame with the summary statistics of the ratio.
    """
    met1 = pd.read_parquet(met1_file)
    met2 = pd.read_parquet(met2_file)

    ratio_sumstat = pd.merge(met1, met2, suffixes=['_1', '_2'], how='inner',
                               on=['base_pair_location', 'chromosome', 'effect_allele', 'other_allele',])


    ratio_sumstat['maf'] = ratio_sumstat['effect_allele_frequency_1'].apply(eaf_to_maf)


    ratio_sumstat['z_1'] = ratio_sumstat['beta_1'] / ratio_sumstat['se_1']
    ratio_sumstat['z_2'] = ratio_sumstat['beta_2'] / ratio_sumstat['se_2']
    
    result = root_scalar(objective, args=(ratio_sumstat,), bracket=[-1, 1], method='brentq')
    correlation = result.root

    ratio_sumstat = ratio_sumstat[ratio_sumstat.maf > maf_threshold]


    se_calc = np.sqrt(
        ratio_sumstat['se_1'] ** 2 + ratio_sumstat['se_2'] ** 2 - 2 * correlation * ratio_sumstat['se_1'] *
        ratio_sumstat['se_2'])
    ratio_sumstat['se'] = se_calc
    ratio_sumstat['beta'] = ratio_sumstat['beta_1'] - ratio_sumstat['beta_2']
    ratio_sumstat['z'] = ratio_sumstat['beta'] / ratio_sumstat['se']

    return ratio_sumstat
