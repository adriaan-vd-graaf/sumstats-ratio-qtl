import pandas as pd
import numpy as np
from sumstats_ratio import main
import tempfile

def test_calculate_ratio():
    # Define the paths to the test files
    met1_met_2_ratio_file = "tests/test_data/paper_validation_files/GCST90200793.parquet.gz"

    full_df = pd.read_parquet(met1_met_2_ratio_file)

    full_df = full_df[['pos_name', 'chromosome', 'position', 'reference_allele', 'effect_allele',
                       'beta_1', 'se_1', 'z_1', 'effect_allele_frequency_1',
                       'beta_2', 'se_2', 'z_2', 'effect_allele_frequency_2',
                       'beta_built', 'se_built_correlation_unit_lambda', 'z_built_correlation_unit_lambda'
                       ]]
    full_df = full_df.rename(columns= {'position': 'base_pair_location',
                                        'reference_allele': 'other_allele',
                                        'beta_built': 'beta',
                                        'se_built_correlation_unit_lambda': 'se',
                                        'z_built_correlation_unit_lambda': 'z'
                                       }
                             )

    full_df = full_df.dropna()
    full_df = full_df.drop_duplicates(subset=['base_pair_location', 'chromosome', 'effect_allele', 'other_allele', ])


    met_1_df = (full_df[['pos_name', 'chromosome', 'base_pair_location', 'other_allele', 'effect_allele', 'beta_1', 'se_1', 'z_1', 'effect_allele_frequency_1']]
                .copy().rename(columns={'position': 'base_pair_location',
                                        'beta_1': 'beta',
                                        'se_1': 'se',
                                        'z_1':'z',
                                        'effect_allele_frequency_1': 'effect_allele_frequency'
                                        })
                )

    met_2_df = (full_df[['pos_name', 'chromosome', 'base_pair_location', 'other_allele', 'effect_allele', 'beta_2', 'se_2', 'z_2',  'effect_allele_frequency_2']]
                .copy().rename(columns={'position': 'base_pair_location',
                                        'reference_allele': 'other_allele',
                                        'beta_2': 'beta',
                                        'se_2': 'se',
                                        'z_2': 'z',
                                        'effect_allele_frequency_2': 'effect_allele_frequency'})
                )

    print(f'{met_1_df.shape=}')
    print(f'{met_2_df.shape=}')

    with tempfile.TemporaryFile() as tempfile1, tempfile.TemporaryFile() as tempfile2:
        met_1_df.to_parquet(tempfile1)
        met_2_df.to_parquet(tempfile2)

        result_df = main.calculate_ratio(tempfile1, tempfile2, maf_threshold=-1.0)
        print(result_df.shape)

    # 1. Check that the returned object is a pandas DataFrame
    assert isinstance(result_df, pd.DataFrame)

    # 2. Check that the DataFrame has the expected columns
    expected_columns = ['se', 'beta', 'z']
    for col in expected_columns:
        assert col in result_df.columns

    # 3. Assertions
    # Compare result_df values ('_calc') with your manual calculations, NOT the file's columns
    assert np.allclose(result_df['beta'],full_df['beta'])
    assert np.allclose(result_df['se'], full_df['se'])
    assert np.allclose(result_df['z'], full_df['z'])