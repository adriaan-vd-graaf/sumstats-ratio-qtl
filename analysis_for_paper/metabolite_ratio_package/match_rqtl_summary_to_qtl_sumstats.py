import pandas as pd
import scipy.stats
import numpy as np

def read_from_df_load_gwas_provide_top_hit(df_line):
    ratio_name = df_line.ratio_name
    met1_name = df_line.met1_name
    met2_name = df_line.met2_name

    accession = df_line.ratio_accession

    file_name = f'/Users/avanderg/sgg_unil_postdoc/metabolite_ratio_app/metabolite-ratio-app/data/rQTL_with_mQTL/{accession}.parquet.gz'
    df = pd.read_parquet(file_name)
    df = df.copy()
    df['p_value'] = 2* scipy.stats.norm.sf(np.abs(df['z']))
    df = df.rename(columns={'position': 'base_pair_location'})

    ratio_top_hit = df[df.p_value.min() == df.p_value]
    met1_top_hit = df[df.pval_1.min() == df.pval_1]
    met2_top_hit = df[df.pval_2.min() == df.pval_2]
    print(accession)
    for call, hit, pname in [(ratio_name, ratio_top_hit, 'p_value'), (met1_name, met1_top_hit, 'pval_1'), (met2_name, met2_top_hit, 'pval_2')]:
        print(f'\t{call}, {hit[pname].iloc[0]}, {hit.pos_name.iloc[0]}, {hit.base_pair_location.iloc[0]}')