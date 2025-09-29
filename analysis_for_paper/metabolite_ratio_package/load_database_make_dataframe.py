
import sqlite3
import pandas as pd
import numpy as np


def load_database_make_dataframe(db_path: str, silent = True):

    # --- Configuration ---
    # Path to your SQLite database file

    # Significance thresholds
    missense_ld_threshold = 0.8
    eqtl_p_alpha_threshold = 0.05 / 137_102
    pqtl_p_alpha_threshold = 0.05/ 21_601
    gtex_p_alpha_threshold = 0.05 / 554_898
    gwas_p_alpha_threshold = 0.05 / (327*5_095)

    # --- Database Connection ---
    if not silent:
        print("--- Starting Data Loading and Enrichment ---")
    conn = sqlite3.connect(db_path)

    # --- Step 1: Query base data (Ratios, Genomic Regions, Metabolites) ---
    if not silent:
        print("Step 1: Querying base data from normalized tables...")

    base_query = """
    SELECT
        r.ratio_accession,
        r.ratio_name,
        r.max_pgain,
        r.num_associated_loci,
        r.num_novel_loci,
        r.llm_response,
        r.top_finngen_mr,
        r.top_mrlink_mr,
        r.top_gtex_mr_gene,
        r.top_gtex_mr_tissue,
        r.top_eqtl_mr_gene,
        r.top_pqtl_mr_gene,
        r.reaction_distance,
        gr.region_id,
        gr.snp_id AS pos_name,
        gr.chromosome,
        gr.position AS bp_position,
        gr.effect_allele,
        gr.reference_allele as other_allele,
        gr.is_novel,
        gr.is_numerator_driven AS numerator_driven,
        gr.is_denominator_driven AS denominator_driven,
        gr.beta,
        gr.se,
        gr.z,
        gr.log_pval, 
        gr.beta_1, 
        gr.se_1, 
        gr.z_1, 
        gr.beta_2, 
        gr.se_2, 
        gr.z_2, 
        gr.maf,
        gr.log_pgain,
        gr.closest_genes,
        gr.infomap_clusters,
        num_met.name AS num_name,
        num_met.metabolite_accession AS num_accession,
        num_met.hmdb_id AS num_hmdb,
        den_met.name AS den_name,
        den_met.metabolite_accession AS den_accession,
        den_met.hmdb_id AS den_hmdb
    FROM
        genomic_regions gr
    JOIN
        ratios r ON gr.ratio_accession = r.ratio_accession
    LEFT JOIN
        metabolites num_met ON r.numerator_metabolite_accession = num_met.metabolite_accession
    LEFT JOIN
        metabolites den_met ON r.denominator_metabolite_accession = den_met.metabolite_accession;
    """
    df_base = pd.read_sql_query(base_query, conn)
    df_base['position'] = df_base['chromosome'].astype(str) + ':' + df_base['bp_position'].astype(str)

    # --- Step 2: Query and aggregate Missense Genes ---
    if not silent:
        print(f"Step 2: Querying missense variants with LD > {missense_ld_threshold}...")

    missense_query = "SELECT region_id, gene_name FROM missense_variants WHERE ld > ?"
    missense_data = pd.read_sql_query(missense_query, conn, params=(missense_ld_threshold,))
    if not missense_data.empty:
        missense_df = missense_data.groupby('region_id')['gene_name'].apply(lambda x: sorted(list(set(x)))).reset_index(name='missense_genes')
        df_base = pd.merge(df_base, missense_df, on='region_id', how='left')
    else:
        df_base['missense_genes'] = [[] for _ in range(len(df_base))]

    # --- Step 3: Query and aggregate significant MR results as tuples ---
    if not silent:
        print(f"Step 3a: Querying eQTLs with p(alpha) < {eqtl_p_alpha_threshold:.2e}...")
    eqtl_query = """
                 SELECT gr.ratio_accession, cr.exposure_name, cr.alpha, cr.se_alpha, cr.p_alpha
                 FROM colocalization_results cr
                          JOIN genomic_regions gr ON cr.region_id = gr.region_id
                 WHERE cr.source_type = 'eQTL' \
                   AND cr.p_alpha < ? \
                 """
    eqtl_data = pd.read_sql_query(eqtl_query, conn, params=(eqtl_p_alpha_threshold,))
    if not eqtl_data.empty:
        # Create mr_result_eqtl column with full tuple (phenotype, alpha, se, p)
        eqtl_agg_results = eqtl_data.drop_duplicates().groupby('ratio_accession').apply(
            lambda x: sorted(list(zip(x['exposure_name'], x['alpha'], x['se_alpha'], x['p_alpha'])))
        ).reset_index(name='mr_result_eqtl')
        df_base = pd.merge(df_base, eqtl_agg_results, on='ratio_accession', how='left')

        # Create eQTL_genes column with only the phenotype names.
        eqtl_agg_genes = eqtl_data.drop_duplicates().groupby('ratio_accession')['exposure_name'].apply(
            lambda x: sorted(list(set(x)))).reset_index(name='eQTL_genes')
        df_base = pd.merge(df_base, eqtl_agg_genes, on='ratio_accession', how='left')
    else:
        df_base['mr_result_eqtl'] = [[] for _ in range(len(df_base))]
        df_base['eQTL_genes'] = [[] for _ in range(len(df_base))]

    if not silent:
        print(f"Step 3b: Querying pQTLs with p(alpha) < {pqtl_p_alpha_threshold:.2e}...")
    pqtl_query = """
                 SELECT gr.ratio_accession, cr.exposure_name, cr.alpha, cr.se_alpha, cr.p_alpha
                 FROM colocalization_results cr
                          JOIN genomic_regions gr ON cr.region_id = gr.region_id
                 WHERE cr.source_type = 'pQTL' \
                   AND cr.p_alpha < ? \
                 """
    pqtl_data = pd.read_sql_query(pqtl_query, conn, params=(pqtl_p_alpha_threshold,))
    if not pqtl_data.empty:
        # Create mr_result_pqtl column with full tuple (phenotype, alpha, se, p)
        pqtl_agg_results = pqtl_data.drop_duplicates().groupby('ratio_accession').apply(
            lambda x: sorted(list(zip(x['exposure_name'], x['alpha'], x['se_alpha'], x['p_alpha'])))
        ).reset_index(name='mr_result_pqtl')
        df_base = pd.merge(df_base, pqtl_agg_results, on='ratio_accession', how='left')

        # Create pQTL_genes column with only the phenotype names.
        pqtl_agg_genes = pqtl_data.drop_duplicates().groupby('ratio_accession')['exposure_name'].apply(
            lambda x: sorted(list(set(x)))).reset_index(name='pQTL_genes')
        df_base = pd.merge(df_base, pqtl_agg_genes, on='ratio_accession', how='left')
    else:
        df_base['mr_result_pqtl'] = [[] for _ in range(len(df_base))]
        df_base['pQTL_genes'] = [[] for _ in range(len(df_base))]

    if not silent:
        print(f"Step 3c: Querying GTEx with p(alpha) < {gtex_p_alpha_threshold:.2e}...")
    gtex_query = """
                 SELECT gr.ratio_accession, cr.exposure_name, cr.tissue, cr.alpha, cr.se_alpha, cr.p_alpha
                 FROM colocalization_results cr
                          JOIN genomic_regions gr ON cr.region_id = gr.region_id
                 WHERE cr.source_type = 'GTEx' \
                   AND cr.p_alpha < ? \
                 """
    gtex_data = pd.read_sql_query(gtex_query, conn, params=(gtex_p_alpha_threshold,))
    if not gtex_data.empty:
        gtex_data['other_name'] = gtex_data['exposure_name'] + ' (' + gtex_data['tissue'] + ')'

        # Create mr_result_gtex column with full tuple (phenotype (tissue), alpha, se, p)
        gtex_agg_results = gtex_data.drop_duplicates().groupby('ratio_accession').apply(
            lambda x: sorted(list(zip(x['other_name'], x['alpha'], x['se_alpha'], x['p_alpha'])))
        ).reset_index(name='mr_result_gtex')
        df_base = pd.merge(df_base, gtex_agg_results, on='ratio_accession', how='left')

        # Create GTEx_genes column with only the (phenotype, tissue) tuples.
        gtex_agg_genes = gtex_data.drop_duplicates().groupby('ratio_accession')[['exposure_name', 'tissue']].apply(
            lambda x: sorted(list(zip(x['exposure_name'], x['tissue'])))).reset_index(name='GTEx_genes')
        df_base = pd.merge(df_base, gtex_agg_genes, on='ratio_accession', how='left')
    else:
        df_base['mr_result_gtex'] = [[] for _ in range(len(df_base))]
        df_base['GTEx_genes'] = [[] for _ in range(len(df_base))]

    if not silent:
        print(f"Step 3d: Querying Complex Traits (FinnGen) with p(alpha) < {gwas_p_alpha_threshold:.2e}...")
    complex_traits_query = """
                           SELECT gr.ratio_accession, cr.exposure_name, cr.alpha, cr.se_alpha, cr.p_alpha
                           FROM colocalization_results cr
                                    JOIN genomic_regions gr ON cr.region_id = gr.region_id
                           WHERE cr.source_type IN ('finngen') \
                             AND cr.p_alpha < ? \
                           """
    complex_traits_data = pd.read_sql_query(complex_traits_query, conn, params=(gwas_p_alpha_threshold,))
    if not complex_traits_data.empty:
        # Create mr_result_finngen column with full tuple (phenotype, alpha, se, p)
        complex_traits_agg_results = complex_traits_data.drop_duplicates().groupby('ratio_accession').apply(
            lambda x: sorted(list(zip(x['exposure_name'], x['alpha'], x['se_alpha'], x['p_alpha'])))
        ).reset_index(name='mr_result_finngen')
        df_base = pd.merge(df_base, complex_traits_agg_results, on='ratio_accession', how='left')

        # Create complex_traits column with only the phenotype names.
        complex_traits_agg_genes = complex_traits_data.drop_duplicates().groupby('ratio_accession')[
            'exposure_name'].apply(
            lambda x: sorted(list(set(x)))).reset_index(name='complex_traits')
        df_base = pd.merge(df_base, complex_traits_agg_genes, on='ratio_accession', how='left')
    else:
        df_base['mr_result_finngen'] = [[] for _ in range(len(df_base))]
        df_base['complex_traits'] = [[] for _ in range(len(df_base))]

    # --- Close DB Connection ---
    conn.close()
    if not silent:
        print("Database connection closed.")

    # --- Step 4: Final DataFrame cleaning for compatibility ---
    if not silent:
        print("Step 4: Finalizing the DataFrame...")
    # Fill NaN values in new list-like columns with empty lists
    list_cols = [
        'missense_genes', 'mr_result_eqtl', 'mr_result_pqtl', 'mr_result_gtex', 'mr_result_finngen',
        'eQTL_genes', 'pQTL_genes', 'GTEx_genes', 'complex_traits'
    ]
    for col in list_cols:
        df_base[col] = df_base[col].apply(lambda d: d if isinstance(d, list) else [])


    # Split 'closest_genes' and 'infomap_clusters' strings into lists of the correct type
    df_base['closest_genes'] = df_base['closest_genes'].apply(lambda x: x.split(',') if isinstance(x, str) else [])
    df_base['infomap_clusters'] = df_base['infomap_clusters'].fillna('')
    df_base['infomap_cluster'] = df_base['infomap_clusters'].apply(lambda x: int(x.split(',')[0]) if x else np.nan)

    # Replicate the metabolite name swap from the original notebook's logic
    df_base['met1_name'] = df_base['den_name']
    df_base['met2_name'] = df_base['num_name']
    df_base['met1_accession'] = df_base['den_accession']
    df_base['met2_accession'] = df_base['num_accession']

    # --- Assign to the final variable name used by the notebook ---
    all_rqtl_df = df_base.copy()

    # --- Display Final Results ---
    if not silent:
        print("\n--- Data Loading and Enrichment Complete ---")
        print(f"Final DataFrame has {all_rqtl_df.shape[0]} rows and {all_rqtl_df.shape[1]} columns.")

    key_cols = [
        'ratio_accession', 'pos_name', 'num_name', 'den_name',
        'eQTL_genes', 'pQTL_genes', 'missense_genes', 'closest_genes', 'complex_traits'
    ]
    if not silent:
        print("\nSample of the final enriched DataFrame:")
        print(all_rqtl_df[key_cols].drop_duplicates('ratio_accession').head())

    return all_rqtl_df