import sys
import pandas as pd
import numpy as np  
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib_venn import venn2

def updating_TE_names():
    dict_names = '/work/users/g/a/gabivla/lab/SV_mosquitos/TEs_Aedes_aegypti/running_pipeline_new_library/consensus_TEs_dict/TEforest_library/consensus_file_dictionary_TEforest.txt'
    name_dict = {}
    with open(dict_names, 'r') as f:
        for line in f:
            if line.startswith("Original:"):
                parts = line.strip().split("|")
                original_name = parts[0].split(":")[1].strip()
                updated_name  = parts[1].split(":")[1].strip()
                name_dict[updated_name] = original_name
    return name_dict

def combine_genotypes_vectorized(s1, s2):
    alt1 = s1 == '1/1'
    alt2 = s2 == '1/1'
    return np.where(alt1 | alt2, '1/1', '0/0')

def compare_retroseq_teforest(df_retroseq, df_teforest, tolerance):
    total_number_retroseq = len(df_retroseq)
    total_number_teforest = len(df_teforest)

    print(f'total number of rows in RetroSeq: {total_number_retroseq}')
    print(f'total number of rows in TEforest: {total_number_teforest}')

    base_cols_retroseq = ['chromosome', 'position_start', 'position_end', 'TE_name', 'strand', 'Unnamed: 0', 'Unnamed: 5']
    base_cols_teforest = ['chromosome', 'position_start', 'position_end', 'TE_name', 'strand', 'score']

    retro_geno_cols  = [c for c in df_retroseq.columns if c not in base_cols_retroseq]
    forest_geno_cols = [c for c in df_teforest.columns if c not in base_cols_teforest]

    # encode genotypes as int8 in one shot (no fragmentation yet)
    df_retroseq = df_retroseq.copy()
    df_teforest  = df_teforest.copy()
    df_retroseq[retro_geno_cols]  = df_retroseq[retro_geno_cols].replace({'0/0': 0, '1/1': 1}).astype(np.int8)
    df_teforest[forest_geno_cols] = df_teforest[forest_geno_cols].replace({'0/0': 0, '1/1': 1}).astype(np.int8)

    # drop these columns if they are there
    df_retroseq = df_retroseq.drop(columns=['Unnamed: 0', 'Unnamed: 5'], errors='ignore')
    df_teforest  = df_teforest.drop(columns=['score'], errors='ignore')

    base_cols = ['chromosome', 'position_start', 'position_end', 'TE_name', 'strand']

    retro_meta  = (df_retroseq[base_cols]
                   .rename(columns={'position_start': 'position_start_retro',
                                    'position_end':   'position_end_retro'})
                   .reset_index(names='retro_idx'))
    forest_meta = (df_teforest[base_cols]
                   .rename(columns={'position_start': 'position_start_forest',
                                    'position_end':   'position_end_forest'})
                   .reset_index(names='forest_idx'))

    # print the number of genotypes 
    print(f'number of genotype columns in RetroSeq: {len(retro_geno_cols)}')
    print(f'number of genotype columns in TEforest: {len(forest_geno_cols)}')

    df_merged = retro_meta.merge(forest_meta, on=['chromosome', 'TE_name', 'strand'], how='inner')

    position_match = (
        (np.abs(df_merged['position_start_retro']  - df_merged['position_start_forest']) <= tolerance) &
        (np.abs(df_merged['position_end_retro']    - df_merged['position_end_forest'])   <= tolerance)
    )
    matching_positions = df_merged.loc[position_match]

    print(f'candidate merged rows: {len(df_merged)}')
    print(f'matching rows after coordinate filter: {len(matching_positions)}')

    shared_samples = [c for c in retro_geno_cols  if c in forest_geno_cols]
    retro_only     = [c for c in retro_geno_cols  if c not in shared_samples]
    forest_only    = [c for c in forest_geno_cols if c not in shared_samples]

    retro_idx  = matching_positions['retro_idx'].values
    forest_idx = matching_positions['forest_idx'].values

    # ── build all columns at once to avoid DataFrame fragmentation ──────────
    meta_dict = {
        'chromosome':            matching_positions['chromosome'].values,
        # Keep both tool coordinates so we can compute final min/max below
        'position_start_retro':  matching_positions['position_start_retro'].values,
        'position_end_retro':    matching_positions['position_end_retro'].values,
        'position_start_forest': matching_positions['position_start_forest'].values,
        'position_end_forest':   matching_positions['position_end_forest'].values,
        'TE_name':               matching_positions['TE_name'].values,
        'strand':                matching_positions['strand'].values,
    }

    geno_dict = {}
    for sample in shared_samples:
        retro_vals  = df_retroseq.iloc[retro_idx][sample].to_numpy(dtype=np.int8)
        forest_vals = df_teforest.iloc[forest_idx][sample].to_numpy(dtype=np.int8)
        geno_dict[sample] = (retro_vals | forest_vals).astype(np.int8)
    for sample in retro_only:
        geno_dict[sample] = df_retroseq.iloc[retro_idx][sample].to_numpy(dtype=np.int8)
    for sample in forest_only:
        geno_dict[f'{sample}_forest'] = df_teforest.iloc[forest_idx][sample].to_numpy(dtype=np.int8)

    # single concat — no fragmentation
    result = pd.DataFrame({**meta_dict, **geno_dict})

    # ── aggregate duplicate coordinate pairs ────────────────────────────────
    group_cols   = ['chromosome', 'position_start_retro', 'position_end_retro',
                    'position_start_forest', 'position_end_forest', 'TE_name', 'strand']
    genotype_cols = [c for c in result.columns if c not in group_cols]

    agg_dict = {col: 'max' for col in genotype_cols}
    result = result.groupby(group_cols, as_index=False).agg(agg_dict)

    # ── Compute final positions: min of starts, max of ends (per diagram) ───
    result['position_start'] = result[['position_start_retro', 'position_start_forest']].min(axis=1)
    result['position_end']   = result[['position_end_retro',   'position_end_forest']].max(axis=1)

    # Drop intermediate coordinate columns — final positions are position_start / position_end
    result = result.drop(columns=['position_start_retro', 'position_end_retro',
                                  'position_start_forest', 'position_end_forest'])

    # back to VCF genotype strings
    all_sample_cols = shared_samples + retro_only + [f'{c}_forest' for c in forest_only]
    for col in all_sample_cols:
        if col in result.columns:
            result[col] = result[col].map({0: '0/0', 1: '1/1'})

    # ── Final column order ───────────────────────────────────────────────────
    ordered_cols = (
        ['chromosome', 'position_start', 'position_end', 'TE_name', 'strand']
        + shared_samples
        + retro_only
        + [f'{c}_forest' for c in forest_only]
    )
    result = result[ordered_cols]

    result.to_csv(f'matching_positions_retroseq_teforest_{tolerance}_{chromosome}.csv', index=False)

    pct_retro   = len(result) / total_number_retroseq * 100
    pct_forest  = len(result) / total_number_teforest * 100

    print(f'Matching positions: {len(result)}')
    print(f'Percentage of RetroSeq positions that match TEforest: {pct_retro:.2f}%')
    print(f'Percentage of TEforest positions that match RetroSeq: {pct_forest:.2f}%')
    print(f'Shared samples combined: {len(shared_samples)}')
    print(f'RetroSeq-only samples:   {len(retro_only)}')
    print(f'TEforest-only samples:   {len(forest_only)}')

    # ── Venn diagram ─────────────────────────────────────────────────────────
    retro_only_tes   = total_number_retroseq - len(result)
    forest_only_tes  = total_number_teforest - len(result)
    shared_tes       = len(result)

    if retro_only_tes == 0 and forest_only_tes == 0 and shared_tes == 0:
        print('Warning: no TEs to plot in Venn diagram — skipping.')
    else:
        fig, ax = plt.subplots(figsize=(8, 6))
        venn2(
            subsets=(retro_only_tes, forest_only_tes, shared_tes),
            set_labels=('RetroSeq', 'TEforest'),
            ax=ax,
        )
        ax.set_title('Overlap of TE insertions between RetroSeq and TEforest')
        fig.savefig(f'venn_diagram_TEs_{tolerance}_{chromosome}.png', bbox_inches='tight')
        plt.close(fig)

# ── Load ──────────────────────────────────────────────────────────────────────
tolerance = int(sys.argv[1])
chromosome = sys.argv[2]

dict_TE_names = updating_TE_names()
base_columns = ['chromosome', 'position_start', 'position_end', 'TE_name', 'strand']

# Load ALL columns (base + genotypes) from both files
df_retroseq = pd.read_csv(
    f"../../make_retroseq_initial_VCF/updated_post_processing/VCF_nonreference_post_processing_retroseq_{tolerance}.csv",
    dtype={
        'chromosome':     'category',
        'TE_name':        'category',
        'strand':         'category',
        'position_start': 'int32',
        'position_end':   'int32',
    }
)
df_retroseq = df_retroseq[df_retroseq['chromosome'] == chromosome].reset_index(drop=True)

df_teforest = pd.read_csv(
    f"../../make_teforest_initial_VCF/running_script_other_tolerance/merging_VCFs/TEforest_VCF_nonreference_TEs_{tolerance}.csv",
    dtype={
        'chromosome':     'category',
        'TE_name':        'category',
        'strand':         'category',
        'position_start': 'int32',
        'position_end':   'int32',
    }
)
df_teforest  = df_teforest[df_teforest['chromosome'] == chromosome].reset_index(drop=True)
df_teforest['TE_name'] = df_teforest['TE_name'].map(dict_TE_names).fillna(df_teforest['TE_name'])

compare_retroseq_teforest(df_retroseq, df_teforest, tolerance)