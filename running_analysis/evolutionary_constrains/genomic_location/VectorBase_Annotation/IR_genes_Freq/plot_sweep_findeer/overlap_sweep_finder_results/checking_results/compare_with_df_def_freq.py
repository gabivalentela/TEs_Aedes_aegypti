import pandas as pd

ref_df_sweeps = pd.read_csv('./reference_top1percent_outliers.csv')
nonref_df_sweeps = pd.read_csv('./nonreference_top1percent_outliers.csv')

# dfs ref and nonref - biggest difference
ref_df_def_freq = pd.read_csv('../../TEs_Aedes_aegypti/Evolutionary_constrains/genomic_location/VectorBase_Annotation/genes_overlap_tes_ref_no_duplicates.csv')
nonref_df_def_freq = pd.read_csv('../../TEs_Aedes_aegypti/Evolutionary_constrains/genomic_location/VectorBase_Annotation/genes_overlap_tes_nonref_no_duplicates.csv')

def split_te_info(df):
    pattern = r'.*(AaegL5_\d)_(\d+)_(\d+)_(.+)_(.+)'

    extracted = df['TE_information'].str.extract(pattern)
    
    df['chromosome'] = extracted[0]
    df['position_start'] = extracted[1].astype(int)
    df['position_end'] = extracted[2].astype(int)
    df['TE_name'] = extracted[3]
    df['strand'] = extracted[4]
    
    return df

ref_df_def_freq = split_te_info(ref_df_def_freq)
nonref_df_def_freq = split_te_info(nonref_df_def_freq)

# Keep all TEs from ref_df_def_freq and add sweep information where available
merged_ref = pd.merge(ref_df_sweeps, ref_df_def_freq, 
                      on=['chromosome', 'position_start', 'position_end', 'TE_name', 'strand'], 
                      how='right')

merged_nonref = pd.merge(nonref_df_sweeps, nonref_df_def_freq, 
                         on=['chromosome', 'position_start', 'position_end', 'TE_name', 'strand'], 
                         how='right')

print("Reference merged data:")
print(merged_ref)
print(f"\nTotal TEs from ref_df_def_freq: {len(ref_df_def_freq)}")
print(f"TEs with sweep data: {merged_ref['sweep_LR'].notna().sum()}")
print(f"TEs without sweep data: {merged_ref['sweep_LR'].isna().sum()}")

print("\n" + "="*60 + "\n")

print("Non-reference merged data:")
print(merged_nonref)
print(f"\nTotal TEs from nonref_df_def_freq: {len(nonref_df_def_freq)}")
print(f"TEs with sweep data: {merged_nonref['sweep_LR'].notna().sum()}")
print(f"TEs without sweep data: {merged_nonref['sweep_LR'].isna().sum()}")

# Save merged data to CSV files
print("\n" + "="*60)
print("Saving merged data to CSV files...")
print("="*60)

# Reference data - all merged
merged_ref.to_csv('reference_merged_all.csv', index=False)
print("✓ Saved: reference_merged_all.csv")

# Reference data - with sweep data only
ref_with_sweeps = merged_ref[merged_ref['sweep_LR'].notna()]
ref_with_sweeps.to_csv('reference_with_sweeps.csv', index=False)
print(f"✓ Saved: reference_with_sweeps.csv ({len(ref_with_sweeps)} TEs)")

# Reference data - without sweep data only
ref_without_sweeps = merged_ref[merged_ref['sweep_LR'].isna()]
ref_without_sweeps.to_csv('reference_without_sweeps.csv', index=False)
print(f"✓ Saved: reference_without_sweeps.csv ({len(ref_without_sweeps)} TEs)")

print()

# Non-reference data - all merged
merged_nonref.to_csv('nonreference_merged_all.csv', index=False)
print("✓ Saved: nonreference_merged_all.csv")

# Non-reference data - with sweep data only
nonref_with_sweeps = merged_nonref[merged_nonref['sweep_LR'].notna()]
nonref_with_sweeps.to_csv('nonreference_with_sweeps.csv', index=False)
print(f"✓ Saved: nonreference_with_sweeps.csv ({len(nonref_with_sweeps)} TEs)")

# Non-reference data - without sweep data only
nonref_without_sweeps = merged_nonref[merged_nonref['sweep_LR'].isna()]
nonref_without_sweeps.to_csv('nonreference_without_sweeps.csv', index=False)
print(f"✓ Saved: nonreference_without_sweeps.csv ({len(nonref_without_sweeps)} TEs)")

print("\n" + "="*60)
print("Summary of saved files:")
print("="*60)
print(f"Reference TEs:")
print(f"  - Total: {len(merged_ref)}")
print(f"  - With sweeps: {len(ref_with_sweeps)}")
print(f"  - Without sweeps: {len(ref_without_sweeps)}")
print(f"\nNon-reference TEs:")
print(f"  - Total: {len(merged_nonref)}")
print(f"  - With sweeps: {len(nonref_with_sweeps)}")
print(f"  - Without sweeps: {len(nonref_without_sweeps)}")