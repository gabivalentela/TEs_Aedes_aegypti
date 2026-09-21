import os
import pandas as pd

list_ref = '../ref'
list_nonref = '../nonref'

list_files_ref = [f for f in os.listdir(list_ref) if f.endswith('.csv')]
list_files_nonref = [f for f in os.listdir(list_nonref) if f.endswith('.csv')]

list_dfs_ref = {}
list_dfs_nonref = {}

for file in list_files_ref:
    df = pd.read_csv(os.path.join(list_ref, file))
    list_dfs_ref[file] = df

for file in list_files_nonref:
    df = pd.read_csv(os.path.join(list_nonref, file))
    list_dfs_nonref[file] = df

concat_df_ref = pd.concat(list_dfs_ref.values(), ignore_index=True)
concat_df_nonref = pd.concat(list_dfs_nonref.values(), ignore_index=True)

# update the dfs
concat_df_ref['type'] = 'ref'
concat_df_nonref['type'] = 'nonref'

# merge the two dataframes on the common columns adding ref and nonref as new column type
merged_df = pd.concat([concat_df_ref, concat_df_nonref], ignore_index=True)
merged_df.to_csv('all_high_frequency_results.csv', index=False)

# remove columns that start with frequency_*
merged_df = merged_df.loc[:, ~merged_df.columns.str.startswith('frequency_')]

# remove duplicate rows based on all columns
merged_df = merged_df.drop_duplicates()

# if gene_ID AAEL023266 make gene_name VGSC
merged_df.loc[merged_df['gene_ID'] == 'AAEL023266', 'gene_name'] = 'VGSC'

# get main family name - get all leters before the first number in the gene_name column
merged_df['main_family'] = merged_df['gene_name'].str.extract(r'([A-Za-z]+)')[0]

# count how many times each main_family appears in the merged_df
main_family_counts = merged_df['main_family'].value_counts()
print(main_family_counts)

# Report the families with the most TE occurrences. Each row represents one TE copy.
family_summary = (
    merged_df.dropna(subset=['main_family'])
    .groupby('main_family')
    .agg(
        unique_genes=('gene_name', 'nunique'),
        te_copies=('TE_name', 'size'),
    )
    .sort_values('te_copies', ascending=False)
    .head(4)
)

print('\nTop 4 gene families by TE count:')
for main_family, summary in family_summary.iterrows():
    print(
        f"{main_family} has {summary['unique_genes']} genes "
        f"and {summary['te_copies']} TEs/copies"
    )
