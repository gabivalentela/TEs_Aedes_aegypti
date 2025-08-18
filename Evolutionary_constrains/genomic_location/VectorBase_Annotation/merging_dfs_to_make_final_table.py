import pandas as pd

nonref = 'genes_overlap_tes_nonref_no_duplicates.csv'
ref = 'genes_overlap_tes_ref_no_duplicates.csv'

df_ref = pd.read_csv(ref)
df_nonref = pd.read_csv(nonref)

# column names
nonref_columns = df_nonref.columns.tolist()
ref_columns = df_ref.columns.tolist()

print("Original nonref columns:", nonref_columns)
print("Original ref columns:", ref_columns)

# Remove 'Unnamed: 0' columns if they exist
if 'Unnamed: 0' in df_ref.columns:
    df_ref = df_ref.drop('Unnamed: 0', axis=1)
if 'Unnamed: 0' in df_nonref.columns:
    df_nonref = df_nonref.drop('Unnamed: 0', axis=1)

# Function to rename columns that start with 'count_1/1_'
def rename_count_columns(df):
    new_columns = {}
    for col in df.columns:
        if col.startswith('count_1/1_'):
            # Extract country name after 'count_1/1_'
            country = col.replace('count_1/1_', '')
            new_columns[col] = f'# genomes with TE in {country}'
    return df.rename(columns=new_columns)

# Rename columns in both dataframes
df_ref = rename_count_columns(df_ref)
df_nonref = rename_count_columns(df_nonref)

print("\nColumns after renaming:")
print("Ref columns:", df_ref.columns.tolist())
print("Nonref columns:", df_nonref.columns.tolist())

# Add source column to identify reference vs nonreference
df_ref['source'] = 'reference'
df_nonref['source'] = 'nonreference'

# Merge the dataframes
merged_df = pd.concat([df_ref, df_nonref], ignore_index=True)

print(f"\nMerged dataframe shape: {merged_df.shape}")
print("Final columns:", merged_df.columns.tolist())
print("\nSource value counts:")
print(merged_df['source'].value_counts())

# Save the merged dataframe
merged_df.to_csv('merged_genes_overlap_tes.csv', index=False)
print("\nMerged dataframe saved as 'merged_genes_overlap_tes.csv'")

# Display first few rows to verify
print("\nFirst 5 rows of merged dataframe:")
print(merged_df.head())