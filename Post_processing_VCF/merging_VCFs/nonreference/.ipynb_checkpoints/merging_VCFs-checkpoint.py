import os
import polars as pl

general_path = '/work/users/g/a/gabivla/lab/SV_mosquitos/VCF_genotypes_TEs/separate_chromosomes_windows_reference/nonreference_files/'
list_files = os.listdir(general_path)

dfs = []
all_columns = set()

# First pass: read files and collect all column names
for file in list_files:
    df = pl.read_csv(os.path.join(general_path, file))
    df = df.with_columns([pl.col(col).cast(pl.Utf8) for col in df.columns])  # Convert all columns to string
    dfs.append(df)
    all_columns.update(df.columns)

# Prioritize metadata columns
priority_cols = ['chromosome', 'position_start', 'position_end', 'TE_name', 'strand']
other_cols = sorted([col for col in all_columns if col not in priority_cols])
final_column_order = priority_cols + other_cols

# Align all DataFrames to have the same columns
aligned_dfs = []
for df in dfs:
    missing_cols = [col for col in final_column_order if col not in df.columns]
    for col in missing_cols:
        df = df.with_columns(pl.lit("0/0").alias(col))
    # Reorder columns
    df = df.select(final_column_order)
    aligned_dfs.append(df)

# Merge (concatenate) all DataFrames vertically
final_df = pl.concat(aligned_dfs, how="vertical")

# Print total rows
print(f'The number of rows the final df has to have: {final_df.shape[0]}')

# Save to CSV
final_df.write_csv("/work/users/g/a/gabivla/lab/SV_mosquitos/VCF_genotypes_TEs/separate_chromosomes_windows_reference/merged_nonreference/VCF_nonreference_TEs_pre_processing.csv")
