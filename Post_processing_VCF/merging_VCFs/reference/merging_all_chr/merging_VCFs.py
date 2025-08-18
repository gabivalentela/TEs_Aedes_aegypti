import os
import polars as pl

general_path = '/work/users/g/a/gabivla/lab/SV_mosquitos/VCF_genotypes_TEs/separate_chromosomes_windows_reference/merged_reference/separate_chromosomes/'
list_files = os.listdir(general_path)
dfs = []
all_columns = set()
for file in list_files:
    df = pl.read_csv(os.path.join(general_path, file))
    df = df.with_columns([pl.col(col).cast(pl.Utf8) for col in df.columns])  # Convert all columns to String
    df = df.drop("")
    dfs.append(df)
    all_columns.update(df.columns)

count_rows = 0
for df in dfs:
    rows = df.shape[0]
    count_rows += rows

print(f'the number of rows the final df has to have: {count_rows}')

priority_cols = ['chromosome', 'position_start', 'position_end', 'TE_name', 'strand']
other_cols = sorted([col for col in all_columns if col not in priority_cols])
all_columns = priority_cols + other_cols

aligned_dfs = []
for df in dfs:
    missing_cols = [col for col in all_columns if col not in df.columns]
    for col in missing_cols:
        df = df.with_columns(pl.lit("0/0").alias(col))
    df = df.select(all_columns)
    aligned_dfs.append(df)

merged_df = pl.concat(aligned_dfs, how="vertical")
print(merged_df)
print(f"✅ Merged DataFrame has {merged_df.shape[0]} rows and {merged_df.shape[1]} columns.")

merged_df.write_csv("/work/users/g/a/gabivla/lab/SV_mosquitos/VCF_genotypes_TEs/separate_chromosomes_windows_reference/merged_reference/VCF_reference_TEs_pre_processing.csv")