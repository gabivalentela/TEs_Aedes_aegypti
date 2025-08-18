import os
import polars as pl

general_path = '/work/users/g/a/gabivla/lab/SV_mosquitos/VCF_genotypes_TEs/separate_chromosomes_windows_reference/reference_final/AaegL5_3'
list_files = os.listdir(general_path)

dfs = []
all_columns = set()

for file in list_files:
    df = pl.read_csv(os.path.join(general_path, file))
    df = df.with_columns([pl.col(col).cast(pl.Utf8) for col in df.columns])
    if "" in df.columns:
        df = df.drop("")
    dfs.append(df)
    all_columns.update(df.columns)

count_rows = sum(df.shape[0] for df in dfs)
print(f'The number of rows the final df has to have: {count_rows}')

priority_cols = ['chromosome', 'position_start', 'position_end', 'TE_name', 'strand']
other_cols = sorted([col for col in all_columns if col not in priority_cols])
final_columns = priority_cols + other_cols

aligned_dfs = []
for df in dfs:
    missing_cols = [col for col in final_columns if col not in df.columns]
    for col in missing_cols:
        df = df.with_columns(pl.lit("0/0").alias(col))
    df = df.select(final_columns)
    aligned_dfs.append(df)

final_df = pl.concat(aligned_dfs, how="vertical")
print(f"✅ Merged DataFrame has {final_df.shape[0]} rows and {final_df.shape[1]} columns.")

final_df.write_csv(
    '/work/users/g/a/gabivla/lab/SV_mosquitos/VCF_genotypes_TEs/separate_chromosomes_windows_reference/merged_reference/separate_chromosomes/VCF_reference_TEs_AaegL5_3_pre_processing.csv'
)
