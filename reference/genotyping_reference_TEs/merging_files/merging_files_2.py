import pandas as pd
import os
import re
from functools import reduce
from collections import defaultdict

list_files = os.listdir('../')
base_columns = ['chromosome', 'position_start', 'position_end', 'strand', 'TE_name', 'score']

# adjust this pattern to match your actual filenames
file_pattern = r'df_merged_onecode_strict_subset_(AaegL5_\d+)_(.+)\.csv'

for file in list_files:
    if file.startswith('results_') and file.endswith('2_either'):
        pattern = r'results_(\d+)_(.+)'
        match = re.match(pattern, file)
        if match:
            number_supporting_reads = match.group(1)
            side_condition = match.group(2)
            print(f'Number of supporting reads: {number_supporting_reads}, Side condition: {side_condition}')
            list_files_results = os.listdir(f'../{file}')

            # group files by chromosome
            chrom_groups = defaultdict(list)
            for results in list_files_results:
                m = re.match(file_pattern, results)
                if not m:
                    print(f'  Skipping unrecognized file: {results}')
                    continue
                chromosome = m.group(1)
                sample_name = m.group(2)
                chrom_groups[chromosome].append((sample_name, results))

            chrom_dfs = []
            for chromosome, samples in chrom_groups.items():
                sample_dfs = []
                for sample_name, results in samples:
                    df = pd.read_csv(f'../{file}/{results}')
                    genotype_cols = [c for c in df.columns if c not in base_columns]
                    df = df.rename(columns={c: f'{c}' for c in genotype_cols})

                    sample_dfs.append(df)

                # merge samples belonging to the SAME chromosome, on base columns
                merged_chrom_df = reduce(
                    lambda left, right: pd.merge(left, right, on=base_columns, how='outer'),
                    sample_dfs
                )
                chrom_dfs.append(merged_chrom_df)

            # stack all chromosomes together (rows, not columns)
            final_df = pd.concat(chrom_dfs, axis=0, ignore_index=True)

            final_df.to_csv(f'merged_results_{number_supporting_reads}_{side_condition}.csv', index=False)