import pandas as pd
import os
import sys
import re

variation = sys.argv[1]  
list_files = os.listdir('../')

list_dfs = []
for file in list_files:
    if file.endswith('.csv'):
        variation_file = re.search('matching_positions_retroseq_teforest_(.+)_A.+.csv', file).group(1)
        if variation_file == variation:
            updated_path = f'../{file}'
            df = pd.read_csv(updated_path)
            list_dfs.append(df)

# Concatenate all DataFrames in the list
merged_df = pd.concat(list_dfs, ignore_index=True)
merged_df.to_csv(f'matching_positions_retroseq_teforest_{variation}_all_chromosomes.csv', index=False)