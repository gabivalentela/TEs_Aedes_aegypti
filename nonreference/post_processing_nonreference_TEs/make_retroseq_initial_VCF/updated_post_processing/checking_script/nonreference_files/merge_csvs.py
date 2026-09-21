import pandas as pd
import os

list_files = os.listdir('./')
list_dfs = []
for file in list_files:
    if file.endswith('.csv'):
        variation = file.split('.')[-2]
        variation = variation.split('_')[-1]
        df = pd.read_csv(file)
        # if TE_name == TE_1463_SINE_SINE_SINE_SINE
        #df['TE_name'] = df[df['TE_name']== 'TE_1463_SINE_SINE_SINE_SINE']
        df = df.loc[df['TE_name'] == 'TE_1463_SINE_SINE_SINE_SINE']
        first_row = df.head(1)
        first_row['variation'] = variation
        list_dfs.append(first_row)
df_merged = pd.concat(list_dfs, ignore_index=True)
df_merged.to_csv('merged_TE_1463_SINE_SINE_SINE_SINE.csv', index=False)