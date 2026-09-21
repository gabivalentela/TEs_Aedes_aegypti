import pandas as pd
import os

def remove_rows_without_TEs(df):
    mask = (df.iloc[:, 5:] == '1/1').any(axis=1)
    df = df[mask]
    return df

list_files = os.listdir('../')

list_dfs = []
for file in list_files:
    if file.endswith('.csv'):
        path_file = os.path.join('../', file)
        df = pd.read_csv(path_file)
        list_dfs.append(df)

df_concat = pd.concat(list_dfs, ignore_index=True)
df_nonref_updated = remove_rows_without_TEs(df_concat)
df_nonref_updated.to_csv('TEforest_VCF_nonreference_TEs.csv', index=False)
        