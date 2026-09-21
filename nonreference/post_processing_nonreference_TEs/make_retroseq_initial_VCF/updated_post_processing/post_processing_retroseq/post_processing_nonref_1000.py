import pandas as pd
import polars as pl
import sys

def remove_rows_without_TEs(df):
    mask = (df.iloc[:, 5:] == '1/1').any(axis=1)
    df = df[mask]
    return df
   
non_reference_VCF = '../nonreference_final/VCF_nonreference_TEs_pre_processing_1000.csv'
df_nonref = pd.read_csv(non_reference_VCF)
print(f'Initial number of rows: {df_nonref.shape[0]}')
df_nonref = df_nonref.sort_values(by=['chromosome', 'position_start'])
df_nonref_updated = remove_rows_without_TEs(df_nonref)
print(f'Number of rows after removing rows without TEs: {df_nonref_updated.shape[0]}')
# check for duplicated rows in this df
duplicated_rows = df_nonref_updated[df_nonref_updated.duplicated()]
print(duplicated_rows)
df_nonref_updated.to_csv('../VCF_nonreference_post_processing_retroseq_1000.csv')
