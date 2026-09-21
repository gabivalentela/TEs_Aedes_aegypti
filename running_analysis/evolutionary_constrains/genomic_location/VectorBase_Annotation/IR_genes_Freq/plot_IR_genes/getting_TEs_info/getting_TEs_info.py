import pandas as pd
import os

general_list = os.listdir('../../')

ref_files = []
nonref_files = []

for file in general_list:
    if file.endswith('_ref.csv'):
        path_reference = '../../' + file
        ref_files.append(path_reference)
    elif file.endswith('_nonref.csv'):
        path_nonreference = '../../' + file
        nonref_files.append(path_nonreference)
        
list_dfs_ref = []
for ref in ref_files:
    df = pd.read_csv(ref)
    list_dfs_ref.append(df)

list_dfs_nonref = []
for nonref in nonref_files:
    df = pd.read_csv(nonref)
    list_dfs_nonref.append(df)

# concat all the dataframes in the list
final_df_ref = pd.concat(list_dfs_ref, ignore_index=True)
final_df_nonref = pd.concat(list_dfs_nonref, ignore_index=True)
# keep only TE_information column
final_df_ref = final_df_ref[['TE_information']]
final_df_nonref = final_df_nonref[['TE_information']]

path_reference = '../../../genes_overlap_tes_ref_no_duplicates.csv'
path_nonreference = '../../../genes_overlap_tes_nonref_no_duplicates_1000.csv'
df_nonreference = pd.read_csv(path_nonreference)
df_reference = pd.read_csv(path_reference)

def get_TE_information(df, TEs_df, type_of_df):
    # get the rows from path_nonreference that have the same TE_information
    subset_df = df[df['TE_information'].isin(TEs_df['TE_information'])]
    # remove first column
    subset_df = subset_df.iloc[:, 1:]
    # save the dataframe to a csv file
    subset_df.to_csv(f'TE_information_all_genes_overlap_{type_of_df}.csv', index=False)  
    
get_TE_information(df_reference, final_df_ref, 'reference')
get_TE_information(df_nonreference, final_df_nonref, 'nonreference')