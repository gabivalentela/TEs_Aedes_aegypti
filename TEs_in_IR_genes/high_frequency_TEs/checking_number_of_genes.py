import os 
import pandas as pd
import re

list_ref = os.listdir("./ref")
list_nonref = os.listdir("./nonref")

list_dfs_ref = []
for ref in list_ref:
    if '.csv' in ref:
        country = re.search('TEs_in_IR_genes_(.+).csv', ref).group(1)
        path = os.path.join("./ref", ref)
        df = pd.read_csv(path)
        df['country'] = country
        # change frequency_* column to just frequency
        df.columns = df.columns.str.replace(r'^frequency_.*', 'frequency', regex=True)
        # merge dataframes
        list_dfs_ref.append(df)

general_df = pd.concat(list_dfs_ref, ignore_index=True)
# get unique gene_name counts
unique_genes = general_df['gene_name'].dropna().unique()

list_cyp_genes = []
# get CYP gene
for unique_gene in unique_genes:
    if 'nan' in unique_gene:
        continue
    elif 'NVY' in unique_gene:
        list_cyp_genes.append(unique_gene)
        print(unique_gene)

#print(f"Number of unique genes in reference-based analysis: {len(unique_genes)}")
print(f"Number of unique NVY genes in reference-based analysis: {len(list_cyp_genes)}")

# get only rows from general_df that had these genes from the list_cyp_genes
cyp_df = general_df[general_df['gene_name'].isin(list_cyp_genes)]
# get unique TE names
unique_te_names = cyp_df['TE_name'].dropna().unique()
# number of all TEs
print(f"Number of all TEs in reference-based analysis: {cyp_df.shape[0]}")
# number of unique TEs names
print(f"Number of unique TEs in NVY insertions: {len(unique_te_names)}")