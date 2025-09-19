import pandas as pd

nonref_df = './genes_overlap_tes_nonref_no_duplicates.csv'
ref_df = './genes_overlap_tes_ref_no_duplicates.csv'

nonref_genes = pd.read_csv(nonref_df)
ref_genes = pd.read_csv(ref_df)

print("Number of non-reference genes:", nonref_genes.shape[0])
print("Number of reference genes:", ref_genes.shape[0])

# how many have something in the "gene_symbol" column?
print("Number of non-reference genes with gene_symbol:", nonref_genes['gene_symbol'].notnull().sum())
print("Number of reference genes with gene_symbol:", ref_genes['gene_symbol'].notnull().sum())