import pandas as pd
import os
path_becca = '/nas/longleaf/home/gabivla/lab/SV_mosquitos/IR_genes/becca_ir_list_coords.txt'
path_list = '/nas/longleaf/home/gabivla/lab/SV_mosquitos/IR_genes/cyp_gst_list_coords.txt'

df_becca = pd.read_csv(path_becca, sep="\t")
df_list = pd.read_csv(path_list, sep="\t", header=None)

df_becca = df_becca[['chrom', 'start', 'end', 'other name', '##gene ID']]

df_list[4] = df_list[3].str.extract(r'.+;Name=(.+);')[0]
df_list[5] = df_list[3].str.extract(r'ID=(.+);')[0]

df_list = df_list.drop(columns=[3])

df_list = df_list.rename(columns={0: 'chrom', 1: 'start', 2: 'end', 4: 'other name', 5: '##gene ID'})

df_merged = pd.concat([df_becca, df_list], ignore_index=True)

# use gene ID to get unique genes
unique_genes = df_merged['##gene ID'].dropna().unique()

ref_files = []
# get files ref
path_ref = './ref'
list_ref = os.listdir(path_ref)
for ref in list_ref:
    if ref.endswith('.csv'):
        df = pd.read_csv(os.path.join(path_ref, ref))
        ref_files.append(df)
        
concat_ref = pd.concat(ref_files, ignore_index=True)
# get unique gene_ID from concat_ref
unique_genes_ref = concat_ref['gene_ID'].dropna().unique()

# get intersection of unique_genes and unique_genes_ref
intersection = set(unique_genes).intersection(set(unique_genes_ref))
print(f"Number of unique genes in the list: {len(unique_genes)}")
print(f"Number of unique genes in reference-based analysis: {len(unique_genes_ref)}")
print(f"Number of unique genes in both the list and reference-based analysis: {len(intersection)}")
# get files nonref
nonref_files = []
path_nonref = './nonref'
list_nonref = os.listdir(path_nonref)
for nonref in list_nonref:
    if nonref.endswith('.csv'):
        df = pd.read_csv(os.path.join(path_nonref, nonref))
        nonref_files.append(df) 
concat_nonref = pd.concat(nonref_files, ignore_index=True)
# get unique gene_ID from concat_nonref
unique_genes_nonref = concat_nonref['gene_ID'].dropna().unique()
# get intersection of unique_genes and unique_genes_nonref
intersection_nonref = set(unique_genes).intersection(set(unique_genes_nonref))
print(f"Number of unique genes in non-reference-based analysis: {len(unique_genes_nonref)}")
print(f"Number of unique genes in both the list and non-reference-based analysis: {len(intersection_nonref)}")
# get intersection of unique_genes_ref and unique_genes_nonref
intersection_ref_nonref = set(unique_genes_ref).intersection(set(unique_genes_nonref))
print(f"Number of unique genes in both reference-based and non-reference-based analysis: {len(intersection_ref_nonref)}")   
