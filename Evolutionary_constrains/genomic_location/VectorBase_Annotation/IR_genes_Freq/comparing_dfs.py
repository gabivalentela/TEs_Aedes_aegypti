import pandas as pd

path_nonref = '../genes_overlap_tes_nonref_no_duplicates.csv'
path_ref = '../genes_overlap_tes_ref_no_duplicates.csv'

df_nonref = pd.read_csv(path_nonref)
df_ref = pd.read_csv(path_ref)

# remove ID= from gene_symbol
df_nonref['gene_symbol'] = df_nonref['gene_symbol'].str.replace('ID=', '', regex=False)
df_ref['gene_symbol'] = df_ref['gene_symbol'].str.replace('ID=', '', regex=False)

# get AAEL ids using LOC info
## get only columns gene_symbol and TE_information
df_nonref = df_nonref[['gene_symbol', 'TE_information']]
df_ref = df_ref[['gene_symbol', 'TE_information']]

# df IR genes
path_becca = '/nas/longleaf/home/gabivla/lab/SV_mosquitos/IR_genes/becca_ir_list_coords.txt'
path_list = '/nas/longleaf/home/gabivla/lab/SV_mosquitos/IR_genes/cyp_gst_list_coords.txt'

df_becca = pd.read_csv(path_becca, sep='\t')
df_list = pd.read_csv(path_list, sep='\t', header=None)

# name columns df_list
df_list.columns = ['chromosome', 'start', 'end', 'TE_information']
df_list['AAEL_code'] = df_list.iloc[:, 3].str.extract(r'ID=(AAEL\d+);')

# get Loc_code from df_b and df_l - these are list get only first item
df_b = df_becca[['##gene ID']]
df_l = df_list[['AAEL_code']]

# compare two dfs
## compare df_b and df_nonref and df_ref
df_b['in_nonref'] = df_b['##gene ID'].isin(df_nonref['gene_symbol'])
df_b['in_ref'] = df_b['##gene ID'].isin(df_ref['gene_symbol'])
df_l['in_nonref'] = df_l['AAEL_code'].isin(df_nonref['gene_symbol'])
df_l['in_ref'] = df_l['AAEL_code'].isin(df_ref['gene_symbol'])

df_l.to_csv('TEs_cyp_gst_list_aael.csv')
df_b.to_csv('TEs_becca_ir_list_aael.csv')

# print only columns with True
print(df_b[df_b['in_nonref']])
print(df_b[df_b['in_ref']])
print(df_l[df_l['in_nonref']])
print(df_l[df_l['in_ref']])


