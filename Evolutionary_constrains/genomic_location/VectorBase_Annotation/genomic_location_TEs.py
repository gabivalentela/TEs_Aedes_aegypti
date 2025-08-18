import pandas as pd
import re

# Function to load TE data
def make_df(path):
    df = pd.read_csv(path)

    # Split TE_information using regex
    pattern = r'(.+_\d)_(\d+)_(\d+)_(.+)' 

    df[['chromosome', 'start', 'end', 'te_name']] = df['TE_information'].str.extract(pattern)

    # Convert start and end to integers
    df[['start', 'end']] = df[['start', 'end']].astype(int)
    return df

# Function to load gene annotation data and modify the Chromosome column
def gene_annotation_info(file_path):
    df = pd.read_csv(file_path, delimiter='\t')
    # Update the content in the 'Chromosome' column to match TE format
    df['Chromosome'] = df['Chromosome'].apply(lambda x: f'AaegL5_{x}')
    return df

def find_overlapping_genes(te_df, gene_annotation_info_df):
    # Create empty columns for gene information
    te_df['gene_symbol'] = None
    te_df['gene_id'] = None
    te_df['gene_name'] = None 

    for index, row in te_df.iterrows():
        te_chromosome = row['chromosome']
        te_start_position = int(row['start'])
        te_end_position = int(row['end'])

        # Initialize variables to store gene info for each TE
        matching_gene_symbol = None
        matching_gene_id = None
        matching_gene_name = None  

        for _, row_df in gene_annotation_info_df.iterrows():
            gene_chromosome = row_df['seqid']
            gene_start_position = int(row_df['start'])
            gene_end_position = int(row_df['end'])
            gene_symbol = row_df['ID']
            gene_id = row_df['type']
            gene_name = row_df['attributes']

            interval_start = gene_start_position - 10000
            interval_end = gene_end_position + 10000

            # Check if TE and gene overlap or are nearby
            if te_chromosome == gene_chromosome and interval_start < te_start_position and te_end_position < interval_end:
                matching_gene_symbol = gene_symbol
                matching_gene_id = gene_id
                matching_gene_name = gene_name 
                break  

        # Update the te_df DataFrame
        te_df.at[index, 'gene_symbol'] = matching_gene_symbol
        te_df.at[index, 'gene_id'] = matching_gene_id
        te_df.at[index, 'gene_name'] = matching_gene_name  

    return te_df

def read_data_annotation(path_annotation):
    data = pd.read_csv(
        path_annotation,
        sep='\t',
        comment='#',  
        header=None, 
    )
    
    data.columns = [
        'seqid', 'source', 'type', 'start', 'end', 'score',
        'strand', 'phase', 'attributes'
    ]
    return data


# Load data
ref_TEs = make_df('../../reference/no_duplicates/final_df_reference.csv')
nonref_TEs = make_df('../../nonreference/no_duplicates/final_df_nonreference.csv')

path_annotation = '/nas/longleaf/home/gabivla/lab/SV_mosquitos/final_scripts/Annotation_VectorBase/VectorBase-54_AaegyptiLVP_AGWG.gff'
df_annotation = read_data_annotation(path_annotation)
# split the 'attributes' column into multiple columns
df_annotation['ID'] = df_annotation['attributes'].str.split(';').str[0]

# Find overlapping genes for reference TEs
ref_TEs_with_genes = find_overlapping_genes(ref_TEs, df_annotation)
ref_TEs_with_genes.to_csv('genes_overlap_tes_ref_no_duplicates.csv')

# Find overlapping genes for non-reference TEs
nonref_TEs_with_genes = find_overlapping_genes(nonref_TEs, df_annotation)
print(nonref_TEs_with_genes)
nonref_TEs_with_genes.to_csv('genes_overlap_tes_nonref_no_duplicates_new.csv')