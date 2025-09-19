import pandas as pd

def print_gene_rows(annotation_vector_base, gene_code):
    # Open GFF file as a pandas DataFrame
    gff_df = pd.read_csv(annotation_vector_base, sep='\t', comment='#', header=None)

    # Filter for rows where the last column (column 8) contains the gene code
    gene_rows = gff_df[gff_df[8].str.contains(gene_code, na=False)]

    # Print the matching rows
    print(gene_rows)

# Example usage
annotation_vector_base = '../../../../../Annotation_VectorBase/VectorBase-54_AaegyptiLVP_AGWG.gff'
print_gene_rows(annotation_vector_base, 'AAEL011934')
