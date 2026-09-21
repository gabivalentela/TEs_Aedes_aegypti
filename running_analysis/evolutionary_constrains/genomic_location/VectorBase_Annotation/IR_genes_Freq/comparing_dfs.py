import pandas as pd
import sys

def interpreting_results(path, type_of_data):
    # Load and clean nonref data
    df = pd.read_csv(path)
    df['gene_symbol'] = df['gene_symbol'].str.replace('ID=', '', regex=False)
    df = df[['gene_symbol', 'TE_information']]

    # Load IR gene lists
    path_becca = '/nas/longleaf/home/gabivla/lab/SV_mosquitos/IR_genes/becca_ir_list_coords.txt'
    path_list = '/nas/longleaf/home/gabivla/lab/SV_mosquitos/IR_genes/cyp_gst_list_coords.txt'

    df_becca = pd.read_csv(path_becca, sep='\t')
    df_list = pd.read_csv(path_list, sep='\t', header=None)

    # Clean df_list
    df_list.columns = ['chromosome', 'start', 'end', 'TE_information']
    df_list['AAEL_code'] = df_list['TE_information'].str.extract(r'ID=(AAEL\d+);')

    # Filter rows where gene is found in nonref, then merge to get TE info
    df_becca_matches = df_becca[df_becca['##gene ID'].isin(df['gene_symbol'])].merge(
        df,
        left_on='##gene ID',
        right_on='gene_symbol',
        how='left'
    )

    df_list_matches = df_list[df_list['AAEL_code'].isin(df['gene_symbol'])].merge(
        df,
        left_on='AAEL_code',
        right_on='gene_symbol',
        how='left'
    )

    # Save results
    df_becca_matches.to_csv(f'TEs_becca_ir_list_aael_{type_of_data}.csv', index=False)
    df_list_matches.to_csv(f'TEs_cyp_gst_list_aael_{type_of_data}.csv', index=False)

    # Print results
    print(f"=== Becca IR list matches ({len(df_becca_matches)} genes) ===")
    print(df_becca_matches)

    print(f"\n=== CYP/GST list matches ({len(df_list_matches)} genes) ===")
    print(df_list_matches)

path_ref = sys.argv[1]
path_nonref = sys.argv[2]

interpreting_results(path_ref, 'ref')
interpreting_results(path_nonref, 'nonref')

