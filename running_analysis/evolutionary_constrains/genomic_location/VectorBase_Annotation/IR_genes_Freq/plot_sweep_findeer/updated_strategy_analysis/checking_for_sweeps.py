import pandas as pd
import os

def overlapping_sweeps_TEs(df_TEs_genes, df_sweep):
    # check for overlap between sweeps and TEs with overlap status
    overlapping_TEs = []
    for index, sweep_row in df_sweep.iterrows():
        sweep_position = sweep_row['location']
        sweep_chrom = sweep_row['chrom']
        sweep_LR = sweep_row['LR']  # Capture LR value
        sweep_pop = sweep_row['pop']  # Capture pop value
        
        overlapping = df_TEs_genes[
            (df_TEs_genes['chromosome'] == sweep_chrom) &
            (df_TEs_genes['start'] - 10000 <= sweep_position) &
            (df_TEs_genes['end'] + 10000 >= sweep_position)
        ]
        
        if not overlapping.empty:
            # Add sweep information to each overlapping TE
            overlapping_copy = overlapping.copy()
            overlapping_copy['sweep_location'] = sweep_position
            overlapping_copy['sweep_chrom'] = sweep_chrom
            overlapping_copy['sweep_LR'] = sweep_LR
            overlapping_copy['sweep_pop'] = sweep_pop
            overlapping_TEs.append(overlapping_copy)
    
    # Combine all overlapping TEs into a single dataframe
    if overlapping_TEs:
        return pd.concat(overlapping_TEs, ignore_index=True)
    else:
        return pd.DataFrame()  # Return empty dataframe if no overlaps
    
def read_data_annotation(path_annotation):
    data = pd.read_csv(
        path_annotation,
        sep='\t',
        comment='#',  # Skip lines starting with '#'
        header=None,  # GFF files don't usually have column headers
    )
    # Assign column names based on GFF3 specification
    data.columns = [
        'seqid', 'source', 'type', 'start', 'end', 'score',
        'strand', 'phase', 'attributes'
    ]
    return data

path_annotation = '../../../../../../../../../SV_mosquitos/TEs_Aedes_aegypti/Annotation_VectorBase/VectorBase-54_AaegyptiLVP_AGWG.gff'
df_annotation = read_data_annotation(path_annotation)
exons = df_annotation[df_annotation['type'] == 'exon']

# get the gene ID from the genes
exons['gene_id'] = exons['attributes'].str.extract(r'gene_id=([^;]+)')

df_TEs_genes = pd.read_csv('./genes_overlap_tes_ref_with_GO_noduplicate.csv')
# extract from gene name
df_TEs_genes['gene_info'] = df_TEs_genes['gene_name'].str.extract(r'ID=([^;]+)')

# open new dfs for sweeps
general_paths_sweeps = './make_new_sweep_files/new_sweeps_files/'
list_sweep_files = os.listdir(general_paths_sweeps)

all_results = []  # Store results from all sweep files

for sweep_file in list_sweep_files:
    if 'Colombia' in sweep_file:
        pass
    else:
        df_sweep = pd.read_csv(general_paths_sweeps + sweep_file, sep='\t')
        print(f"Processing {sweep_file}...")
        # calculating the threshold for LR
        # calculate threshold
        threshold_99 = df_sweep['LR'].quantile(0.99)
        # filtering df_sweep based on threshold
        df_sweep = df_sweep[df_sweep['LR'] >= threshold_99]

        result = overlapping_sweeps_TEs(df_TEs_genes, df_sweep)
        if not result.empty:
            all_results.append(result)
            print(f"Found {len(result)} overlaps in {sweep_file}")

# Combine all results
if all_results:
    final_results = pd.concat(all_results, ignore_index=True)
    
    # Save to CSV
    output_path = 'TEs_overlapping_sweeps.csv'
    final_results.to_csv(output_path, index=False)
    print(f"\nSaved {len(final_results)} total overlaps to {output_path}")
else:
    print("No overlaps found")