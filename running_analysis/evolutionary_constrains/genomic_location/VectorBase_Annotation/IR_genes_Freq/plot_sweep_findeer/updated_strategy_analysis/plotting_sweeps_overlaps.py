import pandas as pd
import matplotlib.pyplot as plt
import os
import re

def plot_sweep_region(exons, df_TEs_genes, list_sweep_files):
    # for row in df_TEs_genes.itertuples():
    for index, TEs in df_TEs_genes.iterrows():
        # information TEs
        TE_chrom = TEs['chromosome']
        TE_start = TEs['start']
        TE_end = TEs['end']
        TE_name = TEs['te_name']
        gene_name = TEs['gene_name']
        
        # get gene_info
        gene_info = TEs['gene_info']
        
        # country with highest percentage - compare all columns starting with {country}_percentage
        max_country = df_TEs_genes.filter(like='_percentage').idxmax(axis=1)[index]
        max_percentage = df_TEs_genes[max_country][index]
        max_country = max_country.replace('_percentage', '')

        # condition to see if gene_info is nan
        if (gene_info != gene_info):
            # use TE position plus minus 10kb
            region_start = TE_start - 10000
            region_end = TE_end + 10000
            
        else:
            # make list with the exon information in gene_info
            exons_in_gene = exons[exons['gene_id'] == gene_info]
            # get information exons 
            exons_information = {}
            for exon_index, exon_row in exons_in_gene.iterrows():
                exon_start = exon_row['start']
                exon_end = exon_row['end']
                exons_information[exon_index] = (exon_start, exon_end)

            # get the region to plot based on the exons
            # min start and max end
            region_start = min([start for start, end in exons_information.values()]) - 10000
            region_end = max([end for start, end in exons_information.values()]) + 10000
                        
        # get the items in the sweep list that have the chromosome
        dict_information_sweeps = {}
        specific_files = [sweep_file for sweep_file in list_sweep_files if TE_chrom in sweep_file]
        for sweep_file in specific_files:
            if 'Colombia' in sweep_file:
                pass
            else:
                df_sweep = pd.read_csv(general_paths_sweeps + sweep_file, sep='\t')
                print(f"Processing {sweep_file}...")
                # get region of the df_sweep
                df_sweep_region = df_sweep[(df_sweep['location'] >= region_start) & (df_sweep['location'] <= region_end)]
                pop_code = re.search(r'AaegL5_\d_no_relatives_(.+)_raw.+', sweep_file).group(1)
                # get the ones overlapping with the region - threshold 99
                # calculate threshold
                threshold_99 = df_sweep['LR'].quantile(0.99)
                df_sweep_t = df_sweep[df_sweep['LR'] >= threshold_99]
                df_sweep_region_threshold = df_sweep_t[(df_sweep_t['location'] >= region_start) & (df_sweep_t['location'] <= region_end)]
                dict_information_sweeps[pop_code] = {'full_region': df_sweep_region, 'threshold_region': df_sweep_region_threshold}

        # CHECK: Only plot if max_country has sweep above threshold INSIDE the TE region (not extended region)
        if max_country not in dict_information_sweeps:
            print(f"Skipping {TE_name}: {max_country} (highest frequency) has no sweep data")
            continue
        
        # Filter threshold_region to only include sweeps INSIDE the TE boundaries
        sweeps_inside_TE = dict_information_sweeps[max_country]['threshold_region'][
            (dict_information_sweeps[max_country]['threshold_region']['location'] >= TE_start) &
            (dict_information_sweeps[max_country]['threshold_region']['location'] <= TE_end)
        ]
        
        if sweeps_inside_TE.empty:
            print(f"Skipping {TE_name}: {max_country} (highest frequency) has no sweep above threshold INSIDE TE")
            continue

        # plot the regions with sweeps - LR values full_region and different color for threshold_region
        plt.figure(figsize=(50, 20))
        for pop_code, sweep_data in dict_information_sweeps.items():
            plt.plot(sweep_data['full_region']['location'], sweep_data['full_region']['LR'], label=f'{pop_code} Full Region', alpha=0.5)
            plt.scatter(sweep_data['threshold_region']['location'], sweep_data['threshold_region']['LR'], color='red', label=f'{pop_code} Threshold Region')
        
        # add the region of the TE
        plt.axvspan(TE_start, TE_end, color='yellow', alpha=0.3, label='TE Region')
        # add the exon regions if available
        if gene_info == gene_info:
            for exon_index, (exon_start, exon_end) in exons_information.items():
                plt.axvspan(exon_start, exon_end, color='pink', alpha=0.3, label='Exon Region' if exon_index == 0 else "")
        # legend outside of the plot
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.xlabel('Location')
        plt.ylabel('LR')
        plt.title(f'{TE_name} on {gene_name} - Max country: {max_country} ({max_percentage:.1f}%)')
        plt.savefig(f'./plots_overlap/sweep_region_TE_{TE_name}_{TE_chrom}.png')
        plt.close()

def overlapping_sweeps_TEs(df_TEs_genes, df_sweep):
    # check for overlap between sweeps and TEs with overlap status
    overlapping_TEs = []
    for index, sweep_row in df_sweep.iterrows():
        sweep_position = sweep_row['location']
        sweep_chrom = sweep_row['chrom']
        sweep_LR = sweep_row['LR']  # Capture LR value
        sweep_pop = sweep_row['pop']  # Capture pop value
        
        #overlapping = df_TEs_genes[
        #    (df_TEs_genes['chromosome'] == sweep_chrom) &
        #    (df_TEs_genes['start'] - 10000 <= sweep_position) &
        #    (df_TEs_genes['end'] + 10000 >= sweep_position)
        #]

        overlapping = df_TEs_genes[
            (df_TEs_genes['chromosome'] == sweep_chrom) &
            (df_TEs_genes['start'] <= sweep_position) &
            (df_TEs_genes['end']  >= sweep_position)
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
exons = exons.copy()
exons['gene_id'] = exons['attributes'].str.extract(r'gene_id=([^;]+)')

# total each country
count_country = {'Brazil': 18, 'Colombia': 24, 'USA': 28, 'Senegal': 20, 'Gabon': 13, 'Kenya': 19}
df_TEs_genes = pd.read_csv('./genes_overlap_tes_ref_with_GO_noduplicate.csv')

# extract from gene name
df_TEs_genes['gene_info'] = df_TEs_genes['gene_name'].str.extract(r'ID=([^;]+)')
# add percentages of each TE - using the columns that start with # and dividing by the total number of TEs
df_TEs_genes = df_TEs_genes.copy()
for country, total in count_country.items():
    df_TEs_genes[f'{country}_percentage'] = df_TEs_genes.filter(like=f'count_1/1_{country}').sum(axis=1) / total * 100

# open new dfs for sweeps
general_paths_sweeps = './make_new_sweep_files/new_sweeps_files/'
list_sweep_files = os.listdir(general_paths_sweeps)

plot_sweep_region(exons, df_TEs_genes, list_sweep_files)

'''all_results = []  # Store results from all sweep files

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
    print("No overlaps found")'''