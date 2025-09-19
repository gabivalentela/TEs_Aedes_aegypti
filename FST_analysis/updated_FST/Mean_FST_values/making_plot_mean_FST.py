# making plot with mean FST values for each combination of countries
import os 
import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def pair_data(list_files_tes, list_files_snps, TEs_data, SNPs_data):
    dict_data = {}
    
    for file in list_files_tes:
        # Construct the full path for the TE file
        path_TE = os.path.join(TEs_data, file)
        
        # Extract country1, country2, and pair type using regex
        country_1 = re.search(r'FST_values_(.+)_(.+)_.+.csv', file).group(1)
        country_2 = re.search(r'FST_values_(.+)_(.+)_.+.csv', file).group(2)
        pair_type = re.search(r'FST_values_(.+)_(.+)_(.+).csv', file).group(3)
        
        # Sort the countries alphabetically to ensure consistent key order
        sorted_countries = tuple(sorted([country_1, country_2]))
        print(f"Processing pair: {sorted_countries} with pair type: {pair_type}")
        
        # Find the corresponding SNP file
        for SNPs in list_files_snps:
            if country_1 in SNPs and country_2 in SNPs:
                # Construct the full path for the SNP file
                snp_path = os.path.join(SNPs_data, SNPs)
                
                # Initialize an entry for this country combination if not already present
                if sorted_countries not in dict_data:
                    dict_data[sorted_countries] = {}
                
                # Add the pair_type data to the dictionary
                dict_data[sorted_countries][pair_type] = [path_TE, snp_path]
    
    return dict_data

def get_mean_FST_values(ref_file):
    df = pd.read_csv(ref_file)
    # get mean value for FST
    mean_FST = df['FST'].mean()
    return mean_FST

def get_mean_FST_values_SNPs(ref_file):
    df = pd.read_csv(ref_file)
    # get mean value for FST
    mean_FST = df['1/1'].mean()
    return mean_FST

TEs_data = '/work/users/g/a/gabivla/lab/SV_mosquitos/FST_analysis/TEs_absence_evidence'
SNPs_data = '/work/users/g/a/gabivla/lab/SV_mosquitos/FST_analysis/SNPs'

# Load data
list_files_tes = os.listdir(TEs_data)
list_files_snps = os.listdir(SNPs_data)
files_SNPs_TEs = pair_data(list_files_tes, list_files_snps, TEs_data, SNPs_data)

dict_mean_FST = {}

for country, values in files_SNPs_TEs.items():
    country_1, country_2 = country
    # get ref files 
    ref_files = values['ref']
    nonref_files = values['nonref'][0]
    nonref_FST = get_mean_FST_values(nonref_files)
    for file in ref_files:
        if 'ref' in file:
            ref_file = file
            ref_FST = get_mean_FST_values(ref_file)
        else:
            SNP_file = file
            snp_FST = get_mean_FST_values_SNPs(SNP_file)

    dict_mean_FST[country] = {'ref': ref_FST, 'nonref': nonref_FST, 'SNPs': snp_FST}

# Get all unique countries
countries = sorted(set([c for pair in dict_mean_FST for c in pair]))

# Desired order of countries
country_order = ['Brazil', 'Colombia', 'USA', 'Gabon', 'Kenya', 'Senegal']
# Initialize empty dataframes for each type with desired order
ref_df = pd.DataFrame(0.0, index=country_order, columns=country_order)
nonref_df = pd.DataFrame(0.0, index=country_order, columns=country_order)
snps_df = pd.DataFrame(0.0, index=country_order, columns=country_order)

# Fill the dataframes
for (c1, c2), values in dict_mean_FST.items():
    for df, key in zip([ref_df, nonref_df, snps_df], ['ref', 'nonref', 'SNPs']):
        df.loc[c1, c2] = values[key]
        df.loc[c2, c1] = values[key] 

# Plot heatmaps
import numpy as np 

# Plot heatmaps
for name, df in zip(['ref', 'nonref', 'SNPs'], [ref_df, nonref_df, snps_df]):
    
    mask = np.triu(np.ones_like(df, dtype=bool))

    plt.figure(figsize=(8, 6))
    sns.heatmap(
        df, annot=True, cmap='plasma', square=True, 
        cbar_kws={'label': f'{name} FST'},
        mask=mask, linewidths=0.5, linecolor='white',
        vmin=0, vmax=0.15  
    )
    plt.title(f'Pairwise FST - {name}')
    plt.xlabel('Country')
    plt.ylabel('Country')
    plt.tight_layout()
    plt.savefig(f'heatmap_{name}.svg')
