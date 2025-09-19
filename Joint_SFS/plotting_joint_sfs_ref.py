import pandas as pd
import matplotlib.pyplot as plt
from utils import create_joint_sfs_subplot
import re
import polars as pl
import os
import numpy as np
import itertools
import allel

def count_TE_occurence(general_df):
    """
    Counts the occurrences of TEs and identifies the top 5 most frequent ones.

    Args:
        general_df (pl.DataFrame): Polars DataFrame containing TE data.

    Returns:
        pl.DataFrame: DataFrame with the top 5 TEs sorted by their occurrence count.
    """
    general_df = general_df.with_columns(pl.col("TE_name").str.to_uppercase())

    df_te = general_df.with_columns(pl.col("TE_name").str.extract(r'RND_\d+_FAMILY_\d+_([A-Z]+(?:_HELITRON)?)', 1).alias("TE_name"))

    grouped = df_te.group_by('TE_name').len()
    sorted_grouped = grouped.sort('len',  descending=True)
    
    return sorted_grouped

def calculate_frequency_TE(df, country_name):
    """
    Calculate the frequency of '1/1' genotypes for each TE insertion.
    
    Parameters:
    df - DataFrame with TE data
    country_name - Name of the country (for reference, not used in calculation)
    
    Returns:
    DataFrame with frequency column added
    """
    base_cols = ['chromosome', 'position_start', 'position_end', 'TE_name', 'strand']
    genotype_cols = [col for col in df.columns if col not in base_cols]
    
    df[f'count_1/1_{country_name}'] = df[genotype_cols].apply(lambda row: (row == '1/1').sum(), axis=1)
    df[f'count_0/0_{country_name}'] = df[genotype_cols].apply(lambda row: (row == '0/0').sum(), axis=1)
    total_samples = df[f'count_1/1_{country_name}'] + df[f'count_0/0_{country_name}']
    df = df.drop(columns=[f'count_0/0_{country_name}'])
    df[f'frequency_{country_name}'] = df[f'count_1/1_{country_name}'] / total_samples
    
    return df

reference_path = '/work/users/g/a/gabivla/lab/SV_mosquitos/VCF_genotypes_TEs/separate_chromosomes_windows_reference/merged_reference/VCF_reference_TEs_post_processing_600_absence_evidence.csv'

neutral_data_path = '/nas/longleaf/home/gabivla/lab/SV_mosquitos/neutral_data/making_SFS/my_script/homozygous_counts_by_country/SFS_files'
df_ref = pd.read_csv(reference_path)

# Read only necessary columns from location file
location_path = '/nas/longleaf/home/gabivla/lab/SV_mosquitos/location_information/mcclintock_output_locations.csv'
location_df = pd.read_csv(location_path, usecols=['Genome_code', 'country'])
location_dict = dict(zip(location_df['Genome_code'], location_df['country']))

base_cols = ['chromosome', 'position_start', 'position_end', 'TE_name', 'strand']

genome_cols_ref = [col for col in df_ref.columns if col not in base_cols]

rename_dict_ref = {col: f"{col}_{location_dict.get(col, '')}" for col in genome_cols_ref}
df_ref = df_ref.loc[:, ~df_ref.columns.str.contains("Unnamed")]
df_ref = df_ref.rename(columns=rename_dict_ref)

countries_reference = set()
for original_col in genome_cols_ref:
    new_col = rename_dict_ref[original_col]
    country = new_col.split('_')[-1]
    countries_reference.add(country)

# reference TEs
df_reference = pl.DataFrame(df_ref)
TE_types_ref = count_TE_occurence(df_reference)
TE_types_ref = TE_types_ref['TE_name'].to_list()
df_ref = df_reference.to_pandas()

population_sizes = {'Brazil': 19, 'Gabon': 13, 'Colombia': 24, 'USA': 28, 'Senegal': 20, 'Kenya': 19}

n_rows = 5  
n_cols = 3  
fig, axes = plt.subplots(n_rows, n_cols, figsize=(30, 30))
axes = axes.flatten()  

df_countries = {}

for i, country in enumerate(countries_reference):

    country_cols = [col for col in df_ref.columns 
                   if col not in base_cols and col.endswith('_' + country)]
    
    if not country_cols:
        continue
    
    country_df = df_ref[base_cols + country_cols]
    
    country_df = country_df.dropna(subset=country_cols, how='all')
    
    print(f"\n--- Data for {country} ---")
    
    frequency_TEs = calculate_frequency_TE(country_df, country)
    df_countries[country] = frequency_TEs

country_list = list(df_countries.keys())
country_pairs = list(itertools.combinations(country_list, 2))

num_pairs = len(country_pairs)
cols = min(5, num_pairs)  
rows = (num_pairs + cols - 1) // cols  

fig, axes = plt.subplots(rows, cols, figsize=(4 * cols, 4 * rows))
axes = np.array(axes) 

axes = axes.flatten()

for i, (country_1, country_2) in enumerate(country_pairs):
    df_country_1 = df_countries[country_1]
    df_country_2 = df_countries[country_2]

    df_country_1['TE_information'] = df_country_1['chromosome'].astype(str) + '_' + df_country_1['position_start'].astype(str) + '_' + df_country_1['position_end'].astype(str) + '_' + df_country_1['TE_name'] + '_' + df_country_1['strand']
    df_country_2['TE_information'] = df_country_2['chromosome'].astype(str) + '_' + df_country_2['position_start'].astype(str) + '_' + df_country_2['position_end'].astype(str) + '_' + df_country_2['TE_name'] + '_' + df_country_2['strand']
    
    df_country_1 = df_country_1.iloc[:, list(range(0)) + list(range(-3, 0))]
    df_country_2 = df_country_2.iloc[:, list(range(0)) + list(range(-3, 0))]
    
    merged_df = df_country_1.merge(df_country_2, on="TE_information", how="outer", suffixes=("_country1", "_country2"))
    merged_df = merged_df[(merged_df[f'count_1/1_{country_1}'] != 0) | (merged_df[f'count_1/1_{country_2}'] != 0)]
    
    cax = create_joint_sfs_subplot(axes[i], country_1, country_2, merged_df)

for j in range(i + 1, len(axes)):
    fig.delaxes(axes[j])

fig.subplots_adjust(right=0.85)
cbar_ax = fig.add_axes([0.87, 0.2, 0.02, 0.6])
fig.colorbar(cax, cax=cbar_ax, label="Number of TEs (log scale)")

plt.tight_layout(rect=[0, 0, 0.85, 1])
plt.savefig('joint_sfs_all_countries_reference.svg')
plt.close()