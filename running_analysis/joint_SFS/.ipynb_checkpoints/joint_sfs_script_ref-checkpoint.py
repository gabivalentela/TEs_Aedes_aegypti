import pandas as pd
import matplotlib.pyplot as plt
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
    # Get columns that contain genotype data (excluding base columns)
    base_cols = ['chromosome', 'position_start', 'position_end', 'TE_name']
    genotype_cols = [col for col in df.columns if col not in base_cols]
    
    # Count '1/1' occurrences in each row
    df[f'count_1/1_{country_name}'] = df[genotype_cols].apply(lambda row: (row == '1/1').sum(), axis=1)
    
    # Calculate frequency (proportion of samples with 1/1)
    total_samples = len(genotype_cols)
    df[f'frequency_{country_name}'] = df[f'count_1/1_{country_name}'] / total_samples
    
    return df

reference_path = '/work/users/g/a/gabivla/lab/SV_mosquitos/VCF_genotypes_TEs/reference/subset/VCF_merging_tes_subset_ref.csv'

neutral_data_path = '/nas/longleaf/home/gabivla/lab/SV_mosquitos/neutral_data/making_SFS/my_script/homozygous_counts_by_country/SFS_files'
df_ref = pd.read_csv(reference_path)

# Read only necessary columns from location file
location_path = '/nas/longleaf/home/gabivla/lab/SV_mosquitos/location_information/mcclintock_output_locations.csv'
location_df = pd.read_csv(location_path, usecols=['Genome_code', 'country'])
location_dict = dict(zip(location_df['Genome_code'], location_df['country']))

# Define base columns that should not be modified
base_cols = ['chromosome', 'position_start', 'position_end', 'TE_name', 'strand']

# Get genome columns at once
genome_cols_ref = [col for col in df_ref.columns if col not in base_cols]

# Create rename dictionary and apply it once
rename_dict_ref = {col: f"{col}_{location_dict.get(col, '')}" for col in genome_cols_ref}
# drop column that is being a problem - only for reference
df_ref = df_ref.loc[:, ~df_ref.columns.str.contains("Unnamed")]
df_ref = df_ref.rename(columns=rename_dict_ref)

# Extract unique countries from the renamed genome columns only - reference
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

# Calculate the number of subplots needed
n_rows = 5  
n_cols = 3  
# Create the figure and axes for subplots
fig, axes = plt.subplots(n_rows, n_cols, figsize=(40, 40))
axes = axes.flatten()  # Flatten the axes array for easy indexing

df_countries = {}

# loop for reference
# Split data by country
for i, country in enumerate(countries_reference):

    # Get all columns for this country from the renamed genome columns
    country_cols = [col for col in df_ref.columns 
                   if col not in base_cols and col.endswith('_' + country)]
    
    if not country_cols:
        continue
    
    # Create a subset with base columns and columns for this country
    country_df = df_ref[base_cols + country_cols]
    
    # Drop rows where all country columns are NA
    country_df = country_df.dropna(subset=country_cols, how='all')
    
    print(f"\n--- Data for {country} ---")
    
    # calculate frequency of each TE for each country
    frequency_TEs = calculate_frequency_TE(country_df, country)
    df_countries[country] = frequency_TEs

country_list = list(df_countries.keys())
# make combination of countries
combinations = itertools.combinations(country_list, 2)
# Save combinations as a list of tuples
country_pairs = list(combinations)
# make df for combination of countries
for i, country_pair in enumerate(country_pairs):
    country_1 = country_pair[0]
    country_2 = country_pair[1]
    df_country_1 = df_countries[country_1]
    df_country_2 = df_countries[country_2]

    # make TE_information column - merge info from chromosome, position_start, position_end, TE_name, strand - with '_'\
    df_country_1['TE_information'] = df_country_1['chromosome'].astype(str) + '_' + df_country_1['position_start'].astype(str) + '_' + df_country_1['position_end'].astype(str) + '_' + df_country_1['TE_name'] + '_' + df_country_1['strand']
    df_country_2['TE_information'] = df_country_2['chromosome'].astype(str) + '_' + df_country_2['position_start'].astype(str) + '_' + df_country_2['position_end'].astype(str) + '_' + df_country_2['TE_name'] + '_' + df_country_2['strand']
    
    # drop genotype columns - after the 5 column and leave country_1/1 and frequency_country
    df_country_1 = df_country_1.iloc[:, list(range(0)) + list(range(-3, 0))]
    df_country_2 = df_country_2.iloc[:, list(range(0)) + list(range(-3, 0))]
    
    # Merge dataframes on base columns
    merged_df = df_country_1.merge(df_country_2, on="TE_information", how="outer", suffixes=("_country1", "_country2"))
    merged_df = merged_df[["TE_information"] + [col for col in merged_df.columns if col != "TE_information"]]

    # if both count_1/1 column are 0 - drop row
    merged_df = merged_df[(merged_df[f'count_1/1_{country_1}'] != 0) | (merged_df[f'count_1/1_{country_2}'] != 0)]

    # use this data to calculate the joint SFS
    # make list with count 1/1 columns
    df1 = merged_df[f'count_1/1_{country_1}'].tolist()
    df2 = merged_df[f'count_1/1_{country_2}'].tolist()

    # Check if the lists have the same length
    if len(df1) != len(df2):
        raise ValueError(f"The lists 'brazil' and 'gabon' must have the same length for file")

    # Number of individuals in each population
    n1 = population_sizes[country_1]  # Brazil
    n2 = population_sizes[country_2]  # Gabon

    # Compute the joint site frequency spectrum
    jsfs = allel.joint_sfs(df1, df2, n1=n1, n2=n2)
    print(jsfs)

    # Plot the joint site frequency spectrum with country names on the axes
    ax = axes[i]
    allel.plot_joint_sfs(jsfs, ax=ax)

    # Set axis labels with country names
    ax.set_xlabel(country_1, fontsize=30)
    ax.set_ylabel(country_2, fontsize=30)
    #change font size
    ax.tick_params(axis='both', which='major', labelsize=30)
    ax.set_title(f"JSFS for {country_1} and {country_2}", fontsize=30)

# Adjust layout for better spacing
plt.tight_layout()
plt.savefig('joint_sfs_reference.png')
