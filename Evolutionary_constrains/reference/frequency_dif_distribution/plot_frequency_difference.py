import pandas as pd
import matplotlib.pyplot as plt
import re
import polars as pl
import os
import numpy as np
from itertools import combinations
import seaborn as sns
import gzip

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
    base_cols = ['chromosome', 'position_start', 'position_end', 'TE_name', 'strand']
    genotype_cols = [col for col in df.columns if col not in base_cols]
    
    # Count '1/1' occurrences in each row
    df[f'count_1/1_{country_name}'] = df[genotype_cols].apply(lambda row: (row == '1/1').sum(), axis=1)
    # count 0/0 occurrences in each row
    df[f'count_0/0_{country_name}'] = df[genotype_cols].apply(lambda row: (row == '0/0').sum(), axis=1)

    # Calculate frequency (proportion of samples with 1/1)
    total_samples = df[f'count_1/1_{country_name}'] + df[f'count_0/0_{country_name}']
    # remove count_0/0_{country_name} column
    df = df.drop(columns=[f'count_0/0_{country_name}'])
    df[f'frequency_{country_name}'] = df[f'count_1/1_{country_name}'] / total_samples
    
    return df

def read_vcf(vcf_file):
    with gzip.open(vcf_file, 'rt') as f:  # 'rt' mode for reading text files
        # Initialize header as None before finding it
        header = None
        # Read through the file and get the header line
        for line in f:
            if line.startswith("#CHROM"):
                header = line.strip("#").strip().split("\t")
                break  # Exit the loop once the header is found
        
        # Read the rest of the VCF file into a DataFrame, skipping comment lines
        df = pd.read_csv(f, sep='\t', comment='#', header=None, names=header)

    # Add a new 'info' column based on 'CHROM' and 'POS'
    df['info'] = df['CHROM'] + '_' + df['POS'].astype(str)
    
    # Keep only the portion before the first ":" in all columns except 'info'
    #df = df.applymap(lambda x: x.split(':')[0] if isinstance(x, str) else x)
    df = df.map(lambda x: x.split(':')[0] if isinstance(x, str) else x)

    # Drop unwanted columns
    columns_to_remove = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'FORMAT', 'INFO']
    df_selected = df.drop(columns=columns_to_remove)
    
    return df_selected

def info_location(loc_info):
    """
    Reads the location information from a CSV file into a Polars DataFrame.

    Args:
        loc_info (str): Path to the CSV file containing location information.

    Returns:
        pl.DataFrame: Polars DataFrame with location information.
    """
    df = pd.read_csv(loc_info)
    return df

def split_by_countries(df, loc_info):
    # Create a dictionary with the genome code as the key and the country as the value
    loc_info_dict = loc_info.set_index('Genome_code')['country'].to_dict()

    # Prepare the dictionary to store results
    country_dict = {}
    
    # Iterate over each column in the dataframe (each column corresponds to a genome)
    for col in df.columns:
        # Get the genome code from the column name (assuming column names are genome codes)
        genome_code = col
        
        # Get the country for the current genome
        country = loc_info_dict.get(genome_code, 'Cali')  # Default to 'Unknown' if the genome code isn't in loc_info_dict
        
        # Add a new column to the dictionary for this country if it doesn't already exist
        if country not in country_dict:
            country_dict[country] = df.loc[:, [col]]  # Add the entire column for this country
        else:
            country_dict[country] = pd.concat([country_dict[country], df.loc[:, [col]]], axis=1)  # Add the column to the existing country data

    return country_dict

def update_SNPs_VCF(df):
    # Exclude the "info" column from genotype processing
    genotype_cols = df.columns.difference(["info"])

    # Count occurrences of each genotype per row
    counts = df[genotype_cols].apply(lambda row: row.value_counts(), axis=1).fillna(0)

    # Get number of `1/1` and `0/0` per row
    count_11 = counts.get("1/1", 0)  # Count occurrences of 1/1
    count_00 = counts.get("0/0", 0)  # Count occurrences of 0/0

    # Filter: Remove rows with only one `1/1` (singleton) or with only 1 or 2 `0/0`
    filtered_df = df[~((count_11 == 1) | (count_00 <= 1))]

    return filtered_df


reference_path = '/work/users/g/a/gabivla/lab/SV_mosquitos/VCF_genotypes_TEs/separate_chromosomes_windows_reference/merged_reference/VCF_reference_TEs_post_processing_600.csv'
neutral_data_path = '/nas/longleaf/home/gabivla/lab/SV_mosquitos/neutral_data/making_SFS/my_script/homozygous_counts_by_country/SFS_files'
df = pd.read_csv(reference_path)
# remove first column
df = df.iloc[:, 1:]

# read VCF SNPs
vcf_SNPs = '/work/users/g/a/gabivla/lab/SV_mosquitos/neutral_data_vcf/updated_AaegL5_full_filtered_snps_110122_rep_map_100kbgenes_centromere_filtered_sites_mask.vcf.gz'
df_SnpS = read_vcf(vcf_SNPs)
df_SnpS = update_SNPs_VCF(df_SnpS)

#location info
loc_info = info_location('/nas/longleaf/home/gabivla/lab/SV_mosquitos/location_information/mcclintock_output_locations.csv')

vcf_files = split_by_countries(df_SnpS, loc_info)

# Read only necessary columns from location file
location_path = '/nas/longleaf/home/gabivla/lab/SV_mosquitos/location_information/mcclintock_output_locations.csv'
location_df = pd.read_csv(location_path, usecols=['Genome_code', 'country'])
location_dict = dict(zip(location_df['Genome_code'], location_df['country']))

# Define base columns that should not be modified
base_cols = ['chromosome', 'position_start', 'position_end', 'TE_name', 'strand']

# Get genome columns at once
genome_cols = [col for col in df.columns if col not in base_cols]

# Create rename dictionary and apply it once
rename_dict = {col: f"{col}_{location_dict.get(col, '')}" for col in genome_cols}
# drop column that is being a problem
df = df.loc[:, ~df.columns.str.contains("Unnamed")]
df = df.rename(columns=rename_dict)

# Extract unique countries from the renamed genome columns only
countries = set()
for original_col in genome_cols:
    new_col = rename_dict[original_col]
    country = new_col.split('_')[-1]
    countries.add(country)

df = pl.DataFrame(df)
TE_types = count_TE_occurence(df)
TE_types = TE_types['TE_name'].to_list()
df = df.to_pandas()

country_pairs_ref = list(combinations(countries, 2))
all_dfs_ref = {}
country_vcf_dict_ref = {}

# Split data by country
for i, country in enumerate(countries):

    # Get all columns for this country from the renamed genome columns
    country_cols = [col for col in df.columns 
                   if col not in base_cols and col.endswith('_' + country)]
    
    if not country_cols:
        continue
    
    # Create a subset with base columns and columns for this country
    country_df = df[base_cols + country_cols]
    
    # Drop rows where all country columns are NA
    country_df = country_df.dropna(subset=country_cols, how='all')
    
    print(f"\n--- Data for {country} ---")
    # calculate frequency
    frequency_df = calculate_frequency_TE(country_df, country)
    all_dfs_ref[country] = frequency_df

for country1, country2 in country_pairs_ref:
    country_vcf_dict_ref[(country1, country2)] = (all_dfs_ref[country1], all_dfs_ref[country2])

all_top10_diffs = {}

# Iterate over country pairs and their corresponding VCF data
for i, (pair, vcf_tuple) in enumerate(country_vcf_dict_ref.items()):
    country1, country2 = pair
    df_country_1, df_country_2 = vcf_tuple

   # make TE_information column - merge info from chromosome, position_start, position_end, TE_name, strand - with '_'\
    df_country_1['TE_information'] = df_country_1['chromosome'].astype(str) + '_' + df_country_1['position_start'].astype(str) + '_' + df_country_1['position_end'].astype(str) + '_' + df_country_1['TE_name'] + '_' + df_country_1['strand']
    df_country_2['TE_information'] = df_country_2['chromosome'].astype(str) + '_' + df_country_2['position_start'].astype(str) + '_' + df_country_2['position_end'].astype(str) + '_' + df_country_2['TE_name'] + '_' + df_country_2['strand']
    
    # drop genotype columns - after the 5 column and leave country_1/1 and frequency_country
    df_country_1 = df_country_1.iloc[:, list(range(0)) + list(range(-3, 0))]
    df_country_2 = df_country_2.iloc[:, list(range(0)) + list(range(-3, 0))]
    
    # Merge dataframes on base columns
    merged_df = df_country_1.merge(df_country_2, on="TE_information", how="outer", suffixes=("_country1", "_country2"))
    merged_df = merged_df[["TE_information"] + [col for col in merged_df.columns if col != "TE_information"]]
    
    # Calculate frequency difference
    merged_df["frequency_diff"] = abs(merged_df[f"frequency_{country1}"] - merged_df[f"frequency_{country2}"])
    
    # Count the number of each frequency difference
    diff_counts = merged_df["frequency_diff"].value_counts(normalize=True).reset_index()
    diff_counts.columns = ["frequency_diff", "percentage"]
    diff_counts["percentage"] *= 100  # Convert to percentage
    diff_counts = diff_counts.sort_values("frequency_diff")

    # get SNP info
    country_1_SNPs = vcf_files[country1]
    country_2_SNPs = vcf_files[country2]
    # calculate frequency for each row - count number of 1/1 and 0/1
    count_11_country1 = country_1_SNPs.apply(lambda row: (row == "1/1").sum(), axis=1)
    count_01_country1 = country_1_SNPs.apply(lambda row: (row == "0/1").sum(), axis=1)
    ## sum counts of 1
    count_1_country_1 = ((count_11_country1*2) + count_01_country1)/(len(country_1_SNPs.columns)*2)
    # country 2
    count_11_country2 = country_2_SNPs.apply(lambda row: (row == "1/1").sum(), axis=1)
    count_01_country2 = country_2_SNPs.apply(lambda row: (row == "0/1").sum(), axis=1)
    ## sum counts of 1
    count_1_country_2 = ((count_11_country2*2) + count_01_country2)/(len(country_2_SNPs.columns)*2)
    
    # add counts to df 
    country_1_SNPs[f'frequency_SNP_{country1}'] = count_1_country_1
    country_2_SNPs[f'frequency_SNP_{country2}'] = count_1_country_2

    # make a df with the frequency_SNPs for each country - based on the row indexes
    df_SNPs = country_1_SNPs.merge(country_2_SNPs, left_index=True, right_index=True)
    df_SNPs_frequency = df_SNPs[[f'frequency_SNP_{country1}', f'frequency_SNP_{country2}']]

    # Calculate frequency difference
    df_SNPs_frequency["frequency_diff"] = abs(df_SNPs_frequency[f"frequency_SNP_{country1}"] - df_SNPs_frequency[f"frequency_SNP_{country2}"])
    
    # Count the number of each frequency difference as a percentage
    diff_counts_SNPs = df_SNPs_frequency["frequency_diff"].value_counts(normalize=True).reset_index()
    diff_counts_SNPs.columns = ["frequency_diff", "percentage"]
    diff_counts_SNPs["percentage"] *= 100  # Convert to percentage
    diff_counts_SNPs = diff_counts_SNPs.sort_values("frequency_diff")
    
    # Create the plot
    plt.figure(figsize=(8, 6))

    # Plot TE frequency difference line in blue
    plt.plot(diff_counts["frequency_diff"], diff_counts["percentage"], linestyle="-", color="b", label="TE Frequency Difference")

    # Plot SNP frequency difference line in red
    plt.plot(diff_counts_SNPs["frequency_diff"], diff_counts_SNPs["percentage"], linestyle="-", color="r", label="SNP Frequency Difference")

    # Labels and title
    plt.xlabel("Frequency Difference")
    plt.ylabel("Count")
    plt.title(f"Frequency Difference Distribution: {country1} vs {country2}")

    # Add legend to distinguish lines
    plt.legend()

    # Save the plot
    plt.savefig(f"./results/frequency_diff_distribution_{country1}_{country2}.png")

    # Create the plot
    plt.figure(figsize=(8, 6))

    # Plot TE frequency difference histogram (blue)
    plt.hist(
        diff_counts["frequency_diff"], 
        weights=diff_counts["percentage"],  # Use percentage as weights
        bins=20,  # Adjust bin size as needed
        alpha=0.5,  # Transparency to see overlapping
        color="b", 
        label="TE Frequency Difference"
    )

    # Plot SNP frequency difference histogram (red)
    plt.hist(
        diff_counts_SNPs["frequency_diff"], 
        weights=diff_counts_SNPs["percentage"],  
        bins=20, 
        alpha=0.5, 
        color="r", 
        label="SNP Frequency Difference"
    )

    # Labels and title
    plt.xlabel("Frequency Difference")
    plt.ylabel("Percentage")
    plt.title(f"Frequency Difference Distribution: {country1} vs {country2}")

    # Add legend to distinguish histograms
    plt.legend()

    # Save the plot
    plt.savefig(f"./results_hist/frequency_diff_distribution_{country1}_{country2}.png")

    # Show the plot
    plt.show()