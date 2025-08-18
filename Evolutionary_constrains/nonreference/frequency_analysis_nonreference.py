import pandas as pd
import matplotlib.pyplot as plt
import re
import polars as pl
import os
import numpy as np
from itertools import combinations
import seaborn as sns

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
    country_name - Name of the country
    
    Returns:
    DataFrame with frequency column added
    """
    # Get columns that contain genotype data (excluding base columns)
    base_cols = ['chromosome', 'position_start', 'position_end', 'TE_name', 'strand']
    genotype_cols = [col for col in df.columns if col not in base_cols]
    
    # Count '1/1' occurrences in each row
    df[f'count_1/1_{country_name}'] = df[genotype_cols].apply(lambda row: (row == '1/1').sum(), axis=1)
    
    # Calculate frequency (proportion of samples with 1/1)
    total_samples = len(genotype_cols)
    df[f'frequency_{country_name}'] = df[f'count_1/1_{country_name}'] / total_samples
    
    return df


nonreference_path = '/work/users/g/a/gabivla/lab/SV_mosquitos/VCF_genotypes_TEs/separate_chromosomes_windows_reference/merged_nonreference/VCF_nonreference_post_processing.csv'
neutral_data_path = '/nas/longleaf/home/gabivla/lab/SV_mosquitos/neutral_data/making_SFS/my_script/homozygous_counts_by_country/SFS_files'
df = pd.read_csv(nonreference_path)
# remove first column
df = df.iloc[:, 1:]

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
    
    # Get TEs with top 10 frequency differences
    top_10_diff = merged_df.sort_values("frequency_diff", ascending=False).head(10)
    all_top10_diffs[(country1, country2)] = top_10_diff

    # plot distribution of the frequency differences
    plt.hist(merged_df["frequency_diff"], bins=50)
    plt.xlabel("Frequency difference")
    plt.ylabel("Frequency")
    plt.title(f"Distribution of frequency differences between {country1} and {country2}")
    plt.savefig(f"./results/frequency_diff_distribution_{country1}_{country2}.png")
    plt.close()

    # Save top 10 differences to a CSV file
    top_10_diff.to_csv(f"./results/top_10_diff_{country1}_{country2}.csv", index=False)

    # Plot frequency differences
    plt.figure(figsize=(30, 15))
    
    # Create bar plot for each country
    bar_width = 0.4
    x_labels = top_10_diff["TE_information"]
    x_indexes = range(len(x_labels))

    plt.bar(x_indexes, top_10_diff[f'frequency_{country1}'], width=bar_width, label=country1, color="blue")
    plt.bar([x + bar_width for x in x_indexes], top_10_diff[f'frequency_{country2}'], width=bar_width, label=country2, color="red")

    plt.xlabel("Transposable Element (TE)")
    plt.ylabel("Frequency")
    plt.title(f"Top 10 TEs with highest frequency differences between {country1} and {country2}")
    plt.xticks([x + bar_width / 2 for x in x_indexes], x_labels, rotation=90)
    plt.legend()
    
    # Save bar plot
    plt.tight_layout()
    plt.savefig(f"./results/top_10_TEs_frequencies_{country1}_{country2}.png")
    plt.close()

# merge all top 10 differences
all_top10_diffs_df = pd.concat(all_top10_diffs.values())
# drop all columns that contain frequency_* in the header
all_top10_diffs_df = all_top10_diffs_df.drop(columns=all_top10_diffs_df.filter(like="frequency_").columns)
# get the information from the TE_information column in the original df - all_dfs_ref
for col in all_top10_diffs_df.columns:
    if col.startswith("count_1/1_"):  # Identify relevant columns
        country = col.split("count_1/1_")[1]  # Extract country name
        
        if country in all_dfs_ref:  # Ensure the country exists in all_dfs_ref
            country_df = all_dfs_ref[country]  # Get the corresponding DataFrame
            # Map TE_information from all_top10_diffs_df to the count_1/1 column in country_df
            fill_values = country_df.set_index("TE_information")["count_1/1_" + country]
            # Fill NaN values with the corresponding values from all_dfs_ref
            all_top10_diffs_df[col] = all_top10_diffs_df[col].fillna(all_top10_diffs_df["TE_information"].map(fill_values))

# make table with frequencies for all countries
all_top10_diffs_df = all_top10_diffs_df.sort_values(by="TE_information").reset_index(drop=True)
all_top10_diffs_df.to_csv('all_top10_diffs_all_countries_nonreference.csv', index=False)

df = all_top10_diffs_df

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Define country groups
group1 = ['count_1/1_Kenya', 'count_1/1_Senegal', 'count_1/1_Gabon']
group2 = ['count_1/1_USA', 'count_1/1_Brazil', 'count_1/1_Colombia']

# Calculate mean frequency for each country in both groups
group1_means = df[group1].mean()
group2_means = df[group2].mean()

# Remove 'count_1/1_' from the labels for plotting
group1_labels = [country.replace('count_1/1_', '') for country in group1]
group2_labels = [country.replace('count_1/1_', '') for country in group2]

# Plot setup for mean frequencies
plt.figure(figsize=(10, 6))

# Plot mean frequencies for Group 1
plt.bar(group1_labels, group1_means.values, color=['blue', 'green', 'orange'], label='Group 1 (Kenya, Senegal, Gabon)', alpha=0.7)

# Plot mean frequencies for Group 2
plt.bar(group2_labels, group2_means.values, color=['brown', 'purple', 'red'], label='Group 2 (USA, Brazil, Colombia)', alpha=0.7)

# Labels and title
plt.xlabel('Country', fontsize=12)
plt.ylabel('Mean Frequency', fontsize=12)
plt.title('Mean NonReference TE Insertion Frequencies by Country', fontsize=14)
plt.legend(title="Country Group")

# Show plot
plt.tight_layout()
plt.savefig('mean_frequencies_nonreference.svg')

# Correlation calculation between all countries
correlation_matrix = df[group1 + group2].corr()

# Rename columns and index in the correlation matrix for plotting
correlation_matrix.index = [x.replace('count_1/1_', '') for x in correlation_matrix.index]
correlation_matrix.columns = [x.replace('count_1/1_', '') for x in correlation_matrix.columns]

# Plot the correlation matrix as a heatmap
plt.figure(figsize=(8, 6))
sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', fmt='.2f', cbar=True, linewidths=0.5)
plt.title('Correlation Matrix of Nonreference TE Frequencies Across Countries', fontsize=14)
plt.tight_layout()
plt.savefig('correlation_matrix_nonreference.svg')

# Print the correlation matrix for reference
print("Correlation Matrix for All Countries:")
print(correlation_matrix)
