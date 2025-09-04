import os
import pandas as pd
from itertools import combinations
import numpy as np
import matplotlib.pyplot as plt 

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
    df['count_1/1'] = df[genotype_cols].apply(lambda row: (row == '1/1').sum(), axis=1)
    
    # Calculate frequency (proportion of samples with 1/1)
    total_samples = len(genotype_cols)
    df['frequency'] = df['count_1/1'] / total_samples
    
    return df

def SFS_data(all_dfs):
    """
    Create a dictionary of DataFrames containing frequency counts for different categories across countries,
    and include the total sum for each category.
    
    Parameters:
    - all_dfs (dict): A nested dictionary where the first level keys are countries and second level keys
                      are categories, with DataFrames as values.
    
    Returns:
    - dict: A dictionary with countries as keys and DataFrames of frequency counts as values,
            including a 'Total' row showing the sum for each category.
    """
    
    dict_country = {}
    for country, dfs in all_dfs.items():
        freq_tables = {}
        for category, df in dfs.items():
            freq_counts = df['frequency_y'].value_counts().sort_index()
            freq_tables[category] = freq_counts
        
        # Combine all frequency count Series into a single DataFrame
        freq_df = pd.concat(freq_tables, axis=1)
        freq_df = freq_df.fillna(0).astype(int)  # Fill missing values with 0
        # Order the df by frequency - smallest to largest
        freq_df = freq_df.sort_index()
        # Add a row with the sum total for each column
        freq_df.loc['Total'] = freq_df.sum()
        # divide each row by the total number of TEs
        freq_df = freq_df.div(freq_df.loc['Total'], axis=1)
        # Remove the 'Total' row from the frequency counts
        freq_df = freq_df.drop('Total')
        
        dict_country[country] = freq_df
    
    return dict_country
            
def plot_SFS(all_dfs, name="SFS"):
    """
    Plot the Site Frequency Spectrum (SFS) for each country with all categories displayed
    side by side as different colored bars.

    Parameters:
    - all_dfs: Dictionary where keys are country names and values are DataFrames
               with frequencies as index and categories as columns.
    - name: Base name for the output SVG files.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    
    colors = {
        'overlaping_rep_exons': 'blue',
        'overlapping_exons_not_rep': 'green',
        'overlapping_rep_not_exon': 'orange',
        'not_overlapping_exon_rep': 'red'
    }

    for country, df in all_dfs.items():
        # Filter out the 'Total' row if it exists
        if 'Total' in df.index:
            plot_df = df.drop('Total')
        else:
            plot_df = df
            
        # Get categories present in this dataset
        categories = [cat for cat in colors.keys() if cat in plot_df.columns]
        
        # Set up the plot
        plt.figure(figsize=(12, 6))
        
        # Calculate bar width based on number of categories
        bar_width = 0.8 / len(categories)
        
        # Set up x positions for each group of bars
        x = np.arange(len(plot_df.index))
        
        # Plot each category as a set of bars
        for i, category in enumerate(categories):
            # Calculate the position for this category's bars
            x_pos = x - (0.4 - bar_width/2) + i * bar_width
            
            plt.bar(x_pos, 
                   plot_df[category].values, 
                   width=bar_width, 
                   color=colors[category], 
                   label=category.replace("_", " ").title())
        
        # Set up plot labels and legend
        plt.xlabel('Frequency')
        plt.ylabel('Count of insertions')
        plt.title(f'Site Frequency Spectrum - {country}')
        
        # Set x-tick positions at the center of each group
        plt.xticks(x, [f"{val:.2f}" for val in plot_df.index], rotation=45)
        
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"{name}_{country}_all_categories.svg")
        plt.close()
        
        # Also create individual plots for each category (optional)
        for category in categories:
            plt.figure(figsize=(10, 5))
            plt.bar(x, plot_df[category].values, width=0.6, color=colors[category])
            plt.xlabel('Frequency')
            plt.ylabel('Count of insertions')
            plt.title(f'{category.replace("_", " ").title()} - {country}')
            plt.xticks(x, [f"{val:.2f}" for val in plot_df.index], rotation=45)
            plt.tight_layout()
            plt.savefig(f"{name}_{country}_{category}.svg")
            plt.close()

repetitive_overlap_ref_files = '/work/users/g/a/gabivla/lab/SV_mosquitos/overlap_analysis_repetitive_regions/'
files = os.listdir(repetitive_overlap_ref_files)
# keep only files from list that have _ref on them
files = [file for file in files if '_ref' in file]
# make df from the first file in the list
df_overlap_repetitive = pd.read_csv(repetitive_overlap_ref_files + files[0])

exons_overlap_nonref_files = '/work/users/g/a/gabivla/lab/SV_mosquitos/overlap_analysis/'
files = os.listdir(exons_overlap_nonref_files)
# keep only files from list that have _ref on them
files = [file for file in files if '_ref' in file]
# make df from the first file in the list
df_overlap = pd.read_csv(exons_overlap_nonref_files + files[0])

# get VCF to check frequency in each country
df = pd.read_csv('/work/users/g/a/gabivla/lab/SV_mosquitos/VCF_genotypes_TEs/separate_chromosomes_windows_reference/merged_reference/VCF_reference_TEs_post_processing_600_absence_evidence.csv')

# fixing df
#df = df.drop(columns=['Unnamed: 0'])

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
df = df.rename(columns=rename_dict)

# Extract unique countries from the renamed genome columns only
countries = set()
for original_col in genome_cols:
    new_col = rename_dict[original_col]
    country = new_col.split('_')[-1]
    countries.add(country)

country_vcf_dict_ref = {}
country_pairs_nonref = list(combinations(countries, 2))

all_dfs_ref = {}
percentage_overlapping = {}
all_dfs_ref_rep = {}
percentage_overlapping_rep = {}

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
    # calculate frequencies
    frequency_df = calculate_frequency_TE(country_df, country)

    # remove rows were frequency == 0
    frequency_df = frequency_df[frequency_df['frequency'] > 0]

    merged_frequency_df = frequency_df.merge(
    df_overlap[['TE_name', 'chromosome', 'position_start', 'position_end', 'attributes']],
    on=['TE_name', 'chromosome', 'position_start', 'position_end'],
    how='left'  # Keeps all rows in frequency_df, adds matching attributes if found
    )

    #if 'attributes' not NaN:
    matched_rows = merged_frequency_df[merged_frequency_df['attributes'].notna()]

    # get rows that have NaN in 'attributes'
    non_matched_rows = merged_frequency_df[merged_frequency_df['attributes'].isna()]
    
    # number of TEs overlapping genes
    TE_overlapping = len(matched_rows)
    total_number_TEs = len(frequency_df)

    percentage_overlapping[country] = {'number_te_overlapping': TE_overlapping, 'total_number_TEs': total_number_TEs, 'percentage': TE_overlapping/total_number_TEs}
    all_dfs_ref[country] = matched_rows, non_matched_rows
    # checking repetitive regions
    merged_frequency_df_repetitive = frequency_df.merge(
    df_overlap_repetitive[['TE_name', 'chromosome', 'position_start', 'position_end', 'rep_start']],
    on=['TE_name', 'chromosome', 'position_start', 'position_end'],
    how='left'  # Keeps all rows in frequency_df, adds matching attributes if found
    )
    #if 'attributes' not NaN:
    matched_rows_rep = merged_frequency_df_repetitive[merged_frequency_df_repetitive['rep_start'].notna()]
    # unique set of TEs overlapping repeats
    matched_rows_rep = matched_rows_rep.drop_duplicates(
        subset=['TE_name', 'chromosome', 'position_start', 'position_end']
    )
    print(f"Unique TEs overlapping repetitive in {country}: {len(matched_rows_rep)}")

    # get rows that have NaN in 'attributes'
    non_matched_rows_rep = merged_frequency_df_repetitive[merged_frequency_df_repetitive['rep_start'].isna()]
    
    # number of TEs overlapping genes
    TE_overlapping_rep = len(matched_rows_rep)
    total_number_TEs_rep = len(merged_frequency_df_repetitive)

    percentage_overlapping_rep[country] = {'number_te_overlapping': TE_overlapping_rep, 'total_number_TEs': total_number_TEs_rep, 'percentage': TE_overlapping_rep/total_number_TEs_rep}

    # Save the merged DataFrame for this country
    all_dfs_ref_rep[country] = matched_rows_rep, non_matched_rows_rep

# Plot the results
# Extract values for plotting
countries = list(percentage_overlapping.keys())
percentages = [percentage_overlapping[country]['percentage'] * 100 for country in countries]  # Convert to percentage
percentages_ref = [percentage_overlapping_rep[country]['percentage'] * 100 for country in countries]  # Convert to percentage

top_genes = {}  # Store the most frequent gene for each country
gene_frequencies = {}  # Store the frequency of the top gene in the population
all_categories = {}

for country, dfs in all_dfs_ref.items():
    df = dfs[0]
    not_overlapping = dfs[1]
    df_rep, not_overlapping_rep = all_dfs_ref_rep[country]

    # get values of TEs overlapping exons and repetitive regions
    overlaping_rep_exons = df.merge(
    df_rep,
    on=["position_start", "position_end", "TE_name"],
    how="inner"
    )
    # get values of TEs overlapping exons and not overlapping repetitive regions
    overlapping_exons_not_rep = df.merge(
    not_overlapping_rep,
    on=["position_start", "position_end", "TE_name"],
    how="inner"
    )
    # get values of TEs overlapping repetitive regions and not overlapping exons
    overlapping_rep_not_exon = df_rep.merge(
    not_overlapping,
    on=["position_start", "position_end", "TE_name"],
    how="inner"
    )
    # get values of TEs not overlapping exons or repetitive regions
    not_overlapping_exon_rep = not_overlapping.merge(
    not_overlapping_rep, 
    on=["position_start", "position_end", "TE_name"], 
    how="inner"
    )

    all_categories[country] = {
        'overlaping_rep_exons': overlaping_rep_exons,
        'overlapping_exons_not_rep': overlapping_exons_not_rep,
        'overlapping_rep_not_exon': overlapping_rep_not_exon,
        'not_overlapping_exon_rep': not_overlapping_exon_rep
    }

categories_examples = SFS_data(all_categories)
#print(categories_examples)
plot_SFS(categories_examples)