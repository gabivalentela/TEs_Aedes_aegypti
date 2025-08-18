import os
import pandas as pd
from itertools import combinations
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

repetitive_overlap_ref_files = '/work/users/g/a/gabivla/lab/SV_mosquitos/overlap_analysis_repetitive_regions/'
files = os.listdir(repetitive_overlap_ref_files)
# keep only files from list that have _ref on them
files = [file for file in files if '_nonref' in file]
# make df from the first file in the list
df_overlap_repetitive = pd.read_csv(repetitive_overlap_ref_files + files[0])

exons_overlap_nonref_files = '/work/users/g/a/gabivla/lab/SV_mosquitos/overlap_analysis/'
files = os.listdir(exons_overlap_nonref_files)
# keep only files from list that have _ref on them
files = [file for file in files if '_nonref' in file]
# make df from the first file in the list
df_overlap = pd.read_csv(exons_overlap_nonref_files + files[0])
df_overlap = df_overlap.groupby(['TE_name', 'chromosome', 'position_start', 'position_end']).first().reset_index()

# get VCF to check frequency in each country
df = pd.read_csv('/work/users/g/a/gabivla/lab/SV_mosquitos/VCF_genotypes_TEs/separate_chromosomes_windows_reference/merged_nonreference/VCF_nonreference_post_processing.csv')

# fixing df
df = df.drop(columns=['Unnamed: 0'])

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

# Define the desired order of countries
desired_order = ['Brazil', 'Colombia', 'USA', 'Gabon', 'Kenya', 'Senegal']

# Filter to only include countries that exist in your data and maintain the desired order
countries_ordered = [country for country in desired_order if country in countries]

# Split data by country
for i, country in enumerate(countries_ordered):

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

    # get rows that have NaN in 'attributes'
    non_matched_rows_rep = merged_frequency_df_repetitive[merged_frequency_df_repetitive['rep_start'].isna()]
    
    # number of TEs overlapping genes
    TE_overlapping_rep = len(matched_rows_rep)
    total_number_TEs_rep = len(merged_frequency_df_repetitive)

    percentage_overlapping_rep[country] = {'number_te_overlapping': TE_overlapping_rep, 'total_number_TEs': total_number_TEs_rep, 'percentage': TE_overlapping_rep/total_number_TEs_rep}

    # Save the merged DataFrame for this country
    all_dfs_ref_rep[country] = matched_rows_rep, non_matched_rows_rep

# Plot the results
countries_plot = countries_ordered  # Use the ordered list
percentages = [percentage_overlapping[country]['percentage'] * 100 for country in countries_plot]
percentages_ref = [percentage_overlapping_rep[country]['percentage'] * 100 for country in countries_plot]

# Create the bar plot - exons
plt.figure(figsize=(10, 6))
plt.bar(countries_plot, percentages, color='m', alpha=0.7, edgecolor='black')

# Add labels above bars
for i, perc in enumerate(percentages):
    plt.text(i, perc + 1, f"{perc:.1f}%", ha='center', fontsize=10, fontweight='bold')

# Labels and title
plt.xlabel("Country")
plt.ylabel("Percentage of TEs Overlapping Genes (%)")
plt.title("Percentage of Transposable Elements (TEs) Overlapping Genes by Country")
plt.xticks(rotation=45)
plt.ylim(0, max(percentages) + 5)

plt.savefig('percentage_TEs_overlapping_genes.svg')

# Create the bar plot - repetitive
plt.figure(figsize=(10, 6))
plt.bar(countries_plot, percentages_ref, color='m', alpha=0.7, edgecolor='black')

# Add labels above bars
for i, perc in enumerate(percentages_ref):
    plt.text(i, perc + 1, f"{perc:.1f}%", ha='center', fontsize=10, fontweight='bold')

# Labels and title
plt.xlabel("Country")
plt.ylabel("Percentage of TEs Overlapping repetitive (%)")
plt.title("Percentage of Transposable Elements (TEs) Overlapping repetitive by Country")
plt.xticks(rotation=45)
plt.ylim(0, max(percentages_ref) + 5)

plt.savefig('percentage_TEs_overlapping_rep.svg')

top_genes = {}  # Store the most frequent gene for each country
gene_frequencies = {}  # Store the frequency of the top gene in the population
TE_frequencies_exon_rep = {}
TE_frequencies_exon_not_rep = {}
TE_frequencies_not_exon_rep = {}
TE_frequencies_not_exon_not_rep = {}

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

    attribute_columns = df['attributes'].str.split(';', expand=True)
    
    gene_counts = attribute_columns[2].value_counts()
    
    # Get the most common gene
    top_gene = gene_counts.idxmax()
    top_gene_count = gene_counts.max()

    # Get the frequency of this gene in the population
    gene_freq = df[df['attributes'].str.contains(top_gene, na=False)]['frequency'].mean()
    # mean frequencies
    # exon and rep
    gene_freq_exon_rep = overlaping_rep_exons['frequency_x'].mean()
    # exon and not rep
    gene_freq_exon_not_rep = overlapping_exons_not_rep['frequency_x'].mean()
    # not exon and rep
    gene_freq_not_exon_rep = overlapping_rep_not_exon['frequency_x'].mean()
    # not exon and not rep
    gene_freq_not_exon_not_rep = not_overlapping_exon_rep['frequency_x'].mean()

    # Store values
    top_genes[country] = (top_gene, top_gene_count)
    gene_frequencies[country] = gene_freq
    TE_frequencies_exon_rep[country] = gene_freq_exon_rep
    TE_frequencies_exon_not_rep[country] = gene_freq_exon_not_rep
    TE_frequencies_not_exon_rep[country] = gene_freq_not_exon_rep
    TE_frequencies_not_exon_not_rep[country] = gene_freq_not_exon_not_rep

df_plot = pd.DataFrame({
    'Exonic and Repetitive': TE_frequencies_exon_rep,
    'Exonic and Non-Repetitive': TE_frequencies_exon_not_rep,
    'Non-Exonic and Repetitive': TE_frequencies_not_exon_rep,
    'Non-Exonic and Non-Repetitive': TE_frequencies_not_exon_not_rep
})

print(df_plot)

ax = df_plot.plot(
    kind='bar',
    figsize=(10, 6),
    color=["#fcbba1", "#fc9272", "#fb6a4a", "#cb181d"],
    edgecolor='black'
)

# Place legend outside the plot with a title
ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), title="TE Type")

# Adjust layout to avoid cutting off labels/legend
plt.tight_layout()
plt.xlabel("Country")
plt.ylabel("Mean TE Frequency")
plt.title("Mean TE Frequency: Overlapping vs. Not Overlapping Genes")
plt.xticks(rotation=45)

plt.savefig('TE_frequency_overlapping_vs_not.svg', bbox_inches='tight')

# Plot the results
fig, ax1 = plt.subplots(figsize=(12, 6))

# Bar plot for gene counts
countries = list(top_genes.keys())
gene_counts = [v[1] for v in top_genes.values()]
ax1.bar(countries, gene_counts, color='b', alpha=0.6)

# Add gene name annotations above bars
for i, country in enumerate(countries):
    gene_name = top_genes[country][0]  # Get the gene name
    gene_name = gene_name.replace('gene_id=', '')  # Remove 'gene_id=' prefix
    ax1.text(i, gene_counts[i] + 2, gene_name, ha="center", fontsize=10, rotation=30, color="black")

ax1.set_xlabel("Country")
ax1.set_ylabel("Number of TEs overlapping gene exons", color="b")
ax1.tick_params(axis="y", labelcolor="b")

# Create a second y-axis
ax2 = ax1.twinx()

# Line plot for TEs overlapping genes
ax2.plot(countries, gene_frequencies.values(), marker="o", color="r", label="TEs overlapping this GENE")

# Line plot for TEs overlapping repetitive elements
#ax2.plot(countries, TE_frequencies_overlap_rep.values(), marker="^", linestyle="--", color="green", label="TEs overlapping repetitive")

# Line plot for TEs not overlapping genes
# ax2.plot(countries, TE_frequencies_no_overlap.values(), marker="s", color="purple", label="TEs not overlapping genes")
# ax2.plot(countries, TE_frequencies_overlap.values(), marker="s", color="orange", label="All TEs overlapping genes")

ax2.set_ylabel("Average TE Frequency", color="r")
ax2.tick_params(axis="y", labelcolor="r")

# Add legends
ax1.legend(loc="upper left")
ax2.legend(loc="upper right")

plt.xticks(rotation=45)
fig.tight_layout()
plt.savefig('genes_overlapping_TEs_and_all_tes.svg')
