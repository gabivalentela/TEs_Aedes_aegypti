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

overlap_nonref_files = '/work/users/g/a/gabivla/lab/SV_mosquitos/overlap_analysis/'
files = os.listdir(overlap_nonref_files)
# keep only files from list that have _ref on them
files = [file for file in files if '_nonref' in file]
# make df from the first file in the list
df_overlap = pd.read_csv(overlap_nonref_files + files[0])
#print(df_overlap)
# groupby TE_name, chromosome, position_start, position_end
# and get the first row of each group
df_overlap = df_overlap.groupby(['TE_name', 'chromosome', 'position_start', 'position_end']).first().reset_index()
#print(df_overlap)


'''attribute_columns = df_overlap['attributes']
# split this column with ';'
attribute_columns = attribute_columns.str.split(';', expand=True)
# count how many time each combination of geneid and exon id appears - column with index 2 and 0
count = attribute_columns.groupby([2, 0]).size()'''

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
    
    # Save the merged DataFrame for this country
    all_dfs_ref[country] = matched_rows, non_matched_rows

# Plot the results
# Extract values for plotting
countries = list(percentage_overlapping.keys())
percentages = [percentage_overlapping[country]['percentage'] * 100 for country in countries]  # Convert to percentage

# Create the bar plot
plt.figure(figsize=(10, 8))
plt.bar(countries, percentages, color='m', alpha=0.7, edgecolor='black')

# Add labels above bars
for i, perc in enumerate(percentages):
    plt.text(i, perc + 1, f"{perc:.1f}%", ha='center', fontsize=10, fontweight='bold')

# Labels and title
plt.xlabel("Country", fontsize=15)
plt.ylabel("Percentage of TEs Overlapping Genes (%)", fontsize=15)
plt.title("Percentage of Transposable Elements (TEs) Overlapping Exons by Country")
plt.xticks(rotation=45, fontsize=15)
plt.yticks(fontsize=12)
plt.ylim(0, max(percentages) + 5)  # Adjust y-axis limit for better visibility

# Show the plot
plt.savefig('percentage_TEs_overlapping_genes_new.svg')

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

top_genes = {}  # Store list of top genes for each country
gene_frequencies = {}  # Store frequency of top genes (averaged)
TE_frequencies_no_overlap = {}
TE_frequencies_overlap = {}

for country, dfs in all_dfs_ref.items():
    df = dfs[0]
    not_overlapping = dfs[1]

    # Split 'attributes' column and extract the gene ID (assuming it's in index 2)
    attribute_columns = df['attributes'].str.split(';', expand=True)

    # Count occurrences of each gene ID
    gene_counts = attribute_columns[2].value_counts()
    max_count = gene_counts.max()

    # Get all genes with the maximum count
    top_gene_ids = gene_counts[gene_counts == max_count].index.tolist()

    # Store gene names and count
    top_genes[country] = [(gene_id, max_count) for gene_id in top_gene_ids]

    # Get the average frequency of these top genes
    freqs = []
    for gene_id in top_gene_ids:
        mask = df['attributes'].str.contains(gene_id, na=False)
        freqs.append(df[mask]['frequency'].mean())
    gene_frequencies[country] = np.mean(freqs)

    # Frequencies of overlapping and not overlapping TEs
    gene_freq_overlapping = df['frequency'].mean()
    gene_freq_not_overlapping = not_overlapping['frequency'].mean()

    TE_frequencies_overlap[country] = gene_freq_overlapping
    TE_frequencies_no_overlap[country] = gene_freq_not_overlapping

# === Plot 1: Mean TE Frequency (Overlapping vs Not Overlapping Genes) ===
df_plot = pd.DataFrame({
    'Overlapping': TE_frequencies_overlap,
    'Not Overlapping': TE_frequencies_no_overlap
})
print(df_plot)
# get the average frequency of TEs overlapping genes and not overlapping genes
df_plot = df_plot.mean().reset_index()
df_plot.columns = ['TE Type', 'Mean Frequency']
print(df_plot)

'''df_plot.plot(kind='bar', figsize=(10, 7), color=['red', 'lightpink'], edgecolor='black')
plt.xlabel("Country", fontsize=15)
plt.xticks(rotation=45, fontsize=15)
plt.yticks(fontsize=15)
plt.ylabel("Mean TE Frequency", fontsize=15)
plt.title("Mean TE Frequency: Overlapping vs. Not Overlapping Exons")
plt.xticks(rotation=45)
plt.legend(title="TE Type")
plt.tight_layout()
plt.savefig('TE_frequency_overlapping_vs_not_new.svg')

# === Plot 2: All Top Genes (multiple genes per country) ===
plot_data = []
for country, gene_list in top_genes.items():
    for gene_id, count in gene_list:
        gene_clean = gene_id.replace('gene_id=', '').strip()
        plot_data.append({'Country': country, 'Gene': gene_clean, 'Count': count})

top_genes_df = pd.DataFrame(plot_data)

plt.figure(figsize=(12, 6))
sns.barplot(data=top_genes_df, x='Country', y='Count', hue='Gene', dodge=True)
plt.xlabel("Country")
plt.ylabel("Number of TEs overlapping gene exons")
plt.title("Genes with the Highest Number of Insertions per Country")
plt.xticks(rotation=45)
plt.legend(title="Gene", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig('all_top_genes_per_country.svg')

# === Plot 3: Combined Bar (Count) + Line (Frequency) ===
fig, ax1 = plt.subplots(figsize=(12, 6))

# Bar plot: count (still shows only one bar per country)
aggregated_counts = {country: gene_list[0][1] for country, gene_list in top_genes.items()}
print(aggregated_counts)
countries = list(aggregated_counts.keys())
counts = list(aggregated_counts.values())
ax1.bar(countries, counts, color='b', alpha=0.6)

# Annotate top gene names (comma-separated)
for i, country in enumerate(countries):
    gene_names = ', '.join([g[0].replace('gene_id=', '').strip() for g in top_genes[country]])
    ax1.text(i, counts[i] + 1, gene_names, ha="center", fontsize=9, rotation=30)

ax1.set_xlabel("Country")
ax1.set_ylabel("Number of TEs overlapping gene exons", color="b")
ax1.tick_params(axis="y", labelcolor="b")

# Line plot: average frequency of top genes
print('Average frequencies of top genes:', gene_frequencies)
ax2 = ax1.twinx()
ax2.plot(countries, list(gene_frequencies.values()), marker="o", color="r", label="Avg freq of top genes")
ax2.set_ylabel("Average TE Frequency", color="r")
ax2.tick_params(axis="y", labelcolor="r")

# Legends
#ax1.legend(["Insertion Count"], loc="upper left")
ax2.legend(loc="upper right")

plt.xticks(rotation=45)
fig.tight_layout()
plt.savefig('combined_top_gene_freq_and_counts_new.svg')'''