import os
import sys
import pandas as pd
import polars as pl
import matplotlib.pyplot as plt
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

def count_TE(general_df):
    # Ensure TE names are consistently formatted
    general_df = general_df.with_columns(pl.col("TE_name").str.to_uppercase())

    possible_classes = ['LINE', 'SINE', 'LTR', 'MITE', 'TIR', 'UNKNOWN', 'HELITRON', 'PENELOPE']

    # Extract clean TE names - BASICALLY CHECK WHICH OF THE POSSIBLE CLASSES IS PRESENT IN THE TE NAME AND KEEP ONLY THAT - ADD OTHER IF NOT IN MAIN CLASSES
    df_te = general_df.with_columns(
        pl.col("TE_name").map_elements(
            lambda x: next((cls for cls in possible_classes if cls in x.upper()), 'OTHER'),
            return_dtype=pl.Utf8
        ).alias("TE_name")
    )
    # Count occurrences of each TE name
    grouped = df_te.group_by('TE_name').agg(pl.len().alias("count"))
    print(f'this is the grouped df: {grouped}')
    # Add percentage column
    total_count = grouped["count"].sum()
    grouped = grouped.with_columns((pl.col("count") / total_count * 100).alias("percentage"))

    # Sort by percentage in descending order
    sorted_grouped = grouped.sort("percentage", descending=True)
    return sorted_grouped

def remove_TEs_not_present(country_df):
    # Identify genome columns
    genome_cols = [col for col in country_df.columns if col not in ['chromosome', 'position_start', 'position_end', 'TE_name', 'strand']]
    
    mask = country_df.select(genome_cols).map_rows(lambda row: '1/1' in row or '0/1' in row).to_series()

    filtered_df = country_df.filter(mask)
    return filtered_df

def count_te_per_chr(country, df, type_data):
    # Ensure TE names are consistently formatted
    df = df.with_columns(pl.col("TE_name").str.to_uppercase())
    
    possible_classes = ['LINE', 'SINE', 'LTR', 'TIR', 'MITE', 'UNKNOWN', 'HELITRON', 'PENELOPE']

    # Extract clean TE names - BASICALLY CHECK WHICH OF THE POSSIBLE CLASSES IS PRESENT IN THE TE NAME AND KEEP ONLY THAT
    df = df.with_columns(
        pl.col("TE_name").map_elements(
            lambda x: next((cls for cls in possible_classes if cls in x.upper()), 'OTHER'),
            return_dtype=pl.Utf8
        ).alias("TE_name")
    )

    # Filter chromosomes starting with 'AaegL'
    df = df.filter(pl.col("chromosome").str.starts_with("AaegL"))

    # Count occurrences of each TE type per chromosome
    grouped = df.group_by(["chromosome", "TE_name"]).agg(pl.len().alias("count"))

    # Add percentage column within each chromosome
    grouped = grouped.with_columns(
    grouped.join(
        grouped.group_by("chromosome").agg(pl.col("count").sum().alias("total_count")),
            on="chromosome"
        ).select((pl.col("count") / pl.col("total_count") * 100).alias("percentage"))
    )

    # Sort results by chromosome and percentage
    grouped = grouped.sort(["chromosome", "percentage"], descending=[False, True])

    print(f"In this country, TE percentages per chromosome are as follows (filtered by AaegL): {country}\n{grouped}")

    # Plot the TE percentages per chromosome
    chromosomes = grouped['chromosome'].to_list()
    te_names = grouped['TE_name'].to_list()
    percentages = grouped['count'].to_list()

    # Organize data by chromosome
    unique_chromosomes = sorted(set(chromosomes))
    te_counts = {chrom: {} for chrom in unique_chromosomes}

    for chrom, te, pct in zip(chromosomes, te_names, percentages):
        if te not in te_counts[chrom]:
            te_counts[chrom][te] = 0
        te_counts[chrom][te] += pct

    # Prepare data for the stacked bar plot
    stacked_data = {}
    for chrom, te_dict in te_counts.items():
        for te, count in te_dict.items():
            if te not in stacked_data:
                stacked_data[te] = []
            stacked_data[te].append(count)

    # Fill missing values with 0 for chromosomes with no data for a specific TE
    for te, counts in stacked_data.items():
        for i in range(len(counts), len(unique_chromosomes)):
            counts.append(0)

    # Create the stacked bar chart
    fig, ax = plt.subplots(figsize=(10, 6))
    bottom_values = np.zeros(len(unique_chromosomes))

    # Generate a color map using the 'Blues' palette
    color_map = plt.get_cmap('turbo', len(stacked_data))
    for i, (te, counts) in enumerate(stacked_data.items()):
        ax.bar(
            unique_chromosomes, 
            counts, 
            bottom=bottom_values, 
            label=te, 
        color=color_map(i / len(stacked_data))
        )
        bottom_values += np.array(counts)

    ax.set_xlabel('Chromosome')
    ax.set_ylabel('Counts')
    ax.set_title(f'TE Counts per Chromosome in {country}')
    ax.legend(title='TE Name', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(rotation=90)
    plt.tight_layout()

    # Save the plot to a file
    plot_filename = f"{type_data}_TE_counts_per_chromosome_{country}.png"
    plt.savefig(plot_filename, bbox_inches='tight')
    plt.close() 

    return grouped
    
def plot_TE_counts(df, key):
    # Fixed country order
    country_order = ["Brazil", "USA", "Colombia", "Senegal", "Gabon", "Kenya"]

    # Calculate total percentage for each TE across countries for sorting
    te_totals = (
        df.group_by('TE_name')
        .agg(pl.col('percentage').sum())
        .sort('percentage', descending=True)
    )

    # Get ordered unique TEs
    unique_tes = te_totals['TE_name'].to_list()
    num_countries = len(country_order)

    color_mapping = {
        "Brazil":   "#60D08F",  # Teal Green
        "USA":      "#785EF0",  # Strong Red
        "Colombia": "#DC267F",  # Orange
        "Senegal":  "#FE6100",  # Blue
        "Gabon":    "#FFB000",  # Vermillion
        "Kenya":    "#648FFF"   # Reddish Purple
    }

    # Plotting logic
    plt.figure(figsize=(12, 8))

    bar_height = 0.8
    positions = np.arange(len(unique_tes))
    width = bar_height / num_countries

    for idx, country in enumerate(country_order):
        country_data = df.filter(pl.col('country') == country)
        te_percentages = dict(zip(country_data['TE_name'].to_list(), country_data['percentage'].to_list()))
        percentages = [te_percentages.get(te, 0) for te in unique_tes]
        offset = positions + (idx - num_countries / 2 + 0.5) * width
        color = color_mapping.get(country, "#7f7f7f")
        plt.barh(offset, percentages, height=width, label=country, color=color)

    plt.yticks(positions, unique_tes, fontsize=30)
    plt.xticks(fontsize=30)
    plt.xlabel('Percentage', fontsize=30)
    plt.title('Top TEs by Percentage per Country', fontsize=22)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=30)
    plt.tight_layout()
    # Save the plot to a file
    plot_filename = f'top_TEs_by_percentage_{key}.svg'
    plt.savefig(plot_filename, bbox_inches='tight')
    plt.close()

def plot_PCA_data(df, base_columns, key):
    # getting genotype columns - all the ones not in base_columns
    genotype_columns = [col for col in df.columns if col not in base_columns]
    # making genotyping df
    genotype_df = df[genotype_columns]
    # print number of columns 
    print(f'Number of genotype columns: {len(genotype_columns)}')

    tie_count = [0]  # mutable counter, captured by closure
    # set seed for reproducibility
    np.random.seed(42)
    # Determine the most frequent genotype for each TE
    def fill_missing(row):
        most_common = row[row != './.'].mode()
        # if there is more than one most common value randomly select one of them
        if len(most_common) > 1:
            most_common = np.random.choice(most_common)
            tie_count[0] += 1
        # if only one most common value, just take that
        else:
            most_common = most_common[0]
        return row.replace({'./.': most_common})

    # Apply the function to fill missing values
    genotype_df = genotype_df.apply(fill_missing, axis=1)

    print(f'Number of TEs with tied genotype counts: {tie_count[0]} out of {len(genotype_df)}')

    genotype_df = genotype_df.T
    # Convert genotype strings to numerical values
    genotype_df = genotype_df.replace({'0/0': 0, '0/1': 1, '1/1': 2})
    genotype_df = genotype_df.infer_objects().astype(int)
    # Extract country names from the row index
    countries = [row.split('_')[-1] for row in genotype_df.index]
    # Extract unique countries
    unique_countries = sorted(set(countries))  
    # Define a color mapping for each country
    color_mapping = {
        "Brazil":   "#60D08F",  # Teal Green
        "USA":      "#785EF0",  # Strong Red
        "Colombia": "#DC267F",  # Orange
        "Senegal":  "#FE6100",  # Blue
        "Gabon":    "#FFB000",  # Vermillion
        "Kenya":    "#648FFF"   # Reddish Purple
    }
    # Standardize the data
    scaler = StandardScaler()
    standardized_data = scaler.fit_transform(genotype_df)
    # Apply PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(standardized_data)
    # Create a scatter plot
    fig, ax = plt.subplots(figsize=(15, 8))
    # Plot each sample, color-coded by country
    for country in unique_countries:
        country_mask = np.array(countries) == country
        country_data = pca_result[country_mask]
        ax.scatter(country_data[:, 0], country_data[:, 1],
                   alpha=0.7,
                   label=country,
                   color=color_mapping.get(country, "#666666"))
    # Labels and title
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.2f}% Variance)", fontsize=30)
    # make font bigger for the x and y ticks
    ax.tick_params(axis='both', which='major', labelsize=30)
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.2f}% Variance)", fontsize=30)
    ax.set_title("PCA of TE Insertions", fontsize=30)
    # Put legend outside the figure on the right
    fig.legend(
    title="Country",
    loc='center left',
    bbox_to_anchor=(1.01, 0.5),
    fontsize=35,
    title_fontsize=35,
    markerscale=3,       # increases the size of marker symbols
    handlelength=3,      # increases the length of line handles (for lines)
    borderpad=1.5        # optional: space around the legend
)
    # Adjust layout to make room for the legend
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    # Save the figure
    plt.savefig(f'PCA_{key}_TEs.svg', bbox_inches='tight')
    plt.close()

def proportion_of_TEs(df, key):
    # calculating TE sizes
    df['TE_size'] = df['position_end'] - df['position_start']

    chromosome_1 = 310827022
    chromosome_2 = 474425716
    chromosome_3 = 409777670
    total_genome_size = chromosome_1 + chromosome_2 + chromosome_3

    total_te_size = df['TE_size'].sum()
    proportion = total_te_size / total_genome_size * 100
    print(f"Total TE genome proportion: {proportion:.2f}%")

    possible_classes = ['LINE', 'SINE', 'LTR', 'MITE', 'TIR', 'UNKNOWN', 'HELITRON', 'PENELOPE', 'OTHER']

    # Convert to Polars
    general_df = pl.DataFrame(df).with_columns(
        pl.col("TE_name").str.to_uppercase()
    )

    # Assign TE class
    df_te = general_df.with_columns(
        pl.col("TE_name").map_elements(
            lambda x: next((cls for cls in possible_classes if cls in x), x),
            return_dtype=pl.Utf8
        ).alias("TE_name")
    )

    grouped = df_te.group_by('TE_name').agg(
        pl.col("TE_size").sum().alias("total_bp")
    )

    # Calculate percentage of genome
    grouped = grouped.with_columns(
        (pl.col("total_bp") / total_genome_size * 100).alias("genome_percentage")
    )

    # Sort
    sorted_grouped = grouped.sort("genome_percentage", descending=True)

    # plotting the proportion of the genome taken up by each TE type
    plt.figure(figsize=(12, 8))
    plt.bar(sorted_grouped['TE_name'], sorted_grouped['genome_percentage'], color='skyblue')
    plt.xlabel('TE Type', fontsize=30)
    plt.ylabel('Percentage of Genome', fontsize=30)
    plt.title('Proportion of Genome Taken Up by Each TE Type', fontsize=30)
    plt.xticks(rotation=45, fontsize=30)
    plt.yticks(fontsize=30)
    plt.tight_layout()
    plt.savefig(f'proportion_of_genome_by_TE_type_{key}.svg', bbox_inches='tight')
    plt.close()

    # drop the TE_size column from the original df
    df = df.drop(columns=['TE_size'])
    return df

# Running functions
reference_path = '/work/users/g/a/gabivla/lab/SV_mosquitos/TE_paper_aedesaegytpi/reference_vcf/genotyping_reference_TEs/post_processing/merged_results_5_2_either_filtered.csv'
non_reference_path = '/work/users/g/a/gabivla/lab/SV_mosquitos/TE_paper_aedesaegytpi/nonreference_vcf/post_processing_nonreference_TEs/post_processing_overlap_df/max_min_updated_strategy/final_df_overlap_retroseq_teforest_nonreference_1000.csv'

# Read only necessary columns from location file
location_path = '/nas/longleaf/home/gabivla/lab/SV_mosquitos/location_information/mcclintock_output_locations.csv'
location_df = pd.read_csv(location_path, usecols=['Genome_code', 'country'])
location_dict = dict(zip(location_df['Genome_code'], location_df['country']))

# making the reference df
df_reference = pd.read_csv(reference_path)

# Define base columns that should not be modified
base_cols_ref = ['chromosome', 'position_start', 'position_end', 'TE_name', 'strand', 'score']

# Get genome columns at once
genome_cols = [col for col in df_reference.columns if col not in base_cols_ref]
# update the genome columns names
for genome in genome_cols:
    if genome.startswith('SRR'):
        genome_updated = genome.replace('_1', '')
    elif genome.startswith('FEMALE') or genome.startswith('MALE'):
        genome_updated = genome.replace('_L001_L002_R1_001', '')
    # update the column name in the dataframe
    df_reference = df_reference.rename(columns={genome: genome_updated})

genome_cols = [col for col in df_reference.columns if col not in base_cols_ref]

# Create rename dictionary and apply it once
rename_dict = {col: f"{col}_{location_dict.get(col, '')}" for col in genome_cols}
df_reference = df_reference.rename(columns=rename_dict)

# Extract unique countries from the renamed genome columns only
countries = set()
for original_col in genome_cols:
    new_col = rename_dict[original_col]
    country = new_col.split('_')[-1]
    countries.add(country)

# making the nonreference df
df_non_reference = pd.read_csv(non_reference_path)

# Define base columns that should not be modified
base_cols_nonref = ['chromosome', 'position_start', 'position_end', 'TE_name', 'strand']

# Get genome columns at once
genome_cols = [col for col in df_non_reference.columns if col not in base_cols_nonref]
# update the genome columns names
for genome in genome_cols:
    if genome.startswith('SRR'):
        genome_updated = genome.replace('_1', '')
    elif genome.startswith('FEMALE') or genome.startswith('MALE'):
        genome_updated = genome.replace('_L001_L002_R1_001', '')
    # update the column name in the dataframe
    df_non_reference = df_non_reference.rename(columns={genome: genome_updated})

genome_cols = [col for col in df_non_reference.columns if col not in base_cols_nonref]

# Create rename dictionary and apply it once
rename_dict = {col: f"{col}_{location_dict.get(col, '')}" for col in genome_cols}
df_non_reference = df_non_reference.rename(columns=rename_dict)

# Extract unique countries from the renamed genome columns only
countries = set()
for original_col in genome_cols:
    new_col = rename_dict[original_col]
    country = new_col.split('_')[-1]
    countries.add(country)

dict_dfs = {'reference': {'df': df_reference, 'base_columns': base_cols_ref}, 'nonreference': {'df': df_non_reference, 'base_columns': base_cols_nonref}}

for key, df_info in dict_dfs.items():
    df = df_info['df']
    base_cols = df_info['base_columns']
    # add function to plot the PCA based on the VCF file
    plot_PCA_data(df, base_cols, key)
    # Plot proportion of the genome taken up by TEs
    df = proportion_of_TEs(df, key)

    dfs = []

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

        country_df_new = pl.DataFrame(country_df)
        # get only rows from VCF where the TE is present in that country at least once
        country_df_new = remove_TEs_not_present(country_df_new)
        # function to calculate the count of each TE type for each country
        counts_TEs_values = count_TE(country_df_new)
        
        # add country name to the counts_TEs df
        counts_TEs = counts_TEs_values.with_columns([pl.lit(f'{country}').alias("country")])
        dfs.append(counts_TEs)
        count_te_per_chr(country, country_df_new, key)

    # Concatenate all DataFrames into one big DataFrame
    df = pl.concat(dfs, how="vertical")
    plot_TE_counts(df, key)