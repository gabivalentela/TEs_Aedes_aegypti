import os
import pandas as pd
import polars as pl
import matplotlib.pyplot as plt
import numpy as np

def count_TE(country, general_df, type_data):
    # Ensure TE names are consistently formatted
    general_df = general_df.with_columns(pl.col("TE_name").str.to_uppercase())

    # Extract clean TE names
    df_te = general_df.with_columns(
        pl.col("TE_name").str.extract(r'RND_\d+_FAMILY_\d+_([A-Z]+(?:_HELITRON)?)', 1).alias("TE_name")
    )

    # Count occurrences of each TE name
    grouped = df_te.group_by('TE_name').agg(pl.len().alias("count"))
    print(f'this is the grouped df: {grouped}')
    # Add percentage column
    total_count = grouped["count"].sum()
    grouped = grouped.with_columns((pl.col("count") / total_count * 100).alias("percentage"))

    # Sort by percentage in descending order
    sorted_grouped = grouped.sort("percentage", descending=True)

    # Get top 10 most present TEs
    df_first_10 = sorted_grouped.head(10)

    print(f'In this country, these are the top 10 TEs by percentage: {country}\n{df_first_10}')

    return df_first_10

def remove_TEs_not_present(country_df):
    # Identify genome columns
    genome_cols = [col for col in country_df.columns if col not in ['chromosome', 'position_start', 'position_end', 'TE_name', 'strand']]
    
    mask = country_df.select(genome_cols).map_rows(lambda row: '1/1' in row).to_series()

    filtered_df = country_df.filter(mask)

    return filtered_df


def count_te_per_chr(country, df, type_data):
    # Ensure TE names are consistently formatted
    df = df.with_columns(pl.col("TE_name").str.to_uppercase())

    # Extract clean TE names
    df = df.with_columns(
        pl.col("TE_name")
        .str.extract(r'RND_\d+_FAMILY_\d+_([A-Z]+(?:_HELITRON)?)', 1)
        .alias("TE_name")
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

    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import numpy as np

    # Create the stacked bar chart
    fig, ax = plt.subplots(figsize=(10, 6))
    bottom_values = np.zeros(len(unique_chromosomes))

    # Generate a color map using the 'Blues' palette
    color_map = cm.get_cmap('turbo', len(stacked_data))  
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
    
def plot_TE_counts(df):
    import matplotlib.pyplot as plt
    import numpy as np

    # Fixed country order
    country_order = ["Brazil", "USA", "Colombia", "Senegal", "Gabon", "Kenya"]

    # Calculate total percentage for each TE across countries for sorting
    te_totals = (
        df.groupby('TE_name')
        .agg(pl.col('percentage').sum())
        .sort('percentage', descending=True)
    )

    # Get ordered unique TEs
    unique_tes = te_totals['TE_name'].to_list()
    num_countries = len(country_order)

    # Define a color mapping for each country
    color_mapping = {
        "Brazil":   "#009E73",  # Teal Green
        "USA":      "#D73027",  # Strong Red
        "Colombia": "#E69F00",  # Orange
        "Senegal":  "#0072B2",  # Blue
        "Gabon":    "#D55E00",  # Vermillion
        "Kenya":    "#CC79A7"   # Reddish Purple
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

    plt.yticks(positions, unique_tes, fontsize=20)
    plt.xticks(fontsize=20)
    plt.xlabel('Percentage', fontsize=20)
    plt.title('Top TEs by Percentage per Country', fontsize=22)
    plt.legend(title='Country', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=20)
    plt.tight_layout()
    # Save the plot to a file
    plot_filename = f'top_TEs_by_percentage_reference.svg'
    plt.savefig(plot_filename, bbox_inches='tight')
    plt.close()

import matplotlib.pyplot as plt
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import pandas as pd

def plot_PCA_data(df):
    # Extract only genotype columns
    genotype_df = df.iloc[:, 6:]  # Assuming first 6 columns are metadata
    
    # Determine the most frequent genotype for each TE
    def fill_missing(row):
        most_common = row[row != './.'].mode()  
        most_common = most_common[0] if not most_common.empty else '0/0' 
        return row.replace({'./.': most_common})

    # Apply the function to fill missing values
    genotype_df = genotype_df.apply(fill_missing, axis=1)
    genotype_df = genotype_df.T
    
    # Convert genotype strings to numerical values
    genotype_df = genotype_df.replace({'0/0': 0, '1/1': 1})
    genotype_df = genotype_df.infer_objects(copy=False).astype(int)
    
    # Extract country names from the row index
    countries = [row.split('_')[-1] for row in genotype_df.index]

    # Extract unique countries
    unique_countries = sorted(set(countries))  
    color_mapping = {
        "Colombia": "#E69F00",  # Orange
        "USA":      "#D73027",  # Strong Red
        "Brazil":   "#009E73",  # Teal Green
        "Kenya":    "#CC79A7",  # Reddish Purple
        "Senegal":  "#0072B2",  # Blue
        "Gabon":    "#D55E00"   # Vermillion
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
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.2f}% Variance)", fontsize=20)
    # make font bigger for the x and y ticks
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.2f}% Variance)", fontsize=20)
    ax.set_title("PCA of TE Insertions", fontsize=16)

    # Put legend outside the figure on the right
    fig.legend(title="Country", loc='center left', bbox_to_anchor=(1.01, 0.5), fontsize=16, title_fontsize=18)

    # Adjust layout to make room for the legend
    plt.tight_layout(rect=[0, 0, 0.85, 1])

    # Save the figure
    plt.savefig('PCA_reference_new.svg', bbox_inches='tight')
    plt.close()

def select_more_than_3(df, cutoff):
    # Ignore the first 5 columns and count occurrences of '1/1'
    mask = (df.iloc[:, 6:] == '1/1').sum(axis=1) >= cutoff
    filtered_df = df[mask]
    return filtered_df

# Running functions
reference_path = '/work/users/g/a/gabivla/lab/SV_mosquitos/VCF_genotypes_TEs/separate_chromosomes_windows_reference/merged_reference/VCF_reference_TEs_post_processing_600_absence_evidence.csv'

# open file as df
df = pd.read_csv(reference_path)

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

# add function to plot the PCA based on the VCF file
plot_PCA_data(df)

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
    print(country_df.columns)

    country_df_new = pl.DataFrame(country_df)
    # get only rows from VCF where the TE is present in that country at least once
    country_df_new = remove_TEs_not_present(country_df_new)
    # function to calculate the count of each TE type for each country
    counts_TEs_values = count_TE(country, country_df_new, 'reference')
    
    # add country name to the counts_TEs df
    counts_TEs = counts_TEs_values.with_columns([pl.lit(f'{country}').alias("country")])
    dfs.append(counts_TEs)
    count_te_per_chr(country, country_df_new, 'ref')

# Concatenate all DataFrames into one big DataFrame
df = pl.concat(dfs, how="vertical")
print(df)
plot_TE_counts(df)