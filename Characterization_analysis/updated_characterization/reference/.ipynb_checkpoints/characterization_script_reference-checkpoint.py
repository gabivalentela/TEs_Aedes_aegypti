import os
import pandas as pd
import polars as pl
import matplotlib.pyplot as plt
import numpy as np

def count_TE(country, general_df, type_data):
    # Ensure TE names are consistently formatted
    general_df = general_df.with_columns(pl.col("TE_name").str.to_uppercase())

    # Extract clean TE names (if using a specific pattern)
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
    
    # Create a Boolean mask: True if any genome column contains '1/1' in a row
    mask = country_df.select(genome_cols).map_rows(lambda row: '1/1' in row).to_series()

    # Apply the filter to keep only rows where the TE is present at least once
    filtered_df = country_df.filter(mask)

    return filtered_df


def count_te_per_chr(country, df, type_data):
    # Ensure TE names are consistently formatted
    df = df.with_columns(pl.col("TE_name").str.to_uppercase())

    # Extract clean TE names (if using a specific pattern)
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
    color_map = cm.get_cmap('turbo', len(stacked_data))  # 'Blues' is the chosen palette

    for i, (te, counts) in enumerate(stacked_data.items()):
        ax.bar(
            unique_chromosomes, 
            counts, 
            bottom=bottom_values, 
            label=te, 
        color=color_map(i / len(stacked_data))  # Scale the colors based on the index
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
    # Calculate total percentage for each TE across countries for sorting
    te_totals = (df
        .groupby('TE_name')
        .agg(pl.col('percentage').sum())
        .sort('percentage', descending=True))

    # Get ordered unique TEs
    unique_tes = te_totals['TE_name'].to_list()
    unique_countries = df['country'].unique().to_list()
    num_countries = len(unique_countries)

    # Define a color mapping for each country (use predefined colors)
    color_mapping = {
        "Colombia": "#97F9F9",  # Blue
        "USA": "#A4DEF9",  # Orange
        "Brazil": "#7CC7F3",  # Green
        "Kenya": "#0C385A",  # Red
        "Senegal": "#4959C1",  # Red
        "Gabon": "#C59FC9"  # Red
        # Add other countries and their corresponding colors as needed
    }

    # Plotting logic
    plt.figure(figsize=(12, 8))

    # Calculate bar positions
    bar_height = 0.8
    positions = np.arange(len(unique_tes))
    width = bar_height / num_countries

    # Plotting for each country with offset positions
    for idx, country in enumerate(unique_countries):
        country_data = df.filter(pl.col('country') == country)
            
        # Create a mapping of TE names to their percentages
        te_percentages = dict(zip(country_data['TE_name'].to_list(), country_data['percentage'].to_list()))
            
        # Get percentages for all TEs (use 0 if TE not present in country)
        percentages = [te_percentages.get(te, 0) for te in unique_tes]
            
        # Calculate offset position for this country's bars
        offset = positions + (idx - num_countries / 2 + 0.5) * width
            
        # Use the predefined color for the country
        color = color_mapping.get(country, "#7f7f7f")  # Default to gray if country is not in the mapping
            
        plt.barh(offset, percentages, height=width, label=country, color=color)

    # Adjust labels and title
    plt.yticks(positions, unique_tes)
    plt.xlabel('Percentage')
    plt.title('Top TEs by Percentage per Country')

    # Add legend to differentiate countries
    plt.legend(title='Country', bbox_to_anchor=(1.05, 1), loc='upper left')

    # Adjust layout to prevent label cutoff
    plt.tight_layout()

    # Save the plot to a file
    plot_filename = f'top_TEs_by_percentage_reference.png'
    plt.savefig(plot_filename, bbox_inches='tight')
    plt.close()

import matplotlib.pyplot as plt
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import pandas as pd
from sklearn.impute import SimpleImputer

def plot_PCA_data(df, cutoff):
    # Extract only genotype columns (exclude metadata columns)
    genotype_df = df.iloc[:, 5:]  # Assuming first 6 columns are metadata
    
    # Print diagnostic information
    print("Original genotype DataFrame shape:", genotype_df.shape)
    
    # Transpose the DataFrame to have one row per TE and one column per sample
    genotype_df = genotype_df.T
    
    # Print diagnostic information after transposing
    print("Transposed genotype DataFrame shape:", genotype_df.shape)
    
    # Convert genotype strings to numerical values
    genotype_df = genotype_df.replace({'0/0': 0, '1/1': 1, './.': np.nan})
    
    # Diagnostic: Check the number of NaN values
    print("Number of NaN values:", genotype_df.isna().sum().sum())
    
    # Diagnostic: Columns with all NaN
    all_nan_columns = genotype_df.columns[genotype_df.isna().all()]
    print("Columns with all NaN values:", len(all_nan_columns))
    print(all_nan_columns)
    
    # Remove columns with all NaN values
    genotype_df = genotype_df.dropna(axis=1, how='all')
    
    print("Shape after removing all-NaN columns:", genotype_df.shape)
    
    # Impute missing values
    # Use median imputation, but first check if we have any columns left
    if genotype_df.shape[1] > 0:
        imputer = SimpleImputer(strategy='median')
        imputed_data = imputer.fit_transform(genotype_df)
        
        # Convert to DataFrame to preserve index/column information
        imputed_df = pd.DataFrame(
            imputed_data, 
            index=genotype_df.index, 
            columns=genotype_df.columns
        )
    
        # Extract country names from the row index
        countries = [row.split('_')[-1] for row in imputed_df.index]
        
        # Extract unique countries
        unique_countries = list(set(countries))
        
        # Create a color map
        country_colors = plt.cm.get_cmap('tab10', len(unique_countries))
        
        # Standardize the data
        scaler = StandardScaler()
        standardized_data = scaler.fit_transform(imputed_df)
        
        # Apply PCA
        pca = PCA(n_components=2)
        pca_result = pca.fit_transform(standardized_data)
        
        # Create a scatter plot
        plt.figure(figsize=(10, 8))
        
        # Plot each sample, color-coded by country
        for i, country in enumerate(unique_countries):
            # Create a mask for the current country
            country_mask = np.array(countries) == country
            # Select data points for the current country
            country_data = pca_result[country_mask]
            # Plot the country's data points
            plt.scatter(country_data[:, 0], country_data[:, 1], 
                        alpha=0.7, 
                        label=country, 
                        color=country_colors(i))
        
        # Labels and title
        plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.2f}% Variance)")
        plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.2f}% Variance)")
        plt.title("PCA of TE Insertions")
        plt.legend(title="Country")
        
        # Save the figure
        plt.savefig(f'PCA_reference_{cutoff}.png')
        plt.close()  # Close the plot to free up memory

def select_missing_genomes(df, cutoff):
    # Ignore the first 5 columns and count occurrences of '0/0'
    mask = (df.iloc[:, 5:] == '0/0').sum(axis=1) >= cutoff
    filtered_df = df[mask]
    return filtered_df

# Running functions
reference_path = '/work/users/g/a/gabivla/lab/SV_mosquitos/VCF_genotypes_TEs/reference/subset/VCF_TEs_reference_subset_100000_strand.csv'

# open file as df
df = pd.read_csv(reference_path)
# fixing df
df = df.drop(columns=['Unnamed: 0'])
df3 = select_missing_genomes(df, 3)
df5 = select_missing_genomes(df, 5)
df7 = select_missing_genomes(df, 7)
df10 = select_missing_genomes(df, 10)

print(f'''
      Number of TEs with at least 1 missing genomes: {df.shape[0]}
      Number of TEs with at least 3 missing genomes: {df3.shape[0]}
      Number of TEs with at least 5 missing genomes: {df5.shape[0]}
      Number of TEs with at least 7 missing genomes: {df7.shape[0]}
      Number of TEs with at least 10 missing genomes: {df10.shape[0]}''')

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
df3 = df3.rename(columns=rename_dict)
df5 = df5.rename(columns=rename_dict)
df7 = df7.rename(columns=rename_dict)
df10 = df10.rename(columns=rename_dict)

# Extract unique countries from the renamed genome columns only
countries = set()
for original_col in genome_cols:
    new_col = rename_dict[original_col]
    country = new_col.split('_')[-1]
    countries.add(country)

# add function to plot the PCA based on the VCF file
plot_PCA_data(df, 1)
plot_PCA_data(df3, 3)
plot_PCA_data(df5, 5)
plot_PCA_data(df7, 7)
plot_PCA_data(df10, 10)

'''dfs = []

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
    counts_TEs_values = count_TE(country, country_df_new, 'reference')
    
    # add country name to the counts_TEs df
    counts_TEs = counts_TEs_values.with_columns([pl.lit(f'{country}').alias("country")])
    dfs.append(counts_TEs)
    count_te_per_chr(country, country_df_new, 'ref')

# Concatenate all DataFrames into one big DataFrame
df = pl.concat(dfs, how="vertical")
plot_TE_counts(df)'''