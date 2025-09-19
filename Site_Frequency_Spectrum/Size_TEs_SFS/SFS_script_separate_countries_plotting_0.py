import pandas as pd
import matplotlib.pyplot as plt
import re
import polars as pl
import os
import numpy as np

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

def split_by_te_type_nonref(top_te, country_df):
    """
    Filters a DataFrame to retain only rows containing the specified top TE across all genomes.

    Args:
        top_te (str): The TE name to filter for. The function will look for this string in the 'TE_name' column.
        country_df (pl.DataFrame): Polars DataFrame containing TE data, with columns including 'TE_name'.

    Returns:
        pl.DataFrame: A DataFrame containing only the rows where the 'TE_name' matches the specified top TE.
    """
    # Drop the 'Ind_code' column if it exists
    country_df = country_df.drop('Ind_code')
    
    # Convert the 'TE_name' column to uppercase and split by '|' to standardize the TE names
    country_df = country_df.with_columns(pl.col('TE_name').str.split('|').list.get(0).str.to_uppercase())
    
    # Create a boolean mask to filter rows where 'TE_name' contains the specified top TE
    mask = country_df[:, 3].str.contains(top_te)
    
    # Use the boolean mask to create a new DataFrame containing only the matching rows
    selected_df = country_df.filter(mask).clone()
    
    return selected_df

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
    base_cols = ['chromosome', 'position_start', 'position_end', 'TE_name', 'strand', 'TE_length', 'TE_length_category']
    genotype_cols = [col for col in df.columns if col not in base_cols]
    
    # Count '1/1' occurrences in each row
    df['count_1/1'] = df[genotype_cols].apply(lambda row: (row == '1/1').sum(), axis=1)
    
    # Calculate frequency (proportion of samples with 1/1)
    total_samples = len(genotype_cols)
    df['frequency'] = df['count_1/1'] / total_samples
    return df

def sfs_data(frequency_of_tes, top_te):
    """
    Generate a DataFrame summarizing the occurrences of TEs based on their frequencies.

    This function calculates the frequency of each Transposable Element (TE) occurrence, 
    counts how many TEs have each frequency, and creates a DataFrame summarizing these counts.

    Parameters:
    frequency_of_tes (pd.DataFrame): DataFrame containing a 'Frequency_of_TE' column 
                                     with the frequencies of TEs.
    top_te (str): The name of the column to use as the TE name in the result DataFrame.

    Returns:
    pd.DataFrame: A DataFrame with the following columns:
                  - 'TE_name': The name of the TE (as provided by `top_te`).
                  - 'occurence': The unique frequencies of TEs.
                  - `top_te`: The count of TEs for each frequency.
    """
    frequency = frequency_of_tes['count_1/1']
    grouped = frequency.value_counts().sort_values(ascending=False)
    # Create a DataFrame with renamed columns
    result_df = pd.DataFrame({'TE_name': top_te, 'occurence': grouped.index, top_te: grouped.values})
    return result_df
    
def sfs_general_data(dataframes):
    """
    Create a DataFrame containing the frequency counts and proportions of all Transposable Elements (TEs) 
    analyzed across multiple input DataFrames.

    Parameters:
    - dataframes (list of pd.DataFrame): A list of DataFrames where each DataFrame contains frequency counts 
      for different TEs. Each DataFrame should have columns 'TE_name' and 'occurence', with the TE names 
      in the first row and occurrence counts in the second column.

    Returns:
    - pd.DataFrame: A DataFrame where each row represents a unique occurrence count, and each column represents 
      the proportion of that occurrence count for each TE.
    """
    # Find the maximum occurrence across all dataframes
    max_occurrence = max(df['occurence'].max() for df in dataframes)
    
    # Create a master DataFrame with all possible occurrences
    #master_df = pd.DataFrame({'occurence': range(1, max_occurrence + 1)})
    master_df = pd.DataFrame({'occurence': range(0, max_occurrence + 1)})
    
    # Process each DataFrame
    for df in dataframes:
        name = df.iloc[0]['TE_name']
        
        # Check for missing occurrences in the DataFrame
        #missing_occurrences = set(range(1, max_occurrence + 1)) - set(df['occurence'])
        missing_occurrences = set(range(0, max_occurrence + 1)) - set(df['occurence'])
        
        # Add missing occurrences with a count of 0
        for missing_occurrence in missing_occurrences:
            new_row = pd.DataFrame({'TE_name': [name], 'occurence': [missing_occurrence], name: [0]})
            df = pd.concat([df, new_row], ignore_index=True)

        # Sort by occurrence to ensure correct order
        df = df.sort_values(by=['occurence'])
        
        # Calculate total count for this TE
        sel_df = df[['occurence', name]]
        total_count = sel_df[name].sum()
        print(f"Total count for {name}: {total_count}")
        
        # Calculate proportions
        sel_df[name] = sel_df[name] / total_count
        
        # Now merge with the master DataFrame
        master_df = pd.merge(master_df, sel_df, on='occurence', how='left')

    print(master_df)
    return master_df

def add_other_information(master_df, neutral_data_path, country):
    """
    Adds neutral information to the master DataFrame based on the specified country.

    Parameters:
    master_df (pd.DataFrame): The DataFrame containing the master SFS data.
    neutral_data_path (str): The path to the directory containing neutral data files.
    country (str): The country to match against the neutral data files.

    Returns:
    pd.DataFrame: The master DataFrame with added neutral information as proportions.
    """
    # Dictionary mapping country names to their corresponding names in the neutral data files
    country_map = {
        "USA": "USA",
        "Colombia": "Rio_Claro",
        "Kenya": "Kenya",
        "Senegal": "Senegal",
        "Brazil": "Brazil",
        "Gabon": "Gabon"
    }

    # Get the corresponding country_name from the dictionary
    country_name_to_match = country_map.get(country, country)
    
    list_of_neutral = os.listdir(neutral_data_path)
    print("Country to match:", country)
    print("Files in neutral data path:", list_of_neutral)
    
    for sfs in list_of_neutral:
        # Use a raw string for the regex pattern
        match = re.search(r'homozygous_counts_(\w+)_rep.+mask_folded.sfs', sfs)  
        if match:
            country_name = match.group(1)
            
            if country_name_to_match == country_name:
                print(f'Matching file: {sfs} -- to this country: {country}')
                neutral = os.path.join(neutral_data_path, sfs)
                # Reading neutral info from sfs file
                neutral_df = pd.read_csv(neutral, sep=r'\s+', header=None, engine='python')
                # removes first column
                #neutral_df = neutral_df.iloc[:, 1:]
                neutral_new_df = pd.melt(neutral_df)
                print(neutral_new_df)
                neutral_new_df = neutral_new_df.rename(columns={"variable": "occurence", "value": "neutral"})
                total_neutral = neutral_new_df['neutral'].sum()
                proportions_neutral = neutral_new_df['neutral'] / total_neutral
                merged_df = pd.concat([master_df, proportions_neutral], axis=1)
                merged_df = merged_df.drop(columns=['occurence'])
                print(merged_df)
                return merged_df
            
def fold_my_sfs(master_df):
    """
    Fold the Site Frequency Spectrum (SFS) data for each column in the input DataFrame 
    except for the 'occurence' column. The folding process involves flipping the SFS 
    and adding it to itself, then retaining only the first half of the resulting array.
    NaN values are removed before the folding process, and all resulting folded SFS arrays 
    are padded with NaNs to ensure they have the same length.

    Parameters:
    - master_df (pd.DataFrame): A DataFrame where each column represents a different 
      TE's frequency counts, and the 'occurence' column represents the occurrence counts.

    Returns:
    - pd.DataFrame: A DataFrame where each column corresponds to the folded SFS of 
      the original columns from the input DataFrame. If the resulting DataFrame would be empty, 
      returns the original DataFrame with the 'occurence' column removed.
    """
    # Create a dictionary to store the folded SFS for each column
    folded_sfs_dict = {}
    max_length = 0
    
    # Iterate through columns, excluding 'occurence'
    for col in master_df.columns:
        if col != 'occurence':
            SFS = master_df[col].dropna().tolist()
            n = len(SFS)
            
            # Folded SFS array
            folded_sfs = np.zeros((n + 1) // 2)
            
            # Sum the extremes
            for i in range((n + 1) // 2):
                if i != n - i - 1:
                    folded_sfs[i] = SFS[i] + SFS[n - i - 1]
                else:
                    folded_sfs[i] = SFS[i]  # Middle element in case of odd length
            
            folded_sfs_dict[col] = folded_sfs
            max_length = max(max_length, len(folded_sfs))

    # Pad arrays with NaN to ensure all arrays have the same size
    for col in folded_sfs_dict:
        if len(folded_sfs_dict[col]) < max_length:
            folded_sfs_dict[col] = np.pad(folded_sfs_dict[col], (0, max_length - len(folded_sfs_dict[col])), constant_values=np.nan)

    # Create a DataFrame from the folded SFS dictionary
    folded_sfs_df = pd.DataFrame(folded_sfs_dict)
    
    # Check if folded_sfs_df is empty, and return the original DataFrame if true
    if folded_sfs_df.empty:
        return master_df.drop(columns=['occurence'])

    return folded_sfs_df

def plot_sfs_master_df(ax, df, suf, country):
    """
    Plots the SFS master DataFrame as a grouped bar plot.

    Parameters:
    ax (matplotlib.axes.Axes): The matplotlib Axes object to plot on.
    df (pd.DataFrame): The DataFrame containing the SFS data.
    suf (str): The suffix for the plot title.
    country (str): The country for the plot title.

    Returns:
    None
    """
    import numpy as np
    import matplotlib.pyplot as plt
    
    # Identify the variable column names dynamically
    variable_columns = [col for col in df.columns]

    # Set the width of the bars
    bar_width = 0.1

    # Create an array of indices for the x-axis positions
    x = np.arange(len(df.index))  # Adding 1 to index to start from 1

    # Create the grouped bar plot dynamically for variable columns
    colormap = plt.get_cmap('Blues', len(variable_columns))  # Use 'Blues' colormap
    
    for i, col in enumerate(variable_columns):
        ax.bar(x + (i - len(variable_columns) / 2) * bar_width, df[col], bar_width, label=col, color=colormap(i), edgecolor='black')

    # Set x-axis labels and title
    ax.set_xlabel('Number of genomes with TE insertion', fontsize=16)  # Adjust the font size
    ax.set_ylabel('Proportion', fontsize=16)  # Adjust the font size
    ax.set_xticks(x)  # Set x-axis labels to the modified 'x' array
    ax.set_title(f"{country} - {suf}", fontsize=18)

    # Adjust the legend font size and location - make legend outside the plot
    legend = ax.legend(loc='upper left')
    for item in legend.get_texts():
        item.set_fontsize(12)  # You can adjust the legend font size here

reference_path = '/work/users/g/a/gabivla/lab/SV_mosquitos/VCF_genotypes_TEs/separate_chromosomes_windows_reference/merged_reference/VCF_reference_TEs_post_processing_not600bp.csv'
neutral_data_path = '/nas/longleaf/home/gabivla/lab/SV_mosquitos/neutral_data/making_SFS/my_script/homozygous_counts_by_country/SFS_files'
df = pd.read_csv(reference_path)

# remove first column and last column
df = df.iloc[:, 1:-1]

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

fig, axs = plt.subplots(3, 3, figsize=(25, 18))

# First, create one figure per TE length category
# Get unique categories (assuming 'short', 'medium', 'long')
te_length_categories = ['small', 'medium', 'large']
outfile = open('600_df_sizes.txt', 'w')

# Create a figure for each category
for te_length_category in te_length_categories:
    print(f"\n--- Processing category: {te_length_category} ---")
    
    # Create one figure per category with 3 rows and 2 columns
    fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(12, 15))
    
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
        
        print(f"  Processing {country}")

        # calculate length of TE
        country_df['TE_length'] = country_df['position_end'] - country_df['position_start']
        country_df['TE_length'] = country_df['TE_length'].astype(int)
            
        # define 3 intervals of TE length based on the quantiles - 3
        #quantiles = country_df['TE_length'].quantile([0, 0.33, 0.66, 1]).tolist()
        # define 4 intervals of TE length based on the quantiles - 4
        #quantiles = country_df['TE_length'].quantile([0, 0.25, 0.5, 0.95, 1]).tolist()
        # Create a new column 'TE_length_category' based on the quantiles
        #country_df['TE_length_category'] = pd.cut(country_df['TE_length'], bins=quantiles, labels=['short', 'medium1', 'medium2','long'], include_lowest=True)

        # Define the bins
        bins = [float('-inf'), 200, 600, float('inf')]
        labels = ['small', 'medium', 'large']

        # Create a new column with the size category
        country_df['TE_length_category'] = pd.cut(country_df['TE_length'], bins=bins, labels=labels, right=False)
       
        # Calculate position for 2 columns, 3 rows layout
        row = i // 2
        col = i % 2

        # keep only rows with the current TE length category
        te_length_df = country_df[country_df['TE_length_category'] == te_length_category]
        # print size of the dataframe for each category - write in outfile
        outfile.write(f"Size of the dataframe for {country} {te_length_category}: {len(te_length_df)}\n")
        
        # Skip if no data for this category in this country
        if len(te_length_df) == 0:
            print(f"    No {te_length_category} data for {country}")
            continue

        sfs_data_ref = {}
        list_of_dataframes_ref = []
        
        # count TE occurrences
        for j, top_te in enumerate(TE_types):
            # Convert to Polars DataFrame here
            pl_df = pl.DataFrame(te_length_df)
            te_content_all_genomes_ref = split_by_te_type_nonref(top_te, pl_df)
            te_content_all_genomes_ref = te_content_all_genomes_ref.to_pandas()
            # calculate frequency for each TE
            frequency_TEs = calculate_frequency_TE(te_content_all_genomes_ref, country)
            # remove TE_length column
            if 'TE_length' in frequency_TEs.columns:
                frequency_TEs = frequency_TEs.drop(columns=['TE_length'])
            if 'TE_length_category' in frequency_TEs.columns:
                frequency_TEs = frequency_TEs.drop(columns=['TE_length_category'])
            sfs_data_counts_ref = sfs_data(frequency_TEs, top_te)
            sfs_data_ref[top_te] = sfs_data_counts_ref
            list_of_dataframes_ref.append(sfs_data_counts_ref)
        
        master_df_ref = sfs_general_data(list_of_dataframes_ref)
        master_df_ref = master_df_ref.drop(columns=['occurence'])
        '''all_info_df_ref = add_other_information(master_df_ref, neutral_data_path, country)
        print(all_info_df_ref)
        folded_ref_dataframe = fold_my_sfs(all_info_df_ref)'''
        
        # Plot on the current axis for this country
        plot_sfs_master_df(axs[row][col], master_df_ref, 'ref', country)

    # Add a title for the whole figure
    fig.suptitle(f'Site Frequency Spectrum - {te_length_category} TEs', fontsize=16)
    
    # Adjust spacing between subplots for better visibility
    plt.subplots_adjust(hspace=0.4, wspace=0.3)  # Increase horizontal and vertical spacing
    plt.tight_layout(rect=[0, 0, 1, 0.97])  # Make room for the suptitle
    
    # Save the plot for this specific category
    plt.savefig(f'600_SiteFrequencySpectrum_reference_{te_length_category}.svg')
    plt.close(fig)  # Close the figure to free memory