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
    country_df = country_df.drop('Ind_code')
    
    country_df = country_df.with_columns(pl.col('TE_name').str.split('|').list.get(0).str.to_uppercase())
    
    mask = country_df[:, 3].str.contains(top_te)
    
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
    base_cols = ['chromosome', 'position_start', 'position_end', 'TE_name', 'strand']
    genotype_cols = [col for col in df.columns if col not in base_cols]
    
    df['count_1/1'] = df[genotype_cols].apply(lambda row: (row == '1/1').sum(), axis=1)
    
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
    
    max_occurrence = max(df['occurence'].max() for df in dataframes)
    
    master_df = pd.DataFrame({'occurence': range(1, max_occurrence + 1)})
    
    for df in dataframes:
        name = df.iloc[0]['TE_name']
        
        missing_occurrences = set(range(1, max_occurrence + 1)) - set(df['occurence'])
        
        for missing_occurrence in missing_occurrences:
            new_row = pd.DataFrame({'TE_name': [name], 'occurence': [missing_occurrence], name: [0]})
            df = pd.concat([df, new_row], ignore_index=True)

        df = df.sort_values(by=['occurence'])
        
        sel_df = df[['occurence', name]]
        total_count = sel_df[name].sum()
        print(f"Total count for {name}: {total_count}")
        
        sel_df[name] = sel_df[name] / total_count
        
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
    country_map = {
        "USA": "USA",
        "Colombia": "Rio_Claro",
        "Kenya": "Kenya",
        "Senegal": "Senegal",
        "Brazil": "Brazil",
        "Gabon": "Gabon"
    }

    country_name_to_match = country_map.get(country, country)
    
    list_of_neutral = os.listdir(neutral_data_path)
    print("Country to match:", country)
    print("Files in neutral data path:", list_of_neutral)
    
    for sfs in list_of_neutral:
        match = re.search(r'homozygous_counts_(\w+)_rep.+mask_folded.sfs', sfs)  
        if match:
            country_name = match.group(1)
            
            if country_name_to_match == country_name:
                print(f'Matching file: {sfs} -- to this country: {country}')
                neutral = os.path.join(neutral_data_path, sfs)
                neutral_df = pd.read_csv(neutral, sep=r'\s+', header=None, engine='python')
                neutral_df = neutral_df.iloc[:, 1:]
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
    folded_sfs_dict = {}
    max_length = 0
    
    for col in master_df.columns:
        if col != 'occurence':
            SFS = master_df[col].dropna().tolist()
            n = len(SFS)
            
            folded_sfs = np.zeros((n + 1) // 2)
            
            for i in range((n + 1) // 2):
                if i != n - i - 1:
                    folded_sfs[i] = SFS[i] + SFS[n - i - 1]
                else:
                    folded_sfs[i] = SFS[i] 
            
            folded_sfs_dict[col] = folded_sfs
            max_length = max(max_length, len(folded_sfs))

    for col in folded_sfs_dict:
        if len(folded_sfs_dict[col]) < max_length:
            folded_sfs_dict[col] = np.pad(folded_sfs_dict[col], (0, max_length - len(folded_sfs_dict[col])), constant_values=np.nan)

    folded_sfs_df = pd.DataFrame(folded_sfs_dict)
    
    if folded_sfs_df.empty:
        return master_df.drop(columns=['occurence'])

    return folded_sfs_df

def plot_binned_sfs(combined_df, suf):
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    combined_df['binned_occurence'] = pd.cut(combined_df['occurence'], bins=10, labels=False)

    country_colors = {
        'Brazil': 'green',
        'USA': 'red',
        'Colombia': 'yellow',
        'Senegal': 'blue',
        'Gabon': 'orange',
        'Kenya': 'pink'
    }

    countries = [c for c in country_colors if c in combined_df['Country'].unique()]
    bar_width = 0.1
    fig, ax = plt.subplots(figsize=(14, 6))

    x_bins = np.arange(10)
    total_countries = len(countries)

    for i, country in enumerate(countries):
        country_df = combined_df[combined_df['Country'] == country]
        binned = country_df.groupby('binned_occurence')[['TEs']].mean().reindex(range(10), fill_value=0)

        offset = (i - total_countries / 2) * bar_width * 1.5
        x_te = x_bins + offset

        ax.bar(x_te, binned['TEs'], width=bar_width, color=country_colors[country], label=f"{country} TEs")

    ax.set_xticks(x_bins)
    ax.set_xticklabels([f"{i/10:.1f}-{(i+1)/10:.1f}" for i in range(10)])
    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=15)
    ax.set_xlabel('Proportion of genomes with insertion (binned)', fontsize=15)
    ax.set_ylabel('Proportion', fontsize=15)
    ax.set_title(f"Binned Site Frequency Spectrum - {suf}")
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
    plt.tight_layout()
    plt.savefig(f'binned_SiteFrequencySpectrum_onlyTEs_{suf}.svg')

nonreference_path = '/work/users/g/a/gabivla/lab/SV_mosquitos/VCF_genotypes_TEs/separate_chromosomes_windows_reference/merged_nonreference/VCF_nonreference_post_processing.csv'

neutral_data_path = '/nas/longleaf/home/gabivla/lab/SV_mosquitos/neutral_data/making_SFS/my_script/homozygous_counts_by_country/SFS_files'
df = pd.read_csv(nonreference_path)
df = df.iloc[:, 1:]

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

df = pl.DataFrame(df)
TE_types = count_TE_occurence(df)
TE_types = TE_types['TE_name'].to_list()
df = df.to_pandas()

fig, axs = plt.subplots(4, 4, figsize=(25, 18))

all_info_by_country = []

for country in countries:
    country_cols = [col for col in df.columns 
                    if col not in base_cols and col.endswith('_' + country)]
    number_of_genomes = len(country_cols)
    
    if not country_cols:
        continue
    
    country_df = df[base_cols + country_cols]
    country_df = country_df.dropna(subset=country_cols, how='all')

    frequency_TEs = calculate_frequency_TE(country_df, country)
    sfs_data_counts_ref = sfs_data(frequency_TEs, 'TEs')
    list_of_dataframes_ref = [sfs_data_counts_ref]
    master_df_ref = sfs_general_data(list_of_dataframes_ref)
    all_info_df_ref = add_other_information(master_df_ref, neutral_data_path, country)
    # make extra column with occurrence from the index
    all_info_df_ref['occurence'] = all_info_df_ref.index
    # divide occurrence by number of genomes
    all_info_df_ref['occurence'] = all_info_df_ref['occurence'] / number_of_genomes
    all_info_df_ref['Country'] = country  # Add a country column
    all_info_by_country.append(all_info_df_ref)

combined_df = pd.concat(all_info_by_country, ignore_index=True)
print(combined_df)
plot_binned_sfs(combined_df, 'nonref')