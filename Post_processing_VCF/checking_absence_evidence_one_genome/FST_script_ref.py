import pandas as pd
from itertools import combinations

def calculate_frequency_TE(df, country_name):
    """
    Calculate the frequency of '1/1' genotypes for each TE insertion.
    
    Parameters:
    df - DataFrame with TE data
    country_name - Name of the country (for reference, not used in calculation)
    
    Returns:
    DataFrame with frequency column added
    """
    # Get columns that contain genotype data
    base_cols = ['chromosome', 'position_start', 'position_end', 'TE_name', 'strand']
    genotype_cols = [col for col in df.columns if col not in base_cols]
    
    # Count '1/1' occurrences in each row
    df['count_1/1'] = df[genotype_cols].apply(lambda row: (row == '1/1').sum(), axis=1)
    df['count_0/0'] = df[genotype_cols].apply(lambda row: (row == '0/0').sum(), axis=1)

    # Calculate frequency
    total_samples = df['count_1/1'] + df['count_0/0']
    df['frequency'] = df['count_1/1'] / total_samples
    
    return df

def calculate_FST(freq_p1, freq_p2):
    # Align DataFrames by CHROM_POS_END
    merged = pd.merge(freq_p1, freq_p2, on=["position_start", 'position_end', 'TE_name', 'chromosome', 'strand'], suffixes=('_p1', '_p2'))
    
    # remove rows were both country_1/1_p1 and country_1/1_p2 are 0
    merged = merged[(merged['count_1/1_p1'] > 0) | (merged['count_1/1_p2'] > 0)]

    # Keep the first 5 columns
    cols_to_keep = list(merged.columns[:5])

    # Add any columns that contain "count_1/1"
    cols_to_keep += [col for col in merged.columns if "frequency_" in col]

    # Subset the DataFrame
    merged = merged[cols_to_keep]

    # Calculate allele frequencies
    p1 = merged['frequency_p1']
    p2 = merged['frequency_p2']

    # Calculate q1 and q2
    q1 = 1 - p1
    q2 = 1 - p2

    # Calculate total allele frequency
    p_t = (p1 + p2) / 2
    q_t = 1 - p_t

    # Calculate expected heterozygosity
    hs_1 = 2 * p1 * q1
    hs_2 = 2 * p2 * q2
    hs = (hs_1 + hs_2) / 2

    ht = 2 * p_t * q_t

    # Calculate FST
    fst = (ht - hs) / ht

    # Return the five specified columns and FST as a DataFrame
    fst_df = pd.DataFrame({
        'chromosome': merged['chromosome'],
        'position_start': merged['position_start'],
        'position_end': merged['position_end'],
        'TE_name': merged['TE_name'],
        'strand': merged['strand'],
        'FST': fst
    })
    
    return fst_df

# Running functions
reference_path = '/work/users/g/a/gabivla/lab/SV_mosquitos/VCF_genotypes_TEs/separate_chromosomes_windows_reference/merged_reference/VCF_reference_TEs_post_processing_600.csv'

# open file as df 
df = pd.read_csv(reference_path)

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
    all_dfs_ref[country] = frequency_df

for country1, country2 in country_pairs_nonref:
    country_vcf_dict_ref[(country1, country2)] = (all_dfs_ref[country1], all_dfs_ref[country2])

# Iterate over country pairs and their corresponding VCF data
for i, (pair, vcf_tuple) in enumerate(country_vcf_dict_ref.items()):
    country1, country2 = pair
    vcf_country_1, vcf_country_2 = vcf_tuple

    # Calculate FST values
    fst_pop1_pop2 = calculate_FST(vcf_country_1, vcf_country_2)
    fst_pop1_pop2.to_csv(f'/work/users/g/a/gabivla/lab/SV_mosquitos/FST_analysis/TEs_updated_frequency/FST_values_{country1}_{country2}_ref.csv', index=False)