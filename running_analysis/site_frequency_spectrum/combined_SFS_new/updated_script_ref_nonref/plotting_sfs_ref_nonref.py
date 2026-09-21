import pandas as pd
from utils import frequency_analysis

def adding_location(df):
    location_path = '/nas/longleaf/home/gabivla/lab/SV_mosquitos/location_information/mcclintock_output_locations.csv'
    location_df = pd.read_csv(location_path, usecols=['Genome_code', 'country'])
    location_dict = dict(zip(location_df['Genome_code'], location_df['country']))

    # Create rename dictionary and apply it once
    rename_dict = {col: f"{col}_{location_dict.get(col, '')}" for col in genome_cols}
    df = df.loc[:, ~df.columns.str.contains("Unnamed")]
    df = df.rename(columns=rename_dict)

    # Extract unique countries from the renamed genome columns only
    countries = set()
    for original_col in genome_cols:
        new_col = rename_dict[original_col]
        country = new_col.split('_')[-1]
        countries.add(country)

    return df, countries

def running_functions_sfs(df, countries, type_df):
    sfs_data_by_country = {}
    high_frequency = {}

    for country in countries:
        # Get all columns for this country from the renamed genome columns
        country_cols = [col for col in df.columns 
                        if col not in base_columns and col.endswith('_' + country)]
        number_of_genomes = len(country_cols)
        
        if not country_cols:
            continue
        
        country_df = df[base_columns + country_cols]
        country_df = country_df.dropna(subset=country_cols, how='all')
        frequency_TEs, high_freq_tes = frequency_analysis.calculate_frequency_TE(country_df, country)
        high_frequency[country] = high_freq_tes
        sfs_data = frequency_analysis.getting_sfs_data(frequency_TEs)
        sfs_data_by_country[country] = sfs_data

    frequency_analysis.plot_binned_sfs(sfs_data_by_country, type_df)
    # make one df with all the high frequency TEs for each country
    high_frequency_df = pd.concat(high_frequency.values(), ignore_index=True)
    # remove duplicates based on chromosome, position_start, position_end, strand, TE_name
    high_frequency_df = high_frequency_df.drop_duplicates(subset=['chromosome', 'position_start', 'position_end', 'strand', 'TE_name'])
    print(f'High frequency TEs for {type_df} dataset:')
    print(high_frequency_df)

def get_high_frequency(df_ref, df_nonref):
    # Get the frequency of TEs in reference and non-reference datasets
    frequency_ref = frequency_analysis.calculate_frequency_TE(df_ref, 'reference')
    frequency_nonref = frequency_analysis.calculate_frequency_TE(df_nonref, 'nonreference')

    # Filter for high-frequency TEs (e.g., frequency > 0.9)
    high_freq_ref = frequency_ref[frequency_ref['frequency'] >= 0.9]
    high_freq_nonref = frequency_nonref[frequency_nonref['frequency'] >= 0.9]

    print("High-frequency TEs in reference dataset:")
    print(high_freq_ref)

    print("\nHigh-frequency TEs in non-reference dataset:")
    print(high_freq_nonref)

reference_path = '/work/users/g/a/gabivla/lab/SV_mosquitos/TE_paper_aedesaegytpi/reference_vcf/genotyping_reference_TEs/post_processing/merged_results_5_2_either_filtered.csv'
nonreference_path = '/work/users/g/a/gabivla/lab/SV_mosquitos/TE_paper_aedesaegytpi/nonreference_vcf/post_processing_nonreference_TEs/post_processing_overlap_df/max_min_updated_strategy/final_df_overlap_retroseq_teforest_nonreference_1000.csv'

df_ref = pd.read_csv(reference_path)
df_nonref = pd.read_csv(nonreference_path)
# drop score column from df_ref
df_ref = df_ref.drop(columns=['score'])
#get_high_frequency(df_ref, df_nonref)
# base columns 
base_columns = ['chromosome', 'position_start', 'position_end', 'strand', 'TE_name']
# update genome names in df_ref to match df_nonref
genome_cols = [col for col in df_ref.columns if col not in base_columns]

# update the genome columns names
for genome in genome_cols:
    if genome.startswith('SRR'):
        genome_updated = genome.replace('_1', '')
    elif genome.startswith('FEMALE') or genome.startswith('MALE'):
        genome_updated = genome.replace('_L001_L002_R1_001', '')
    # update the column name in the dataframe
    df_ref = df_ref.rename(columns={genome: genome_updated})

genome_cols = [col for col in df_ref.columns if col not in base_columns]

updated_df_ref, countries_ref = adding_location(df_ref)
updated_df_nonref, countries_nonref = adding_location(df_nonref)

# running it for ref and nonref
print(f'running it for reference')
running_functions_sfs(updated_df_ref, countries_ref, 'reference')
print(f'running it for nonreference')
running_functions_sfs(updated_df_nonref, countries_nonref, 'nonreference')