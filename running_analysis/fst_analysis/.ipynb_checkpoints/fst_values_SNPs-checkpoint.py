import pandas as pd
from itertools import combinations
import numpy as np
import gzip
import matplotlib.pyplot as plt

def read_vcf(vcf_file):
    with gzip.open(vcf_file, 'rt') as f:  # 'rt' mode for reading text files
        # Initialize header as None before finding it
        header = None
        # Read through the file and get the header line
        for line in f:
            if line.startswith("#CHROM"):
                header = line.strip("#").strip().split("\t")
                break  # Exit the loop once the header is found
        
        # Read the rest of the VCF file into a DataFrame, skipping comment lines
        df = pd.read_csv(f, sep='\t', comment='#', header=None, names=header)

    # Add a new 'info' column based on 'CHROM' and 'POS'
    df['info'] = df['CHROM'] + '_' + df['POS'].astype(str)
    
    # Keep only the portion before the first ":" in all columns except 'info'
    #df = df.applymap(lambda x: x.split(':')[0] if isinstance(x, str) else x)
    df = df.map(lambda x: x.split(':')[0] if isinstance(x, str) else x)

    # Drop unwanted columns
    columns_to_remove = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'FORMAT', 'INFO']
    df_selected = df.drop(columns=columns_to_remove)
    
    return df_selected

def info_location(loc_info):
    """
    Reads the location information from a CSV file into a Polars DataFrame.

    Args:
        loc_info (str): Path to the CSV file containing location information.

    Returns:
        pl.DataFrame: Polars DataFrame with location information.
    """
    df = pd.read_csv(loc_info)
    return df

def split_by_countries(df, loc_info):
    # Create a dictionary with the genome code as the key and the country as the value
    loc_info_dict = loc_info.set_index('Genome_code')['country'].to_dict()

    # Prepare the dictionary to store results
    country_dict = {}
    
    # Iterate over each column in the dataframe (each column corresponds to a genome)
    for col in df.columns:
        # Get the genome code from the column name (assuming column names are genome codes)
        genome_code = col
        
        # Get the country for the current genome
        country = loc_info_dict.get(genome_code, 'Cali')  # Default to 'Unknown' if the genome code isn't in loc_info_dict
        
        # Add a new column to the dictionary for this country if it doesn't already exist
        if country not in country_dict:
            country_dict[country] = df.loc[:, [col]]  # Add the entire column for this country
        else:
            country_dict[country] = pd.concat([country_dict[country], df.loc[:, [col]]], axis=1)  # Add the column to the existing country data

    return country_dict

def calculate_al_freq(df):
    # get the number of samples
    number_of_samples = df.shape[1] - 1
    # calculate frequency of first allele - 1/1
    df['1/1'] = df.iloc[:, 1:].apply(lambda x: (x == '1/1').sum() / number_of_samples, axis=1)
    # calculate frequency of second allele - 0/1
    df['0/0'] = df.iloc[:, 1:].apply(lambda x: (x == '0/0').sum() / number_of_samples, axis=1)

    return df

def calculate_FST(freq_p1, freq_p2):
    p1 = freq_p1['1/1']
    p2 = freq_p2['1/1']

    # calculate q1 and q2
    q1 = 1 - p1
    q2 = 1 - p2

    # calculate total allele frequency
    p_t = (p1 + p2) / 2
    q_t = 1 - p_t

    # calculate expected heterozygosity
    # first calculate expected heterozygosity for the two populations
    # pop1
    hs_1 = 2 * p1 * q1
    # pop2
    hs_2 = 2 * p2 * q2
    # then take the mean of this
    hs = (hs_1 + hs_2) / 2

    # next calculate expected heterozygosity for the metapopulations
    ht = 2 * p_t * q_t

    # calculate fst
    fst = (ht - hs) / ht

    # return output
    return fst

vcf_SNPs = '/work/users/g/a/gabivla/lab/SV_mosquitos/neutral_data_vcf/updated_AaegL5_full_filtered_snps_110122_rep_map_100kbgenes_centromere_filtered_sites_mask.vcf.gz'
#vcf_SNPs = './subset_VCF/output_file.vcf.gz'

df = read_vcf(vcf_SNPs)

#location info
loc_info = info_location('/nas/longleaf/home/gabivla/lab/SV_mosquitos/location_information/mcclintock_output_locations.csv')

vcf_files = split_by_countries(df, loc_info)

country_vcf_dict = {}

country_pairs = list(combinations(vcf_files.keys(), 2))

outfile_nonref = open('FST_values_SNPs.txt', 'a')

for country1, country2 in country_pairs:
    if country1 == 'Cali' or country2 == 'Cali':
        pass
    else:
        country_vcf_dict[(country1, country2)] = (vcf_files[country1], vcf_files[country2])

# Create a figure with subplots
fig, axes = plt.subplots(4, 4, figsize=(40, 30))
i = 0 
for pair, vcf_tuple in country_vcf_dict.items():
    country1, country2 = pair
    vcf_country_1, vcf_country_2 = vcf_tuple
        
    #calculate frequency of alleles
    counts_1 = calculate_al_freq(vcf_country_1)
    counts_2 = calculate_al_freq(vcf_country_2)

    fst_pop1_pop2 = calculate_FST(counts_1, counts_2)
    fst_pop1_pop2.to_csv(f'/work/users/g/a/gabivla/lab/SV_mosquitos/FST_data/SNPs/SNP_FST_values_{country1}_{country2}.csv')
    # Get the FST values as a list or numpy array
    fst_values = fst_pop1_pop2.values

    # Compute the histogram counts and bin edges
    counts, bins = np.histogram(fst_values, bins=20, range=(0, 1))

    # Convert counts to percentages
    total_count = sum(counts)
    percentages = (counts / total_count) * 100

    # Plot the histogram as percentages
    axes[i // 4, i % 4].bar(bins[:-1], percentages, width=np.diff(bins), edgecolor='black', alpha=0.7)

    # Label the plot
    axes[i // 4, i % 4].set_title(f'FST Histogram: {country1} vs {country2}')
    axes[i // 4, i % 4].set_xlabel('FST value')
    axes[i // 4, i % 4].set_ylabel('Percentage')
    axes[i // 4, i % 4].set_xlim(0, 1)
    axes[i // 4, i % 4].set_ylim(0, 100)
    
    fst_value = np.mean(fst_pop1_pop2)
        
    outfile_nonref.write(f'FST value between {country1} and {country2}: {fst_value}\n')
    
    i += 1
# Adjust layout
plt.tight_layout()

# Save the figure
plt.savefig('fst_histograms_SNPs.png')