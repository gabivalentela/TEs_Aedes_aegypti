import re
import numpy as np
import pandas as pd

def parse_slurm_fst(slurm_file):

    fst_data = {}

    pattern = (
        r'^Genome-wide weighted FST\s*'
        r'\((.+?)\s+vs\s+(.+?)\):\s*([0-9.eE+-]+)'
    )

    with open(slurm_file, 'r') as file:

        for line in file:
            if 'Unknown' in line:
                continue

            match = re.search(pattern, line.strip())

            if match:

                country_1, country_2, fst_value = match.groups()

                country_1 = country_1.strip()
                country_2 = country_2.strip()

                fst_value = float(fst_value)

                # Make pair independent of order
                pair = tuple(sorted([country_1, country_2]))

                fst_data[pair] = fst_value
    return fst_data

reference = '../reference/slurm-458158.out'
non_reference = '../nonreference/slurm-458321.out'
SNPs_data = '../SNPs/slurm-477482.out'

snp_fst = parse_slurm_fst(SNPs_data)
reference_fst = parse_slurm_fst(reference)
non_reference_fst = parse_slurm_fst(non_reference)

African_countries = ['Gabon', 'Kenya', 'Senegal']   
American_countries = ['Brazil', 'Colombia', 'USA']

# get african combinations - where both countries are in African_countries
african_combinations = {pair: fst for pair, fst in reference_fst.items() if pair[0] in African_countries and pair[1] in African_countries}
# get american combinations - where both countries are in American_countries
american_combinations = {pair: fst for pair, fst in reference_fst.items() if pair[0] in American_countries and pair[1] in American_countries}

# get the mean FST values for each continent
mean_fst_africa = np.mean(list(african_combinations.values()))
mean_fst_america = np.mean(list(american_combinations.values()))

# print the results
print('reference')
print(f"Mean FST for African countries: {mean_fst_africa}")
print(f"Mean FST for American countries: {mean_fst_america}")

# get african combinations - where both countries are in African_countries
african_combinations = {pair: fst for pair, fst in non_reference_fst.items() if pair[0] in African_countries and pair[1] in African_countries}
# get american combinations - where both countries are in American_countries
american_combinations = {pair: fst for pair, fst in non_reference_fst.items() if pair[0] in American_countries and pair[1] in American_countries}

# get the mean FST values for each continent
mean_fst_africa = np.mean(list(african_combinations.values()))
mean_fst_america = np.mean(list(american_combinations.values()))

print('non-reference')
print(f"Mean FST for African countries: {mean_fst_africa}")
print(f"Mean FST for American countries: {mean_fst_america}")