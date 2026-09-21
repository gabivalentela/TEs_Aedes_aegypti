import pandas as pd
import os
import re

# Your existing code for loading data
list_sweep_path = '/overflow/dschridelab/users/rrlove/aedes/results/sweepfinder/outputs/merged/'
ref_df = '/work/users/g/a/gabivla/lab/SV_mosquitos/VCF_genotypes_TEs/separate_chromosomes_windows_reference/merged_reference/VCF_reference_TEs_post_processing_600_absence_evidence.csv'
nonref_df = '/work/users/g/a/gabivla/lab/SV_mosquitos/VCF_genotypes_TEs/separate_chromosomes_windows_reference/merged_nonreference/VCF_nonreference_post_processing.csv'
script_running = './find_overlap_sweep_results.py'

# dictionary to store country -> chromosome -> file path mappings
sweep_information = {}
# open list files
with open('./list_files_sweepfinder.txt', 'r') as f:
    file_lines = f.read().splitlines()
    for file in file_lines:
        country = re.search(r'.+no_relatives_(.+)_raw', file).group(1)
        chrom = re.search(r'(AaegL5_\d+)_no_relatives.+', file).group(1)
        file = file.replace('\t', '')
        
        # Initialize country dict if it doesn't exist
        if country not in sweep_information:
            sweep_information[country] = {}
        
        # Add chromosome to the country's dict
        sweep_information[country][chrom] = f'{list_sweep_path}{file}'

for country, files in sweep_information.items():
    for chrom, file in files.items():
        print(f'Country: {country}, Chromosome: {chrom}, File: {file}')
        # run script - ref
        command = f'python {script_running} {ref_df} {file} ref'
        print(f'{command}')
        os.system(f'sbatch -p general -n 1 --mem=20G -t 10-20:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="{command}"')
        # run script - nonref
        command = f'python {script_running} {nonref_df} {file} nonref'
        print(f'{command}')
        os.system(f'sbatch -p general -n 1 --mem=20G -t 10-20:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="{command}"')