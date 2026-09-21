import pandas as pd
import matplotlib.pyplot as plt
import re
import os

general_path = '/overflow/dschridelab/users/rrlove/aedes/results/sweepfinder/outputs/merged/'
list_files = os.listdir(general_path)

information_sweeps = {}
for file in list_files:
    if 'no_relatives' in file and file.startswith('AaegL5'):
        chromosome = re.search(r'(AaegL5_\d)_no_relatives_(.+)_raw', file).group(1)
        country = re.search(r'(AaegL5_\d)_no_relatives_(.+)_raw', file).group(2)
        if 'trimmed' in file:
            full_path = os.path.join(general_path, file)
            information_sweeps[(chromosome, country, 'trimmed')] = full_path
        else:
            full_path = os.path.join(general_path, file)
            information_sweeps[(chromosome, country, 'not_trimmed')] = full_path

# Organize data by chromosome
chromosomes = ['AaegL5_1', 'AaegL5_2', 'AaegL5_3']
data_by_chromosome = {chrom: {} for chrom in chromosomes}
position_limits = {'AaegL5_1': (13191810, 310502858), 'AaegL5_2': (536150, 469057317), 'AaegL5_3': (1789994, 397127140)}

# Collect all data first
for key, values in information_sweeps.items():
    chromosome, country, file_type = key
    if file_type == 'not_trimmed':
        not_trimmed = pd.read_csv(values, sep='\t')
        #not_trimmed = pd.read_csv(information_sweeps[(chromosome, country, 'not_trimmed')], sep='\t')
        # make trimmed new version
        trimmed_new = not_trimmed[(not_trimmed['location'] >= position_limits[chromosome][0]) & (not_trimmed['location'] <= position_limits[chromosome][1])]
        trimmed_new = trimmed_new.reset_index(drop=True)
        # make tsv file with trimmed - path ./sweep_files/
        output_path = f'./new_sweeps_files//{chromosome}_no_relatives_{country}_raw_new_trimmed.tsv'
        trimmed_new.to_csv(output_path, sep='\t', index=False)
        
       