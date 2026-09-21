import pandas as pd
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

# make df with chromosomes - AaegL5_1, AaegL5_2, AaegL5_3 as columns and countries as rows
df = pd.DataFrame(
    columns=['AaegL5_1', 'AaegL5_2', 'AaegL5_3'],
    index = ['Brazil', 'Colombia', 'USA', 'Senegal', 'Gabon', 'Kenya']
)

# Set tolerance for LR comparison
tolerance = 1

# add information to the dataframe
for key, values in information_sweeps.items():
    chromosome, country, file_type = key
    if file_type == 'trimmed':
        trimmed = pd.read_csv(values, sep='\t')
        not_trimmed = pd.read_csv(information_sweeps[(chromosome, country, 'not_trimmed')], sep='\t')

        # Find where LR value changes from previous value for the first time (with tolerance)
        first_change_idx = None
        for i in range(1, len(not_trimmed)):
            if abs(not_trimmed.iloc[i]['LR'] - not_trimmed.iloc[i-1]['LR']) > tolerance:
                first_change_idx = i
                break
        
        # check the cutoff edges
        ## beginning chromosome
        if first_change_idx is not None:
            first_position = not_trimmed.iloc[first_change_idx]['location']
        else:
            # If no significant change found, use the first position from trimmed
            first_position = trimmed.iloc[0]['location']
        
        # ending chromosome
        last_position = trimmed.iloc[-1]['location']

        # the cutoff edges
        #print(f'For {chromosome} and {country} ---> The min position trimmed {first_position} is {first_position} and the max position trimmed is {last_position}')
        
        # add information to the dataframe
        df.loc[country, chromosome] = f'{first_position}-{last_position}'

print('cutoff for each chromosome and country:')
print(df)
# print the min value for each chromosome
for chromosome in df.columns:
    min_cutoff = float('inf')  # Start with infinity
    for country in df.index:
        if pd.notna(df.loc[country, chromosome]):
            positions = df.loc[country, chromosome].split('-')
            min_position = int(positions[0])
            if min_position < min_cutoff:
                min_cutoff = min_position
    print(f'Min cutoff for {chromosome} is {min_cutoff}')