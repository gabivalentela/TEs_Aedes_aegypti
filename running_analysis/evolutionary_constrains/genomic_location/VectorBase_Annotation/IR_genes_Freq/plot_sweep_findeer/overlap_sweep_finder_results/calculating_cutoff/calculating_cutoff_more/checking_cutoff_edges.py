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
        
        # Find where LR value changes from next value for the last time (starting from end)
        last_change_idx = None
        for i in range(len(not_trimmed) - 2, -1, -1):
            if abs(not_trimmed.iloc[i]['LR'] - not_trimmed.iloc[i+1]['LR']) > tolerance:
                last_change_idx = i
                break
        
        # Find first position where LR < 50
        first_lr_below_100_idx = None
        for i in range(len(not_trimmed)):
            if not_trimmed.iloc[i]['LR'] < 10:
                first_lr_below_100_idx = i
                break
        
        # Find last position where LR < 100 (starting from end)
        last_lr_below_100_idx = None
        for i in range(len(not_trimmed) - 1, -1, -1):
            if not_trimmed.iloc[i]['LR'] < 10:
                last_lr_below_100_idx = i
                break
        
        # check the cutoff edges
        ## beginning chromosome
        if first_change_idx is not None:
            first_position_change = not_trimmed.iloc[first_change_idx]['location']
        else:
            # If no significant change found, use the first position from trimmed
            first_position_change = trimmed.iloc[0]['location']
        
        if first_lr_below_100_idx is not None:
            first_position_lr100 = not_trimmed.iloc[first_lr_below_100_idx]['location']
        else:
            first_position_lr100 = None
        
        # Take the maximum of the two starting positions (most conservative)
        if first_position_lr100 is not None:
            first_position = max(first_position_change, first_position_lr100)
        else:
            first_position = first_position_change
        
        ## ending chromosome
        if last_change_idx is not None:
            last_position_change = not_trimmed.iloc[last_change_idx]['location']
        else:
            # If no significant change found, use the last position from trimmed
            last_position_change = trimmed.iloc[-1]['location']
        
        if last_lr_below_100_idx is not None:
            last_position_lr100 = not_trimmed.iloc[last_lr_below_100_idx]['location']
        else:
            last_position_lr100 = None
        
        # Take the minimum of the two ending positions (most conservative)
        if last_position_lr100 is not None:
            last_position = min(last_position_change, last_position_lr100)
        else:
            last_position = last_position_change
        
        # add information to the dataframe
        df.loc[country, chromosome] = f'{first_position}-{last_position}'

print('cutoff for each chromosome and country:')
print(df)

# print the max value for each chromosome (start position)
print('\nMax start cutoff for each chromosome:')
for chromosome in df.columns:
    max_cutoff = 0
    for country in df.index:
        if pd.notna(df.loc[country, chromosome]):
            positions = df.loc[country, chromosome].split('-')
            start_position = int(positions[0])
            if start_position > max_cutoff:
                max_cutoff = start_position
    print(f'Max start cutoff for {chromosome} is {max_cutoff}')

# print the min value for each chromosome (end position)
print('\nMin end cutoff for each chromosome:')
for chromosome in df.columns:
    min_cutoff = float('inf')
    for country in df.index:
        if pd.notna(df.loc[country, chromosome]):
            positions = df.loc[country, chromosome].split('-')
            end_position = int(positions[1])
            if end_position < min_cutoff:
                min_cutoff = end_position
    print(f'Min end cutoff for {chromosome} is {min_cutoff}')