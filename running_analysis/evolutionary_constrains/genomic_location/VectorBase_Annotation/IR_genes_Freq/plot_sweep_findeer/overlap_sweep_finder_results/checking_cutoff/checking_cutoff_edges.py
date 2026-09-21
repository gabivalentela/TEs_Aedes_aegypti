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

# add information to the dataframe
for key, values in information_sweeps.items():
    chromosome, country, file_type = key
    if file_type == 'trimmed':
        trimmed = pd.read_csv(values, sep='\t')
        not_trimmed = pd.read_csv(information_sweeps[(chromosome, country, 'not_trimmed')], sep='\t')
        print(values)
        print(trimmed)

        # check the cutoff edges
        ## beginning chromosome
        first_postion = trimmed.iloc[0]
        first_postion = first_postion['location']
        # ending chromosome
        last_position = trimmed.iloc[-1]
        last_position = last_position['location']

        # the cutoff edges
        #print(f'For {chromosome} and {country} ---> The min position trimmed {first_postion} is {first_postion} and the max position trimmed is {last_position}')
        
        # add information to the dataframe
        df.loc[country, chromosome] = f'{first_postion}-{last_position}'
print('cutoff for each chromosome and country:')
print(df)