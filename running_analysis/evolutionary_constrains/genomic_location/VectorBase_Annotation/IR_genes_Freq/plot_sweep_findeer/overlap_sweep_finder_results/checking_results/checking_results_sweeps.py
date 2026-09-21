import os
import pandas as pd

# list the files
files = os.listdir('../results/')

reference_data = {}
nonreference_data = {}

for file in files:
    if file.endswith('.csv') and '_ref' in file:
        # read the csv file
        df = pd.read_csv(f'../results/{file}')
        # store in reference_data dictionary
        reference_data[file] = df
    elif file.endswith('.csv'):
        # read the csv file
        df = pd.read_csv(f'../results/{file}')
        # store in nonreference_data dictionary
        nonreference_data[file] = df

# concatenate all dataframes in reference_data
reference_concat = pd.concat(reference_data.values(), ignore_index=True)

# concatenate all dataframes in nonreference_data
nonreference_concat = pd.concat(nonreference_data.values(), ignore_index=True)

# groupby 'chromosome, position_start, position_end, 'TE_name', 'strand' -> for the other columns merge the values with ';' separator
reference_grouped = reference_concat.groupby(['chromosome', 'position_start', 'position_end', 'TE_name', 'strand'], as_index=False).agg(lambda x: ';'.join(x.astype(str)))
nonreference_grouped = nonreference_concat.groupby(['chromosome', 'position_start', 'position_end', 'TE_name', 'strand'], as_index=False).agg(lambda x: ';'.join(x.astype(str)))

# save the grouped dataframes to csv
reference_grouped.to_csv('../results/reference_grouped.csv', index=False)
nonreference_grouped.to_csv('../results/nonreference_grouped.csv', index=False)