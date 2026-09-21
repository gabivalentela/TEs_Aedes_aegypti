import pandas as pd
import os
import re

general_path = '../make_new_sweep_files/sweep_files/'

list_files = os.listdir(general_path)
data_list = []

for file in list_files:
    complete_path = os.path.join(general_path, file)
    # country
    country = re.search(r'../make_new_sweep_files/sweep_files/(AaegL5_\d)_no_relatives_(.+)_raw_new_trimmed.tsv', complete_path).group(2)
    # chromosome
    chromosome = re.search(r'../make_new_sweep_files/sweep_files/(AaegL5_\d)_no_relatives_(.+)_raw_new_trimmed.tsv', complete_path).group(1)
    
    # open as df
    df = pd.read_csv(complete_path, sep='\t')
    # calculate threshold
    threshold_99 = df['LR'].quantile(0.95)
    
    # add to list
    data_list.append({
        'country': country,
        'chromosome': chromosome,
        'LR_99_threshold': threshold_99
    })

# Create dataframe from list
final_df = pd.DataFrame(data_list)

# Pivot to get countries as rows and chromosomes as columns
final_df_pivot = final_df.pivot(index='country', columns='chromosome', values='LR_99_threshold')

print(final_df_pivot)