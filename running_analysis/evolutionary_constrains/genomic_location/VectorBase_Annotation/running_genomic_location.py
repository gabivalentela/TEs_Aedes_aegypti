import os
import pandas as pd

script = './genomic_location_TEs.py'
list_files = os.listdir('../../nonreference/no_duplicates/')
path_ref = '../../reference/no_duplicates/final_df_reference.csv'
for file in list_files:
    if file.endswith('1000.csv'):
        variation = file.split('_')[3]
        variation = variation.split('.')[0]
        path_file = os.path.join('../../nonreference/no_duplicates/', file)
        command = f'python {script} {path_file} {path_ref} {variation}'
        print(f'sbatch -p general --mem=20G -t 9-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="{command}"')
        os.system(f'sbatch -p general --mem=20G -t 9-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="{command}"')
