import os
import itertools
import re

ref_script = './running_final_functions.py'

# list files 
general_path = '/work/users/g/a/gabivla/lab/SV_mosquitos/VCF_genotypes_TEs/separate_chromosomes_windows_reference/reference_initial_files/'

files = os.listdir(general_path)
for file in files:
    if '_3' in file:
        path = f'''{general_path}{file}'''
        list_files = os.listdir(path)
        for window in list_files:
            print(window)
            window_path = f'''{path}/{window}'''
            chr_window = re.search('VCF_TEs_reference_(.+).csv', window).group(1)
            print(f'sbatch -p general -n 1 --cpus-per-task=4 --mem=100G -t 10-00:00:00 --mail-type=END,FAIL --mail-user=gabivla@email.unc.edu --wrap="python {ref_script} {window_path} {chr_window}"')
            os.system(f'sbatch -p general -n 1 --cpus-per-task=4 --mem=100G -t 10-00:00:00 --mail-type=END,FAIL --mail-user=gabivla@email.unc.edu --wrap="python {ref_script} {window_path} {chr_window}"')
            
        