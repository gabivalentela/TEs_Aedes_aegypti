import os
import re

general_path = '../../combine_retroseq_TEforest/max_min_final_positions/merging_chromosomes'
list_files = os.listdir(general_path)

for file in list_files:
    if file.endswith('.csv'):
        path_file = os.path.join(general_path, file)
        print(path_file)
        variation = re.search(r'.+teforest_(.+)_all_chromosomes.csv', file).group(1)
        command = f'python checking_VCF_nonref.py {path_file} {variation}'
        print(f'sbatch -p general --mem=20G -t 9-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="{command}"')
        os.system(f'sbatch -p general --mem=20G -t 9-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="{command}"')