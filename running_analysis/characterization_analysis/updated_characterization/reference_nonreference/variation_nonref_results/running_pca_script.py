import os
import re

list_nonref_files = '../../../../../nonreference_vcf/post_processing_nonreference_TEs/post_processing_overlap_df/max_min_updated_strategy'
list_nonref_files = [f for f in os.listdir(list_nonref_files) if f.endswith('.csv')]
for file in list_nonref_files:
    variation =  file.split('_')[-1].replace('.csv', '')
    full_path = os.path.join('../../../../../nonreference_vcf/post_processing_nonreference_TEs/post_processing_overlap_df/max_min_updated_strategy', file)
    command = 'python characterization_script_nonreference_reference.py ' + full_path + ' ' + variation
    print(f'sbatch -p general --mem=100G -t 9-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="{command}"')
    os.system(f'sbatch -p general --mem=100G -t 9-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="{command}"')