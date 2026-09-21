import os

list_files = os.listdir('../merging_files/')
for file in list_files:
    if file.endswith('.csv'):
        file_path = os.path.join('../merging_files/', file)
        command = f'python post_processing_script.py {file_path}'
        print(f'sbatch -p general --mem=20G -t 9-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="{command}"')
        os.system(f'sbatch -p general --mem=20G -t 9-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="{command}"')