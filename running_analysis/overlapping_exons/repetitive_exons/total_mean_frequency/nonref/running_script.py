import os

variation = ['1000']
for variable in variation:
    command = f'python non_ref_overlap.py {variable}'
    print(f'sbatch -p general --mem=100G -t 9-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="{command}"')    
    os.system(f'sbatch -p general --mem=100G -t 9-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="{command}"')