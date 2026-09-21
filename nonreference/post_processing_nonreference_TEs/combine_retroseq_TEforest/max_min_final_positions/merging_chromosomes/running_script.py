import os

script = 'merging_chromosomes.py'
variation = ['100', '500','1000']

for var in variation:
    command = f'python {script} {var}'
    print(f'sbatch -p general --mem=100G -t 9-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="{command}"')
    os.system(f'sbatch -p general --mem=100G -t 9-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="{command}"')