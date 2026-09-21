import os

script_python = './combine_teforest_retroseq_vcf.py'
variation_values = ['100', '500', '1000']
chromosome = ['AaegL5_1', 'AaegL5_2', 'AaegL5_3']
#variation_values = ['100']

for variation in variation_values:
    for chrom in chromosome:
        command = f'python {script_python} {variation} {chrom}'
        #print(f'Running command: {command}')
        print(f'sbatch -p general --mem=300G -t 9-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="{command}"')
        os.system(f'sbatch -p general --mem=300G -t 9-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="{command}"')