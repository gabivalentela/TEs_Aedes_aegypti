import os
import itertools

script_nonref = './updated_making_VCF_retroseq.py'

chromosomes = ['AaegL5_1', 'AaegL5_2', 'AaegL5_3']
variation_values = [100, 500, 1000]
# Display the list of country pairs
for chromosome in chromosomes:
    if 'AaegL5' in chromosome:
        for variation in variation_values:
            print(f'sbatch -p general -n 24 --mem=50G -t 9-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="python {script_nonref} {variation} {chromosome}"')
            os.system(f'sbatch -p general -n 24 --mem=50G -t 9-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="python {script_nonref} {variation} {chromosome}"')