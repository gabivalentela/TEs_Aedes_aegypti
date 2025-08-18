import os
import itertools

script_nonref = 'making_vcf_genotyping_tes.py'

path_chromsome = '../../../separate_chromosomes/making_list_of_chromosomes/chromosomes.txt'
chromosomes = open(path_chromsome, 'r').read().splitlines()

# Display the list of country pairs
for chromosome in chromosomes:
    if 'AaegL5' in chromosome:
        variation = 1000
        print(f'sbatch -p general -n 24 --mem=50G -t 9-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="python {script_nonref} {variation} {chromosome}"')
        os.system(f'sbatch -p general -n 24 --mem=50G -t 9-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="python {script_nonref} {variation} {chromosome}"')