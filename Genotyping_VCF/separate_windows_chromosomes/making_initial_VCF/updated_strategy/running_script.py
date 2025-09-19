import os
import itertools

script_ref = 'making_vcf_file_for_each_window.py'

path_chromsome = '../../../separate_chromosomes/making_list_of_chromosomes/chromosomes.txt'

chromosomes = open(path_chromsome, 'r').read().splitlines()

# Display the list of country pairs
for chromosome in chromosomes:
    if 'AaegL5' in chromosome:
        variation = 1000
        print(f'sbatch -p general -n 1 --cpus-per-task=150 --mem=200G -t 10-00:00:00 --mail-type=END,FAIL --mail-user=gabivla@email.unc.edu --wrap="python {script_ref} {variation} {chromosome}"')
        os.system(f'sbatch -p general -n 1 --cpus-per-task=150 --mem=200G -t 10-00:00:00 --mail-type=END,FAIL --mail-user=gabivla@email.unc.edu --wrap="python {script_ref} {variation} {chromosome}"')
        
    else:
        pass

