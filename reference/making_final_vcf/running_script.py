import os 

script = './optimized_making_vcf.py'
chromosomes = ['AaegL5_1', 'AaegL5_2', 'AaegL5_3']
for chrom in chromosomes:
    command = f'python {script} {chrom}'
    print(f'sbatch -p general --mem=150G -t 9-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="{command}"')
    os.system(f'sbatch -p general --mem=150G -t 9-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="{command}"')