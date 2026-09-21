import os

script_path = './optmized_version_vcf_teforest.py'
chromosomes = ["AaegL5_1", "AaegL5_2", "AaegL5_3"]
tolerance = ['100', '500', '1000']
for chromosome in chromosomes:
    for tol in tolerance:
        command = f"python {script_path} {chromosome} {tol}"
        print(f'sbatch -p general --mem=150G -t 9-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="{command}"')    
        os.system(f'sbatch -p general --mem=150G -t 9-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="{command}"')