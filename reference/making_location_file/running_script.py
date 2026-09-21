import os 

list_tolerance = [10, 50, 100]

for tolerance in list_tolerance:
    command = f"python making_gff3_file.py {tolerance}"
    print(f'sbatch -p general --mem=50G -t 9-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="{command}"')
    os.system(f'sbatch -p general --mem=50G -t 9-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="{command}"')