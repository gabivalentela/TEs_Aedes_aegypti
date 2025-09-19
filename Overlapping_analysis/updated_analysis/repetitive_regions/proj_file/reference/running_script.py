import os 

script_ref ='./overlap_TEs_exons_ref.py'

list_of_countries = ['USA', 'Colombia', 'Brazil', 'Gabon', 'Senegal', 'Kenya']

for country in list_of_countries:
    #reference
    print(f'sbatch -p general -n 1 --mem=100G -t 7-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="python {script_ref} {country}"')
    os.system(f'sbatch -p general -n 24 --mem=100G -t 7-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="python {script_ref} {country}"')
    