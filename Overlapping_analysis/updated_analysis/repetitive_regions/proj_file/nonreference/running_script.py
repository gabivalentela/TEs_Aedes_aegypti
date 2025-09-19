import os 

script_nonref = './overlap_analysis_nonreference.py'

list_of_countries = ['USA', 'Colombia', 'Brazil', 'Gabon', 'Senegal', 'Kenya']

for country in list_of_countries:
    #nonreference
    print(f'sbatch -p general -n 1 --mem=100G -t 7-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="python {script_nonref} {country}"')
    os.system(f'sbatch -p general -n 24 --mem=100G -t 7-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="python {script_nonref} {country}"')