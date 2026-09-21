import os 

script_ref ='./overlap_TEs_exons_ref.py'
#script_nonref = '/nas/longleaf/home/gabivla/lab/SV_mosquitos/final_scripts/Overlapping_analysis/updated_analysis/nonreference/overlap_analysis_nonreference.py'

#list_of_countries = ['USA', 'Colombia', 'Brazil', 'Gabon', 'Senegal', 'Kenya']
list_of_countries = ['USA', 'Colombia']

for country in list_of_countries:
    #reference
    print(f'sbatch -p general -n 1 --mem=10G -t 7-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="python {script_ref} {country}"')
    os.system(f'sbatch -p general -n 24 --mem=10G -t 7-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="python {script_ref} {country}"')
    