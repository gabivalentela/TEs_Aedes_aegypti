import os 

script_ref ='./overlap_analysis_nonreference.py'
#script_nonref = '/nas/longleaf/home/gabivla/lab/SV_mosquitos/final_scripts/Overlapping_analysis/updated_analysis/nonreference/overlap_analysis_nonreference.py'
list_files = os.listdir('/work/users/g/a/gabivla/lab/SV_mosquitos/TE_paper_aedesaegytpi/nonreference_vcf/post_processing_nonreference_TEs/post_processing_overlap_df/max_min_updated_strategy/')
list_files = [file for file in list_files if file.endswith('_1000.csv')]
list_of_countries = ['USA', 'Colombia', 'Brazil', 'Gabon', 'Senegal', 'Kenya']
for file in list_files:
    if file.endswith('_1000.csv'):
        for country in list_of_countries:
            variation = file.split('_')[-1].split('.')[0]
            command = f'python {script_ref} {country} /work/users/g/a/gabivla/lab/SV_mosquitos/TE_paper_aedesaegytpi/nonreference_vcf/post_processing_nonreference_TEs/post_processing_overlap_df/max_min_updated_strategy/{file} {variation}'
            print(f'sbatch -p general -n 1 --mem=10G -t 7-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="{command}"')
            os.system(f'sbatch -p general -n 24 --mem=10G -t 7-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="{command}"')