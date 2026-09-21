import os

script = './frequency_analysis_reference.py'
list_dir = os.listdir('/work/users/g/a/gabivla/lab/SV_mosquitos/TE_paper_aedesaegytpi/reference_vcf/genotyping_reference_TEs/post_processing/')
for file in list_dir:
    if file.endswith('.csv'):
        variation = file.split('_')[-1].split('.')[0]
        command = f'python {script} /work/users/g/a/gabivla/lab/SV_mosquitos/TE_paper_aedesaegytpi/reference_vcf/genotyping_reference_TEs/post_processing/merged_results_5_2_either_filtered.csv'
        print(f'sbatch -p general --mem=20G -t 9-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="{command}"')
        os.system(f'sbatch -p general --mem=20G -t 9-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="{command}"')