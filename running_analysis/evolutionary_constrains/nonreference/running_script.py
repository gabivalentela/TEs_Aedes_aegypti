import os

script = './frequency_analysis_nonreference.py'
list_dir = os.listdir('/work/users/g/a/gabivla/lab/SV_mosquitos/TE_paper_aedesaegytpi/nonreference_vcf/post_processing_nonreference_TEs/post_processing_overlap_df/max_min_updated_strategy')
for file in list_dir:
    if file.endswith('_1000.csv'):
        variation = file.split('_')[-1].split('.')[0]
        command = f'python {script} /work/users/g/a/gabivla/lab/SV_mosquitos/TE_paper_aedesaegytpi/nonreference_vcf/post_processing_nonreference_TEs/post_processing_overlap_df/max_min_updated_strategy/{file} {variation}'
        print(f'sbatch -p general --mem=20G -t 9-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="{command}"')
        os.system(f'sbatch -p general --mem=20G -t 9-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="{command}"')