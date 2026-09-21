import os
import pandas as pd

vcf_reference = '/work/users/g/a/gabivla/lab/SV_mosquitos/TEs_Aedes_aegypti/running_pipeline_new_library/new_run_reference_mcclintock/making_final_vcf/merge_chromosomes/combined_reference_tes.vcf'
df = pd.read_csv(vcf_reference, sep='\t', comment='#', nrows=5)
# make columns into list
columns = df.columns.tolist()
remove_col = ['chrom', 'start', 'end', 'strand', 'TE_name', 'score']
columns = [col for col in columns if col not in remove_col]

chromosomes = ['AaegL5_1', 'AaegL5_2', 'AaegL5_3']
#chromosomes = ['AaegL5_1']

n_cores = 6 
#list_supporting_reads = [1, 2]
number_supporting_reads_presence = 2
number_supporting_reads_absence = 5
support_side = ['either']
    
for chrom in chromosomes:
    for genome in columns:
        for sides_support_TEs in support_side:
            command = f'python ./genotyping_dfs_reference_tes_more_strict.py {chrom} {genome} {number_supporting_reads_presence} {number_supporting_reads_absence} {sides_support_TEs}'
            sbatch_cmd = (
                        f'sbatch -p general --mem=10G -t 9-00:00:00 '
                        f'--cpus-per-task={n_cores} '
                    f'--mail-type=FAIL --mail-user=gabivla@email.unc.edu '
                    f'--wrap="{command}"'
                )
                
            print(sbatch_cmd)
            os.system(sbatch_cmd)