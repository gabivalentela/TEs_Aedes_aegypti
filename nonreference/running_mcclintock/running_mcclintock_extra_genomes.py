import os
import re

def find_fastq_folders(fastq_general_path):
    list_of_dir = os.listdir(fastq_general_path)
    list_of_path_folders = []
    for folder in list_of_dir:
        if folder.startswith('SRR'):
            folder_path = f'{fastq_general_path}/{folder}'
            list_files = os.listdir(folder_path)
            list_of_path_folders.append(folder_path)
    return list_of_path_folders

def get_fastq_files(fastq_path):
    fastq_files_list = os.listdir(fastq_path)
    list_fastq = []
    for fastq in fastq_files_list:
        fastq_files_path = f'{fastq_path}/{fastq}'
        list_fastq.append(fastq_files_path)
    return list_fastq
            
def run_mcclintock(fastq_path_general, files_fastq, mcclintock_path, genome_path, consesus_te_path):
    fastq_dict = {}

    for fastq in files_fastq:
        name_match = re.search(f'{fastq_path_general}/(.+)/.+', fastq)
        if name_match:
            name = name_match.group(1).replace('fastq_', '')
            if name not in fastq_dict:
                fastq_dict[name] = {'fastq_1': '', 'fastq_2': ''}

            if '_1.fastq.gz' in fastq:
                fastq_dict[name]['fastq_1'] = fastq
            elif '_2.fastq.gz' in fastq:
                fastq_dict[name]['fastq_2'] = fastq

    for name, fastqs in fastq_dict.items():
        fastq_1 = fastqs['fastq_1']
        fastq_2 = fastqs['fastq_2']
        if fastq_1 and fastq_2:
            print(f'{name}')
            print(f'sbatch -p general -n 24 --mem=200G -t 9-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="python {mcclintock_path} -p 24 -r {genome_path} -c {consesus_te_path} -1 {fastq_1}  -2 {fastq_2} -o /users/g/a/gabivla/new_library_pipeline_nonreference/outputs/output_{name} -m trimgalore,map_reads,retroseq,teflon,temp2"')
            os.system(f'sbatch -p general -n 24 --mem=200G -t 7-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="python {mcclintock_path} -p 24 -r {genome_path} -c {consesus_te_path} -1 {fastq_1}  -2 {fastq_2} -o /users/g/a/gabivla/new_library_pipeline_nonreference/outputs/output_{name} -m trimgalore,map_reads,retroseq,teflon,temp2"')
            
            
def main():
    output_path_general = '/users/g/a/gabivla/new_library_pipeline_nonreference/outputs'
    fastq_path_general = '/users/g/a/gabivla/trimed_reads/extra_genomes'
    
    fastq_folder_path = find_fastq_folders(fastq_path_general)
    mcclintock_path = '/users/g/a/gabivla/New_library_test/mcclintock/mcclintock.py'
    genome_path = '/work/users/g/a/gabivla/data_mosquito/data_running_mcclintock/reference_genome/VectorBase-54_AaegyptiLVP_AGWG_Genome.fasta'
    consesus_te_path = '/users/g/a/gabivla/New_library_test/Daron_2024_paper/aedes_aegypti_TElibrary/modified_consensus_updated.fasta'
    
    for fastq_path in fastq_folder_path:
        fastq_files = get_fastq_files(fastq_path)
        run_mcclintock(fastq_path_general, fastq_files, mcclintock_path, genome_path, consesus_te_path)
    
if __name__ == "__main__":
    main()