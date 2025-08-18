import os
import re

def find_fastq_folders(fastq_general_path):
    list_of_dir = os.listdir(fastq_general_path)
    list_of_path_folders = []
    for folder in list_of_dir:
        if folder.startswith('genomes'):
            folder_path = f'{fastq_general_path}/{folder}'
            list_files = os.listdir(folder_path)
            for subfolder in list_files:
                subfolder_path = f'{fastq_general_path}/{folder}/{subfolder}'
                list_of_path_folders.append(subfolder_path)
    return list_of_path_folders[26:]  

def get_fastq_files(fastq_path):
    fastq_files_list = os.listdir(fastq_path)
    list_fastq = []
    for fastq in fastq_files_list:
        fastq_files_path = f'{fastq_path}/{fastq}'
        list_fastq.append(fastq_files_path)
    return list_fastq

def run_mcclintock(fastq_path_general, files_fastq, mcclintock_path, genome_path, consesus_te_path):
    for fastq in files_fastq:
        print(fastq)
        name = re.search(f'{fastq_path_general}/.+/(.+)/.+', fastq).group(1)
        if '_1.fastq.gz' in fastq:
            fastq_1 = fastq
        elif '_2.fastq.gz' in fastq:
            fastq_2 = fastq
            print(f'{name}')
            print(f'sbatch -p general -n 24 --mem=100G -t 7-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="python {mcclintock_path} -p 24 -r {genome_path} -c {consesus_te_path} -1 {fastq_1}  -2 {fastq_2} -o /work/users/g/a/gabivla/lab/SV_mosquitos/running_mcclintock/US_genomes/output_{name} -m trimgalore,map_reads,temp2,ngs_te_mapper2"')
            os.system(f'sbatch -p general -n 24 --mem=100G -t 7-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="python {mcclintock_path} -p 24 -r {genome_path} -c {consesus_te_path} -1 {fastq_1}  -2 {fastq_2} -o /work/users/g/a/gabivla/lab/SV_mosquitos/running_mcclintock/US_genomes/output_{name} -m trimgalore,map_reads,temp2,ngs_te_mapper2"')
            
def main():
    output_path_general = '/work/users/g/a/gabivla/lab/SV_mosquitos/running_mcclintock/'
    fastq_path_general = '/work/users/g/a/gabivla/data_mosquito/data_running_mcclintock/SRR_genomes_US'
    
    fastq_folder_path = find_fastq_folders(fastq_path_general)
    print(len(fastq_folder_path))
    mcclintock_path = '/work/users/g/a/gabivla/mcclintock/mcclintock.py'
    genome_path = '/work/users/g/a/gabivla/data_mosquito/data_running_mcclintock/reference_genome/VectorBase-54_AaegyptiLVP_AGWG_Genome.fasta'
    consesus_te_path = '/work/users/g/a/gabivla/data_mosquito/data_running_mcclintock/consensus_TE_seq/consensusTEs.fasta'
    
    for fastq_path in fastq_folder_path:
        fastq_files = get_fastq_files(fastq_path)
        run_mcclintock(fastq_path_general, fastq_files, mcclintock_path, genome_path, consesus_te_path)

if __name__ == "__main__":
    main()
