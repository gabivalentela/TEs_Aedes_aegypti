import os
import re

def find_fastq_folders(fastq_general_path):
    list_of_dir = os.listdir(fastq_general_path)
    list_of_path_folders = []

    for folder in list_of_dir:
        if folder == 'fastq_v6':
            folder_path = os.path.join(fastq_general_path, folder)
            list_files = os.listdir(folder_path)
            
            for subfolder in list_files:
                subfolder_path = os.path.join(folder_path, subfolder)
                list_of_path_folders.append(subfolder_path)
                
    return list_of_path_folders

'''def find_fastq_folders(fastq_general_path, skip_count=22, run_count=5):
    list_of_dir = os.listdir(fastq_general_path)
    list_of_path_folders = []
    count = 0
    
    # Skip the first `skip_count` folders
    for folder in list_of_dir:
        if folder == 'fastq_v4':
            folder_path = os.path.join(fastq_general_path, folder)
            list_files = os.listdir(folder_path)
            
            for subfolder in list_files:
                if count < skip_count:
                    count += 1
                    continue
                
                if count >= skip_count + run_count:
                    break
                
                subfolder_path = os.path.join(folder_path, subfolder)
                list_of_path_folders.append(subfolder_path)
                count += 1
                
    return list_of_path_folders'''

def get_fastq_files(fastq_path):
    fastq_files_list = os.listdir(fastq_path)
    list_fastq = []
    for fastq in fastq_files_list:
        fastq_files_path = f'{fastq_path}/{fastq}'
        list_fastq.append(fastq_files_path)
    return(list_fastq)

def run_mcclintock(fastq_path_general, files_fastq, mcclintock_path, genome_path, consesus_te_path):
    for fastq in files_fastq:
        #print(fastq)
        name = re.search(f'{fastq_path_general}/.+/(.+)/.+', fastq).group(1)
        if '_1.fastq.gz' in fastq:
            fastq_1 = fastq
        elif '_2.fastq.gz' in fastq:
            fastq_2 = fastq
            print(f'{name}')
            print(f'sbatch -p general -n 24 --mem=100G -t 9-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="python {mcclintock_path} -p 24 -r {genome_path} -c {consesus_te_path} -1 {fastq_1}  -2 {fastq_2} -o /work/users/g/a/gabivla/lab/SV_mosquitos/running_mcclintock/general_srr/output_{name} -m trimgalore,map_reads,temp2,ngs_te_mapper2"')
            os.system(f'sbatch -p general -n 24 --mem=100G -t 9-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="python {mcclintock_path} -p 24 -r {genome_path} -c {consesus_te_path} -1 {fastq_1}  -2 {fastq_2} -o /work/users/g/a/gabivla/lab/SV_mosquitos/running_mcclintock/general_srr/output_{name} -m trimgalore,map_reads,temp2,ngs_te_mapper2"')

def main():
    output_path_general = '/work/users/g/a/gabivla/lab/SV_mosquitos/running_mcclintock/'
    fastq_path_general = '/work/users/g/a/gabivla/data_mosquito/data_running_mcclintock/SRR_run_genome_test/multiple_srr_genomes'
    
    fastq_folder_path = find_fastq_folders(fastq_path_general)
    fastq_folder_path = fastq_folder_path[8:16]
    print(len(fastq_folder_path))
    
    mcclintock_path = '/work/users/g/a/gabivla/mcclintock/mcclintock.py'
    genome_path = '/work/users/g/a/gabivla/data_mosquito/data_running_mcclintock/reference_genome/VectorBase-54_AaegyptiLVP_AGWG_Genome.fasta'
    consesus_te_path = '/work/users/g/a/gabivla/data_mosquito/data_running_mcclintock/consensus_TE_seq/consensusTEs.fasta'
    for fastq_path in fastq_folder_path:
        fastq_files = get_fastq_files(fastq_path)
        run_mcclintock(fastq_path_general, fastq_files, mcclintock_path, genome_path, consesus_te_path)

if __name__ == "__main__":
    main()