import os
import re

def find_fastq_folders(fastq_general_path):
    list_of_dir = os.listdir(fastq_general_path)
    list_of_dir.sort()
    list_of_path_folders = []
    for folder in list_of_dir:
        folder_path = f'{fastq_general_path}/{folder}'
        list_of_path_folders.append(folder_path)
    
    return list_of_path_folders  

def get_fastq_files(fastq_path):
    fastq_files_list = os.listdir(fastq_path)
    list_fastq = []
    for fastq in fastq_files_list:
        fastq_files_path = f'{fastq_path}/{fastq}'
        list_fastq.append(fastq_files_path)
    return(list_fastq)

def run_mcclintock(fastq_path_general, files_fastq, mcclintock_path, genome_path, consesus_te_path, gff_path, taxonomy_path, output_path_general):
    fastq_dict = {}
    output_general_path = '/work/users/g/a/gabivla/lab/SV_mosquitos/new_library_pipeline_reference'

    for fastq in files_fastq:
        name_match = re.search(f'{fastq_path_general}/(.+)/.+', fastq)
        if name_match:
            name = name_match.group(1).replace('fastq_', '')
            if name not in fastq_dict:
                fastq_dict[name] = {'fastq_1': '', 'fastq_2': ''}

            if 'L001_L002_R1_001.fastq' in fastq:
                fastq_dict[name]['fastq_1'] = fastq
            elif 'L001_L002_R2_001.fastq' in fastq:
                fastq_dict[name]['fastq_2'] = fastq

    for name, fastqs in fastq_dict.items():
        fastq_1 = fastqs['fastq_1']
        fastq_2 = fastqs['fastq_2']
        if fastq_1 and fastq_2:
            print(f'{name}')    
            #print(f'sbatch -p general -n 24 --mem=200G -t 7-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="python {mcclintock_path} -p 24 -r {genome_path} -c {consesus_te_path} -1 {fastq_1}  -2 {fastq_2} -o /work/users/g/a/gabivla/lab/SV_mosquitos/new_library_pipeline_reference/outputs/output_{name} -m trimgalore,map_reads,temp2"')
            #os.system(f'sbatch -p general -n 24 --mem=200G -t 7-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="python {mcclintock_path} -p 24 -r {genome_path} -c {consesus_te_path} -1 {fastq_1}  -2 {fastq_2} -o /work/users/g/a/gabivla/lab/SV_mosquitos/new_library_pipeline_reference/outputs/output_{name} -m trimgalore,map_reads,temp2"')
            cmd = (
                f'sbatch -p general -n 24 --mem=100G -t 9-00:00:00 '
                f'--mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu '
                f'--wrap="python {mcclintock_path} -p 24 '
                f'-r {genome_path} '
                f'-c {consesus_te_path} '
                f'-g {gff_path} '
                f'-t {taxonomy_path} '
                f'-1 {fastq_1} '
                f'-2 {fastq_2} '
                f'-o {output_path_general}/output_{name} '
                f'-m trimgalore,map_reads,temp2"'
            )
            print(cmd)
            #os.system(cmd)

def main():
    #output_path_general = '/users/g/a/gabivla/new_library_pipeline_reference/outputs'
    output_path_general = '/work/users/g/a/gabivla/lab/SV_mosquitos/new_library_pipeline_reference_rerun/outputs'
    
    #fastq_path_general = '/work/users/g/a/gabivla/data_mosquito/data_running_mcclintock/flycross/fastq_flycross'
    fastq_path_general = '/users/g/a/gabivla/flycross_fastq'
    
    fastq_folder_path = find_fastq_folders(fastq_path_general)
    fastq_folder_path = fastq_folder_path[5:10]
    print(len(fastq_folder_path))
    mcclintock_path = '/users/g/a/gabivla/New_library_test/mcclintock/mcclintock.py'
    genome_path = '/work/users/g/a/gabivla/data_mosquito/data_running_mcclintock/reference_genome/VectorBase-54_AaegyptiLVP_AGWG_Genome.fasta'
    consesus_te_path = '/users/g/a/gabivla/New_library_test/Daron_2024_paper/aedes_aegypti_TElibrary/consensus.fasta'
    gff_path = '/work/users/g/a/gabivla/lab/SV_mosquitos/TEs_Aedes_aegypti/running_pipeline_new_library/new_run_reference_mcclintock/making_location_file/results/sorted_files/all_chromosomes_TEs_50_sorted.gff3'
    taxonomy_path = '/work/users/g/a/gabivla/lab/SV_mosquitos/TEs_Aedes_aegypti/running_pipeline_new_library/new_run_reference_mcclintock/making_location_file/results/sorted_files/all_chromosomes_taxonomy_50_sorted.tsv'
    for fastq_path in fastq_folder_path:
        fastq_files = get_fastq_files(fastq_path)
        run_mcclintock(fastq_path_general, fastq_files, mcclintock_path, genome_path, consesus_te_path, gff_path, taxonomy_path, output_path_general)

if __name__ == "__main__":
    main()