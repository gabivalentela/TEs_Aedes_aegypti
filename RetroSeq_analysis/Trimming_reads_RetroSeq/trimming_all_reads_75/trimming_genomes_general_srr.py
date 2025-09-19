import os

def running_trimming_script(list_genomes):
    for genome in list_genomes:
        fastq1 = genome.get('fastq1')
        fastq2 = genome.get('fastq2')
        new_name1 = genome.get('new_name_1')
        new_name2 = genome.get('new_name_2')

        # Check if both fastq1 and fastq2 files are present
        if not fastq1 or not fastq2:
            print("Warning: One or both SRR files are missing for a genome.")
            print(f"fastq1: {fastq1}, fastq2: {fastq2}")
            continue  # Skip to the next genome if files are missing

        # value 75
        print(f'fastp -w 24 -i {fastq1} -I {fastq2} -o {new_name1}_75.fq -O {new_name2}_75.fq --max_len1 75 --max_len2 75')
        os.system(f'sbatch -p general -n 24 --mem=20G -t 7-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="fastp -w 24 -i {fastq1} -I {fastq2} -o {new_name1}_75.fq -O {new_name2}_75.fq --max_len1 75 --max_len2 75"')

def get_list_of_genomes(general_path, output_basic):
    list_of_genomes = []

    initial_dir = os.listdir(general_path)
    for initial in initial_dir:
        if initial == 'genome_fasta_file' or initial == 'fastq_SRR6666_2':
            pass
        else:
            directories = os.listdir(f'{general_path}/{initial}')
            
            for directory in directories:
                # Ensure output directory exists
                output_dir = f'{output_basic}/{directory}'
                os.makedirs(output_dir, exist_ok=True)
                
                directory_genome = {}  # Initialize an empty dictionary for each directory
                files = os.listdir(f'{general_path}/{initial}/{directory}')
                
                for file in files:
                    if '_1' in file:
                        directory_genome['fastq1'] = f'{general_path}/{initial}/{directory}/{file}'
                        directory_genome['new_name_1'] = f'{output_dir}/{file}'
                        
                    elif '_2' in file:
                        directory_genome['fastq2'] = f'{general_path}/{initial}/{directory}/{file}'
                        directory_genome['new_name_2'] = f'{output_dir}/{file}'

                # Check if both fastq1 and fastq2 keys exist in dictionary
                if 'fastq1' not in directory_genome or 'fastq2' not in directory_genome:
                    print(f"Warning: Missing SRR files in {directory}. Expected both _1 and _2 files.")
                
                list_of_genomes.append(directory_genome)
    
    return list_of_genomes

# module load fastp
## List of genome information

general_path = '/work/users/g/a/gabivla/data_mosquito/data_running_mcclintock/SRR_run_genome_test/multiple_srr_genomes'
output_basic = '/users/g/a/gabivla/trimed_reads/general_srr'
list_of_genomes = get_list_of_genomes(general_path, output_basic)

running_trimming_script(list_of_genomes)
