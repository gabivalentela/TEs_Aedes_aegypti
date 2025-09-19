import os

def running_trimming_script(list_genomes):
    for genome in list_genomes:
        fastq1 = genome['fastq1']
        fastq2 = genome['fastq2']
        new_name1 = genome['new_name_1']
        new_name2 = genome['new_name_2']

        # value 75
        print(f'fastp -w 24 -i {fastq1} -I {fastq2} -o {new_name1}_75.fq -O {new_name2}_75.fq --max_len1 75 --max_len2 75')
        os.system(f'sbatch -p general -n 24 --mem=200G -t 7-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="fastp -w 24 -i {fastq1} -I {fastq2} -o {new_name1}_75.fq -O {new_name2}_75.fq --max_len1 75 --max_len2 75"')

def get_list_of_genomes(general_path, output_basic):
    
    list_of_genomes = []
    
    directories = os.listdir(general_path)
    
    for directory in directories:
        # Ensure output directory exists
        output_dir = f'{output_basic}/{directory}'
        os.makedirs(output_dir, exist_ok=True)
        
        directory_genome = {}  # Initialize an empty dictionary for each directory
        files = os.listdir(f'{general_path}/{directory}')
        
        for file in files:
            if 'R1' in file:
                directory_genome['fastq1'] = f'{general_path}/{directory}/{file}'
                directory_genome['new_name_1'] = f'{output_dir}/{file}'
                
            elif 'R2' in file:
                directory_genome['fastq2'] = f'{general_path}/{directory}/{file}'
                directory_genome['new_name_2'] = f'{output_dir}/{file}'
        
        list_of_genomes.append(directory_genome)
    
    return list_of_genomes

# module load fastp
##list of genome information

general_path = '/users/g/a/gabivla/flycross_fastq'
output_basic = '/users/g/a/gabivla/trimed_reads/flycross'
list_of_genomes = get_list_of_genomes(general_path, output_basic)

running_trimming_script(list_of_genomes)