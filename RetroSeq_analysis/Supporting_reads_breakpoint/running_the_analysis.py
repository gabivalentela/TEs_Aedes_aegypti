import os

def find_information(general_path, breakpoint_value):
    information_dict = {}
    list_files = os.listdir(general_path)
    
    for directory in list_files:
        information_dict[directory] = {}

        list_files_dir = os.listdir(f'{general_path}/{directory}')
        #sample_name = list_files_dir[0]
        #ref_name = list_files_dir[1]
        
        sample_name = list_files_dir[-1]
        ref_name = list_files_dir[0]
        
        bam_file = f'{general_path}/{directory}/{sample_name}/intermediate/mapped_reads/{sample_name}.sorted.bam'
        reference_fasta = '/work/users/g/a/gabivla/data_mosquito/data_running_mcclintock/reference_genome/VectorBase-54_AaegyptiLVP_AGWG_Genome.fasta'
        retroseq_out = f'{general_path}/{directory}/{sample_name}/results/retroseq/unfiltered/{sample_name}.call'
        out_dir = f'{general_path}/{directory}/{sample_name}/results/retroseq/edited_results/breakpoint_{breakpoint_value}'
        retroseq_discovery = f'{general_path}/{directory}/{sample_name}/results/retroseq/unfiltered/{sample_name}.discovery'

        # Store information in the dictionary
        information_dict[directory]['sample_name'] = sample_name
        information_dict[directory]['ref_name'] = ref_name
        information_dict[directory]['reference_fasta'] = reference_fasta
        information_dict[directory]['retroseq_out'] = retroseq_out
        information_dict[directory]['out_dir'] = out_dir
        information_dict[directory]['bam_file'] = bam_file
        information_dict[directory]['retroseq_discovery'] = retroseq_discovery
        
    return information_dict

def create_directory_if_not_exists(directory_path):
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)
        print(f"Directory {directory_path} created.")
    else:
        print(f"Directory {directory_path} already exists.")

#general_path = '/work/users/g/a/gabivla/lab/SV_mosquitos/running_mcclintock/flycross'
#general_path = '/work/users/g/a/gabivla/lab/SV_mosquitos/running_mcclintock/general_srr'
#general_path = '/work/users/g/a/gabivla/lab/SV_mosquitos/running_mcclintock/extra_genomes'
general_path = '/work/users/g/a/gabivla/lab/SV_mosquitos/running_mcclintock/US_genomes'

breakpoint_values = [1, 6]

for breakpoint_value in breakpoint_values:
    # Retrieve dictionary with information from directories
    dict_information = find_information(general_path, breakpoint_value)
    first_ten_genomes = {key: dict_information[key] for key in dict_information.keys()}

    for directory, info in first_ten_genomes.items():
        create_directory_if_not_exists(info['out_dir'])

        if breakpoint_value == 1:
            # Paths and parameters
            vcf_options = ["option1", "option2"]
            config_path = "/nas/longleaf/home/gabivla/lab/SV_mosquitos/breakpoint_analysis/breakpoint_files_generating/retroseq_post_scripts/retroseq_post_1.py"
            mcc_path = "/nas/longleaf/home/gabivla/lab/SV_mosquitos/retroseq_checking_TEs/comparing_call_vcf/making_vcf_file"
                
            retroseq_out = info['retroseq_out']
            reference_fasta = info['reference_fasta']
            out_dir = info['out_dir']
            ref_name = info['ref_name']
            sample_name = info['sample_name']
            bam_file = info['bam_file']
            retroseq_discovery = info['retroseq_discovery']
            
            print(f"Processing directory: {directory}")
            print(f'sbatch -p general -n 24 --mem=20G -t 7-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="python one_job_each_genome.py {retroseq_out} {reference_fasta} {out_dir} {ref_name} {sample_name} {vcf_options} {config_path} {mcc_path} {bam_file} {retroseq_discovery} {breakpoint_value}"')
            os.system(f'sbatch -p general -n 24 --mem=20G -t 7-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="python one_job_each_genome.py {retroseq_out} {reference_fasta} {out_dir} {ref_name} {sample_name} {vcf_options} {config_path} {mcc_path} {bam_file} {retroseq_discovery} {breakpoint_value}"')

        elif breakpoint_value == 6:
            # Paths and parameters
            vcf_options = ["option1", "option2"]
            config_path = "/nas/longleaf/home/gabivla/lab/SV_mosquitos/breakpoint_analysis/breakpoint_files_generating/retroseq_post_scripts/retroseq_post_6.py"
            mcc_path = "/nas/longleaf/home/gabivla/lab/SV_mosquitos/retroseq_checking_TEs/comparing_call_vcf/making_vcf_file"
                
            retroseq_out = info['retroseq_out']
            reference_fasta = info['reference_fasta']
            out_dir = info['out_dir']
            ref_name = info['ref_name']
            sample_name = info['sample_name']
            bam_file = info['bam_file']
            retroseq_discovery = info['retroseq_discovery']
            
            print(f"Processing directory: {directory}")
            print(f'sbatch -p general -n 24 --mem=20G -t 7-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="python one_job_each_genome.py {retroseq_out} {reference_fasta} {out_dir} {ref_name} {sample_name} {vcf_options} {config_path} {mcc_path} {bam_file} {retroseq_discovery} {breakpoint_value}"')
            os.system(f'sbatch -p general -n 24 --mem=20G -t 7-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="python one_job_each_genome.py {retroseq_out} {reference_fasta} {out_dir} {ref_name} {sample_name} {vcf_options} {config_path} {mcc_path} {bam_file} {retroseq_discovery} {breakpoint_value}"')
