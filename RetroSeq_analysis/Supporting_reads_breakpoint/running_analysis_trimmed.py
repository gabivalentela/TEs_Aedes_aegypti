import os

# Run module load bedtools before running it

def find_information(general_path, breakpoint_value):
    information_dict = {}
    list_files = os.listdir(general_path)

    #list_of_genomes = ["output_SRR6768001", "output_SRR6768012", "output_SRR6768020"]
    #list_of_genomes = ["output_MALE_6-M6_S28", "output_FEMALE_12-F12_S4", "output_MALE_5-M5_S27", "output_FEMALE_21-F21_S14", "output_FEMALE_8-F8_S21"]
    #list_of_genomes = ["output_FEMALE_1-F1_S1", "output_FEMALE_3-F3_S16", "output_FEMALE_22-F22_S15"]
    #list_of_genomes = ["output_SRR6768003", "output_SRR6768006", "output_SRR6768007", "output_SRR6768008", "output_SRR6768010", "output_SRR6768011", "output_SRR6768011", "output_SRR6768018", "output_SRR6768019", "output_SRR6768023"]
    #list_of_genomes = ['output_SRR6768020']
    #list_of_genomes = ['output_SRR11006851', 'output_SRR11006763']
    #list_of_genomes = ['output_SRR6768004', 'output_SRR6768020', 'output_SRR6768026', 'output_SRR6768027']
    #list_of_genomes = ['output_SRR6768022']
    #list_of_genomes = ['output_SRR6768005']
    #list_of_genomes = ['output_SRR11006836', 'output_SRR11006760']
    #list_of_genomes = ['output_SRR11006685', 'output_SRR11006829', 'output_SRR11006849']
    #list_of_genomes = ['output_SRR11006674', 'output_SRR11006681', 'output_SRR11006757', 'output_SRR11006828', 'output_SRR11006758', 'output_SRR11006752']
    #list_of_genomes = ['output_SRR11006827', 'output_SRR11006668', 'output_SRR11006770', 'output_SRR11006826', 'output_SRR11006751']
    #list_of_genomes = ["output_SRR11006672", "output_SRR11006839", "output_SRR11006852", "output_SRR11006825", "output_SRR11006759"]
    #list_of_genomes = ['output_SRR11006680', 'output_SRR11006753', 'output_SRR11006764', 'output_SRR11006821', 'output_SRR11006830', 'output_SRR11006843', 'output_SRR11006846', 'output_SRR11006854']
    #list_of_genomes = ['output_SRR11006667', 'output_SRR11006673', 'output_SRR11006683', 'output_SRR11006762', 'output_SRR11006765', 'output_SRR11006766', 'output_SRR11006767', 'output_SRR11006831', 'output_SRR11006834']
    list_of_genomes = ['output_SRR6768005']
    
    for directory in list_files:
        if 'output_' in directory and directory in list_of_genomes:
            information_dict[directory] = {}
            list_files_dir = os.listdir(f'{general_path}/{directory}')

            sample_name = None
            ref_name = None

            for file in list_files_dir:
                if 'VectorBase' in file:
                    ref_name = file
                elif '_75' in file:
                    sample_name = file

            if not sample_name or not ref_name:
                print(f"Warning: Missing expected files in directory {directory}")
                continue

            bam_file = f'{general_path}/{directory}/{sample_name}/intermediate/mapped_reads/{sample_name}.sorted.bam'
            reference_fasta = '/work/users/g/a/gabivla/data_mosquito/data_running_mcclintock/reference_genome/VectorBase-54_AaegyptiLVP_AGWG_Genome.fasta'
            retroseq_out = f'{general_path}/{directory}/{sample_name}/results/retroseq/unfiltered/{sample_name}.call'
            out_dir = f'breakpoint_{breakpoint_value}'
            #out_dir = f'{general_path}/{directory}/{sample_name}/results/retroseq/edited_results/breakpoint_{breakpoint_value}'
            retroseq_discovery = f'{general_path}/{directory}/{sample_name}/results/retroseq/unfiltered/{sample_name}.discovery'
            
            # Store information in the dictionary
            information_dict[directory] = {
                'sample_name': sample_name,
                'ref_name': ref_name,
                'reference_fasta': reference_fasta,
                'retroseq_out': retroseq_out,
                'out_dir': out_dir,
                'bam_file': bam_file,
                'retroseq_discovery': retroseq_discovery
            }
    
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
#general_path = '/work/users/g/a/gabivla/lab/SV_mosquitos/running_mcclintock/US_genomes'
#general_path = '/work/users/g/a/gabivla/lab/SV_mosquitos/trimming_reads_fastqs/new_test_fastp/worst_inserts_genomes'
#general_path=  '/users/g/a/gabivla/running_mcclintock_trimmed_reads/flycross'
#general_path = '/users/g/a/gabivla/running_mcclintock_trimmed_reads/extra_genomes'
general_path = '/users/g/a/gabivla/running_mcclintock_trimmed_reads/US_genomes'
#general_path = '/users/g/a/gabivla/running_mcclintock_trimmed_reads/general_srr'
#general_path = '/work/users/g/a/gabivla/lab/SV_mosquitos/running_mcclintock_trimmed/general_srr'

breakpoint_values = [1, 6]

for breakpoint_value in breakpoint_values:
    # Retrieve dictionary with information from directories
    dict_information = find_information(general_path, breakpoint_value)

    #print(dict_information)
    # Get only the first 10 genomes
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