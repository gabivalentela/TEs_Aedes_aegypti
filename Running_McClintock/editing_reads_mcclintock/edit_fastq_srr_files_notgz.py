import os
import gzip
import re
    
#get srr files and edit the fastq content

def list_srr_files(path):
    list_of_srr_files = []
    list_files = os.listdir(path)
    for file in list_files:
        if 'SRR' in file:
            srr_file = file
            list_of_srr_files.append(f'{path}/{srr_file}')
            
    #to test first I am only running the 5 first genomes
    #n = 6
    #sub_list = list_of_srr_files[:n]
    return(list_of_srr_files)


def retrieve_srr_file(srr_dir):
    print(srr_dir)
    #files = os.listdir(srr_dir)
    fastq_dict = {}
    if '.fastq' in srr_dir:
        if '_1' in srr_dir:
            fastq_1 = f'{srr_dir}'
            fastq_dict.update({'1': fastq_1})
        elif '_2' in srr_dir:
            fastq_2 = f'{srr_dir}'
            fastq_dict.update({'2': fastq_2})
                
    return(fastq_dict)

def edit_fastq_files(fastq_dict):
    for type_fastq, fastq_path in fastq_dict.items():
        if '1' in type_fastq:
            name = re.search(r'.+/original_(\w+)\_\w\.fastq', fastq_path).group(1)
            output_file = f'/users/g/a/gabivla/extra_genomes_srr/{name}/{name}_1.fastq.gz'
            try:
                with open(fastq_path, 'r') as input_fastq, gzip.open(output_file, 'wt') as output_fastq:
                    for line in input_fastq:
                        if line.startswith('@') or line.startswith('+'):
                            line = line.rstrip()
                            pattern = re.search(r'(.+)\.(\w+)\s+\w+\s+length=(\d+)', line)
                            if pattern:
                                code = pattern.group(1)
                                position = pattern.group(2)
                                length = pattern.group(3)
                                output_fastq.write(f'{code}.{position} {position}/{type_fastq}\n')
                            else:
                                output_fastq.write(line + '\n')
                        else:
                            output_fastq.write(line)

            except Exception as e:
                print(f"Error processing {fastq_path}: {str(e)}")

        elif '2' in type_fastq:
            name = re.search(r'.+/original_(\w+)\_\w\.fastq', fastq_path).group(1)
            output_file = f'/users/g/a/gabivla/extra_genomes_srr/{name}/{name}_2.fastq.gz'
            
            try:
                with open(fastq_path, 'r') as input_fastq, gzip.open(output_file, 'wt') as output_fastq:
                    for line in input_fastq:
                        if line.startswith('@') or line.startswith('+'):
                            line = line.rstrip()
                            pattern = re.search(r'(.+)\.(\w+)\s+\w+\s+length=(\d+)', line)
                            if pattern:
                                code = pattern.group(1)
                                position = pattern.group(2)
                                length = pattern.group(3)
                                output_fastq.write(f'{code}.{position} {position}/{type_fastq}\n')
                            else:
                                output_fastq.write(line + '\n')
                        else:
                            output_fastq.write(line)

            except Exception as e:
                print(f"Error processing {fastq_path}: {str(e)}")

                        
def main():
    #path = '/work/users/g/a/gabivla/data_mosquito/data/fastq_files'
    #I am only using this because SRR11006666 had a problem but I can run other genomes by adjusting the list length on line 16
    
    #path = '/work/users/g/a/gabivla/data_mosquito/data_running_mcclintock/SRR_run_genome_test/multiple_srr_genomes/fastq_v4/SRR11006820'
    
    paths = ['/users/g/a/gabivla/extra_genomes_srr/SRR11006666', '/users/g/a/gabivla/extra_genomes_srr/SRR11006824', '/users/g/a/gabivla/extra_genomes_srr/SRR11006830', '/users/g/a/gabivla/extra_genomes_srr/SRR11006847']
    
    for path in paths:
        list_srr = list_srr_files(path)
        for folder in list_srr:
            fastq_files_dict = retrieve_srr_file(folder)
            print(fastq_files_dict)
            edit_fastq_files(fastq_files_dict)

if __name__ == "__main__":
    main()
