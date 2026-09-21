import os
script = './comparing_dfs.py'
list_files = os.listdir('../')
ref_file = '../genes_overlap_tes_ref_no_duplicates.csv'
nonref_file = '../genes_overlap_tes_nonref_no_duplicates_1000.csv'
command = f'sbatch -p general -n 1 --mem=10G -t 7-00:00:00 --mail-type=BEGIN,END,FAIL --mail-user=gabivla@email.unc.edu --wrap="python {script} {ref_file} {nonref_file}"'
print(f'{command}')
os.system(command)