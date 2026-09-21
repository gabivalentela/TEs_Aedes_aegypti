import os
import sys
import re
import pandas as pd

def getting_temp2_output(general_path):
    temp2_files = {}
    if not os.path.exists(general_path):
        return temp2_files
    for output_dir in os.listdir(general_path):
        updated_path = os.path.join(general_path, output_dir)
        genome_code = re.search('output_(.+)', output_dir).group(1)
        for inner_file in os.listdir(updated_path):
            if inner_file.endswith('1'):
                genome_code = inner_file
        new_path = os.path.join(updated_path, genome_code, 'results', 'temp2')
        if os.path.exists(new_path):
            for temp2_file in os.listdir(new_path):
                if temp2_file.endswith('_nonredundant.bed'):
                    temp2_files[genome_code] = os.path.join(new_path, temp2_file)
    return temp2_files

MERGE_KEYS = ['chrom', 'start', 'end', 'strand', 'TE_name', 'score']

def load_sample(genome_code, path, chrom):
    """Load and preprocess a single sample's BED file."""
    df = pd.read_csv(
        path, sep='\t', header=None, skiprows=1,
        usecols=[0, 1, 2, 3, 4, 5],                   # avoid reading extra cols
        dtype={'chrom': 'category', 'strand': 'category'},
        names=['chrom', 'start', 'end', 'name', 'score', 'strand']
    )
    df = df[df['chrom'] == chrom]
    # parse name column once
    name_split = df['name'].str.split('|', n=2, expand=True)
    df['TE_name'] = name_split[0]
    df['type']    = name_split[1]
    df = df.loc[df['type'] == 'reference', MERGE_KEYS].copy()
    df[genome_code] = '1/1'
    return df

def making_vcf(temp2_files, chrom):
    # load all samples
    dfs = [load_sample(gc, path, chrom) for gc, path in temp2_files.items()]

    # single concat + groupby instead of iterative merges
    combined = pd.concat(dfs, ignore_index=True)

    genome_cols = list(temp2_files.keys())
    merged_df = (
        combined
        .groupby(MERGE_KEYS, observed=True)[genome_cols]
        .first()                        # each sample only has one row per locus
        .reset_index()
        .fillna('0/0')
        .sort_values(['chrom', 'start', 'end'])
        .reset_index(drop=True)
    )
    return merged_df

chrom = sys.argv[1]
general_path   = "/work/users/g/a/gabivla/lab/SV_mosquitos/new_library_pipeline_reference_rerun/outputs/" 
general_path_2 = "/users/g/a/gabivla/new_library_pipeline_reference/outputs/"

temp2_files = getting_temp2_output(general_path)
temp2_files.update(getting_temp2_output(general_path_2))

merged_df = making_vcf(temp2_files, chrom)

merged_df.to_csv(f'output_complete_{chrom}.vcf', sep='\t', index=False)