import pandas as pd
import polars as pl
from itertools import combinations
import os
import re
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

def count_TE_occurence(general_df):
    """
    Counts the occurrences of TEs and identifies the top 5 most frequent ones.

    Args:
        general_df (pl.DataFrame): Polars DataFrame containing TE data.

    Returns:
        pl.DataFrame: DataFrame with the top 5 TEs sorted by their occurrence count.
    """
    general_df = general_df.with_columns(pl.col("TE_name").str.to_uppercase())

    df_te = general_df.with_columns(pl.col("TE_name").str.extract(r'RND_\d+_FAMILY_\d+_([A-Z]+(?:_HELITRON)?)', 1).alias("TE_name"))

    grouped = df_te.group_by('TE_name').len()
    sorted_grouped = grouped.sort('len',  descending=True)
    
    return sorted_grouped

def calculate_frequency_TE(df, country):
    """
    For each TE insertion (row), calculate the derived allele frequency
    across all sample columns.
    Genotype encoding:
        1/1        -> 2 derived alleles
        0/1 / 1/0  -> 1 derived allele  (heterozygous, phasing-agnostic)
        0/0        -> 0 derived alleles
        anything else (missing, .'/.' etc.) -> excluded from denominator
    """
    base_cols = ['chromosome', 'position_start', 'position_end', 'TE_name', 'strand']
    meta_cols     = [c for c in base_cols if c in df.columns] + ["te_class", "te_size", "size_class"]
    genotype_cols = [c for c in df.columns if c not in meta_cols]
    df = df.copy()
    # calculate the frequencies
    df['count_1_1'] = (df[genotype_cols] == "1/1").sum(axis=1)
    df['count_0_1'] = (df[genotype_cols] == "0/1").sum(axis=1)
    df['count_0_0'] = (df[genotype_cols] == "0/0").sum(axis=1)
    df['called'] = (df['count_1_1'] + df['count_0_1'] + df['count_0_0'])*2
    df[f'frequency_{country}'] = (2 * df['count_1_1'] + df['count_0_1']) / (df['called'])
    df['frequency_0'] = (2 * df['count_0_0'] + df['count_0_1']) / (df['called'])
    # remove rows were all rows are missing (no called genotypes)
    df = df[df['called'] > 0]
    # remove rows where all genotypes are fixed (either all 1/1 or all 0/0)
    df = df[df[f'frequency_{country}'] < 1]
    df = df[df['frequency_0'] < 1]
    # drop the intermediate count columns
    df = df.drop(columns=['count_1_1', 'count_0_1', 'count_0_0', 'called', 'frequency_0'])

    return df

def ir_gene_regions(df, IR_genes, country):
    frequency_column = df[f'frequency_{country}']
    df = df.rename(columns={"chromosome": "chrom", "position_start": "start", "position_end": "end"})

    df[['start', 'end']] = df[['start', 'end']].apply(pd.to_numeric, errors='coerce')
    IR_genes[['start', 'end']] = IR_genes[['start', 'end']].apply(pd.to_numeric, errors='coerce')

    merged = df.merge(IR_genes, on="chrom", how="left")

    merged = merged[
        (merged["start_x"] >= merged["start_y"]) & 
        (merged["end_x"] <= merged["end_y"])
    ]
    
    merged = merged[["chrom", "start_x", "end_x", "TE_name", "other name", "##gene ID", f'frequency_{country}']]
    merged = merged.rename(columns={"start_x": "start", "end_x": "end", "other name": "gene_name", "##gene ID": "gene_ID"})
    
    return merged

'''nonreference_path = '/work/users/g/a/gabivla/lab/SV_mosquitos/TE_paper_aedesaegytpi/nonreference_vcf/post_processing_nonreference_TEs/post_processing_overlap_df/max_min_updated_strategy/final_df_overlap_retroseq_teforest_nonreference_1000.csv'
df_ref = pd.read_csv(nonreference_path)

path_becca = '/nas/longleaf/home/gabivla/lab/SV_mosquitos/IR_genes/becca_ir_list_coords.txt'
path_list = '/nas/longleaf/home/gabivla/lab/SV_mosquitos/IR_genes/cyp_gst_list_coords.txt'

df_becca = pd.read_csv(path_becca, sep="\t")
df_list = pd.read_csv(path_list, sep="\t", header=None)

df_becca = df_becca[['chrom', 'start', 'end', 'other name', '##gene ID']]

df_list[4] = df_list[3].str.extract(r'.+;Name=(.+);')[0]
df_list[5] = df_list[3].str.extract(r'ID=(.+);')[0]

df_list = df_list.drop(columns=[3])

df_list = df_list.rename(columns={0: 'chrom', 1: 'start', 2: 'end', 4: 'other name', 5: '##gene ID'})

df_merged = pd.concat([df_becca, df_list], ignore_index=True)

location_path = '/nas/longleaf/home/gabivla/lab/SV_mosquitos/location_information/mcclintock_output_locations.csv'
location_df = pd.read_csv(location_path, usecols=['Genome_code', 'country'])
location_dict = dict(zip(location_df['Genome_code'], location_df['country']))

base_columns = ['chromosome', 'position_start', 'position_end', 'strand', 'TE_name']
# update genome names in df_ref to match df_nonref
genome_cols = [col for col in df_ref.columns if col not in base_columns]

# update the genome columns names
for genome in genome_cols:
    if genome.startswith('SRR'):
        genome_updated = genome.replace('_1', '')
    elif genome.startswith('FEMALE') or genome.startswith('MALE'):
        genome_updated = genome.replace('_L001_L002_R1_001', '')
    # update the column name in the dataframe
    df_ref = df_ref.rename(columns={genome: genome_updated})

genome_cols = [col for col in df_ref.columns if col not in base_columns]

# Create rename dictionary and apply it once
rename_dict = {col: f"{col}_{location_dict.get(col, '')}" for col in genome_cols}
# drop column that is being a problem
df = df_ref.loc[:, ~df_ref.columns.str.contains("Unnamed")]
df = df.rename(columns=rename_dict)

# Extract unique countries from the renamed genome columns only
countries = set()
for original_col in genome_cols:
    new_col = rename_dict[original_col]
    country = new_col.split('_')[-1]
    countries.add(country)

df = pl.DataFrame(df)
TE_types = count_TE_occurence(df)
TE_types = TE_types['TE_name'].to_list()
df = df.to_pandas()

country_pairs_ref = list(combinations(countries, 2))
all_dfs_ref = {}
country_vcf_dict_ref = {}

for i, country in enumerate(countries):

    country_cols = [col for col in df.columns 
                   if col not in base_columns and col.endswith('_' + country)]
    
    if not country_cols:
        continue
    
    country_df = df[base_columns + country_cols]
    
    country_df = country_df.dropna(subset=country_cols, how='all')
    
    print(f"\n--- Data for {country} ---")
    frequency_df = calculate_frequency_TE(country_df, country)
    frequency_df = frequency_df[frequency_df[f'frequency_{country}'] > 0.5]
    
    df_ir_genes = ir_gene_regions(frequency_df, df_merged, country)
    df_ir_genes.to_csv(f'./TEs_in_IR_genes_{country}.csv', index=False)'''
    
from matplotlib.colors import ListedColormap

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re
import seaborn as sns
from matplotlib.colors import ListedColormap

list_file = os.listdir('./')

csv_files = [file for file in list_file if file.endswith('.csv')]

for file in csv_files:
    df = pd.read_csv(file)
    country_name = re.search('TEs_in_IR_genes_(.+).csv', file).group(1)
    mean_frequency = df[f'frequency_{country_name}'].mean()
    print(f'{country_name} mean frequency reference')
    print(mean_frequency)

    df['gene_family'] = df['gene_name'].str.extract(r'([A-Za-z]+)')
    df['TE_name'] = df['TE_name'].str.upper()
    df['TE_name'] = df['TE_name'].str.extract(
        r'(?:^|_)(LINE|SINE|LTR|TIR|MITE)(?:_|$)', expand=False
    ).fillna('UNKNOWN')
    te_counts = df.groupby(['gene_family', 'TE_name']).size().reset_index(name='count')

    df_pivot = te_counts.pivot(index='TE_name', columns='gene_family', values='count')

    df_pivot = df_pivot.fillna(0).astype(int)

    pastel_palette = sns.color_palette("pastel", len(df_pivot.columns))

    plt.figure(figsize=(15, 6))
    df_pivot.T.plot(kind='bar', width=0.8, color=pastel_palette)

    plt.xlabel("Gene Family", fontsize=15)
    plt.ylabel("Count", fontsize=15)
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.title(f"{country_name}", fontsize=20)
    plt.legend(title="TE Name", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(rotation=45, ha='right')
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    plt.tight_layout()
    plt.savefig(f'./TE_Distribution_{country_name}_50.svg')
    plt.close()
