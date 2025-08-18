import os
import pandas as pd
from itertools import combinations
import matplotlib.pyplot as plt 
import polars as pl

def calculate_frequency_TE(df, country_name):
    """
    Calculate the frequency of '1/1' genotypes for each TE insertion.
    
    Parameters:
    df - DataFrame with TE data
    country_name - Name of the country (for reference, not used in calculation)
    
    Returns:
    DataFrame with frequency column added
    """
    # Get columns that contain genotype data (excluding base columns)
    base_cols = ['chromosome', 'position_start', 'position_end', 'TE_name', 'strand']
    genotype_cols = [col for col in df.columns if col not in base_cols]
    
    # Count '1/1' occurrences in each row
    df['count_1/1'] = df[genotype_cols].apply(lambda row: (row == '1/1').sum(), axis=1)
    
    # Calculate frequency (proportion of samples with 1/1)
    total_samples = len(genotype_cols)
    df['frequency'] = df['count_1/1'] / total_samples
    
    return df

def read_data_annotation(path_annotation):
    data = pd.read_csv(
        path_annotation,
        sep='\t',
        comment='#',
        header=None,
    )
    # Assign column names based on GFF3 specification
    data.columns = [
        'seqid', 'source', 'type', 'start', 'end', 'score',
        'strand', 'phase', 'attributes'
    ]
    return data

import pandas as pd
from mygene import MyGeneInfo
import matplotlib.pyplot as plt
from collections import Counter
import seaborn as sns

def extract_go_terms(gene_list):
    import pandas as pd
    from mygene import MyGeneInfo
    
    mg = MyGeneInfo()
    go_dict = {}
    
    for gene in gene_list:
        try:
            gene_info = mg.getgene(gene, fields='go', species=7159)
            
            if gene_info and 'go' in gene_info:
                go_terms = []
                
                for go_type in ['BP', 'CC', 'MF']:
                    if go_type in gene_info['go']:
                        terms = gene_info['go'][go_type]
                        if isinstance(terms, list):
                            go_terms.extend([term['term'] for term in terms])
                        elif isinstance(terms, dict):
                            go_terms.append(terms['term'])
                
                go_dict[gene] = "; ".join(go_terms) if go_terms else None
            else:
                print(f"no GO terms {gene}")
                go_dict[gene] = None
        except Exception as e:
            print(f"Error getting GO terms {gene}: {e}")
            go_dict[gene] = None
    
    return go_dict

def analyze_go_terms(gene_ids, country):

    none_elements = 0 
    all_go_terms = []
    individual_terms = []
    
    for gene_id in gene_ids:
        gene_id = gene_id.replace('gene_id=', '')
        go_terms_dict = extract_go_terms([gene_id])
        
        if gene_id not in go_terms_dict or go_terms_dict[gene_id] is None:
            none_elements += 1
            print(f"no GO terms {gene_id}")
        else:
            go_terms = go_terms_dict[gene_id]
            all_go_terms.append(go_terms)
            
            if go_terms:
                terms = go_terms.split("; ")
                individual_terms.extend(terms)
    
    print(f"no GO terms: {none_elements}")
    print(f"no GO terms: {len(all_go_terms)}")
    
    term_counts = Counter(individual_terms)
    
    unique_go_terms = set(individual_terms)
    print(f"unique Go terms: {len(unique_go_terms)}")
    
    print("\nmost common GO terms:")
    for term, count in term_counts.most_common(20):
        print(f"{term}: {count}")
    
    plot_go_term_frequencies(term_counts, country)
    
    return {
        'none_count': none_elements,
        'gene_with_go_count': len(all_go_terms),
        'unique_terms': unique_go_terms,
        'term_counts': term_counts
    }

def plot_go_term_frequencies(term_counts, country,top_n=20):
    most_common = term_counts.most_common(top_n)
    terms = [item[0] for item in most_common]
    counts = [item[1] for item in most_common]
    
    df = pd.DataFrame({'GO Term': terms, 'Count': counts})
    
    plt.figure(figsize=(12, 8))
    
    sns.barplot(x='Count', y='GO Term', data=df, palette='viridis')
    
    plt.title('Frequency GO terms - 20 most commons', fontsize=16)
    plt.xlabel('Number of Occurrences', fontsize=12)
    plt.ylabel('GO terms', fontsize=12)
    
    plt.tight_layout()
    
    plt.savefig(f'go_term_frequencies_{country}.svg', dpi=300, bbox_inches='tight')
    
    return plt

overlap_nonref_files = '/work/users/g/a/gabivla/lab/SV_mosquitos/overlap_analysis/'
files = os.listdir(overlap_nonref_files)
# keep only files from list that have _ref on them
files = [file for file in files if '_nonref' in file]
# make df from the first file in the list
df_overlap = pd.read_csv(overlap_nonref_files + files[0])
'''attribute_columns = df_overlap['attributes']
# split this column with ';'
attribute_columns = attribute_columns.str.split(';', expand=True)
# count how many time each combination of geneid and exon id appears - column with index 2 and 0
count = attribute_columns.groupby([2, 0]).size()'''

# get VCF to check frequency in each country
df = pd.read_csv('/work/users/g/a/gabivla/lab/SV_mosquitos/VCF_genotypes_TEs/separate_chromosomes/merged_nonreference/VCF_nonreference_post_processing.csv')

# fixing df
df = df.drop(columns=['Unnamed: 0'])

# Read only necessary columns from location file
location_path = '/nas/longleaf/home/gabivla/lab/SV_mosquitos/location_information/mcclintock_output_locations.csv'
location_df = pd.read_csv(location_path, usecols=['Genome_code', 'country'])
location_dict = dict(zip(location_df['Genome_code'], location_df['country']))

# Define base columns that should not be modified
base_cols = ['chromosome', 'position_start', 'position_end', 'TE_name', 'strand']

# Get genome columns at once
genome_cols = [col for col in df.columns if col not in base_cols]

# Create rename dictionary and apply it once
rename_dict = {col: f"{col}_{location_dict.get(col, '')}" for col in genome_cols}
df = df.rename(columns=rename_dict)

# Extract unique countries from the renamed genome columns only
countries = set()
for original_col in genome_cols:
    new_col = rename_dict[original_col]
    country = new_col.split('_')[-1]
    countries.add(country)

country_vcf_dict_ref = {}
country_pairs_nonref = list(combinations(countries, 2))

all_dfs_ref = {}
percentage_overlapping = {}

# Split data by country
for i, country in enumerate(countries):

    # Get all columns for this country from the renamed genome columns
    country_cols = [col for col in df.columns 
                   if col not in base_cols and col.endswith('_' + country)]
    
    if not country_cols:
        continue

    # Create a subset with base columns and columns for this country
    country_df = df[base_cols + country_cols]
    
    # Drop rows where all country columns are NA
    country_df = country_df.dropna(subset=country_cols, how='all')
    
    print(f"\n--- Data for {country} ---")
    # calculate frequencies
    frequency_df = calculate_frequency_TE(country_df, country)
    # remove rows were frequency == 0
    frequency_df = frequency_df[frequency_df['frequency'] > 0]
    merged_frequency_df = frequency_df.merge(
    df_overlap[['TE_name', 'chromosome', 'position_start', 'position_end', 'attributes']],
    on=['TE_name', 'chromosome', 'position_start', 'position_end'],
    how='left' 
    )

    #if 'attributes' not NaN:
    matched_rows = merged_frequency_df[merged_frequency_df['attributes'].notna()]
    
    # get rows that have NaN in 'attributes'
    non_matched_rows = merged_frequency_df[merged_frequency_df['attributes'].isna()]
    
    # number of TEs overlapping genes
    TE_overlapping = len(matched_rows)
    total_number_TEs = len(frequency_df)

    percentage_overlapping[country] = {'number_te_overlapping': TE_overlapping, 'total_number_TEs': total_number_TEs, 'percentage': TE_overlapping/total_number_TEs}
    
    # Save the merged DataFrame for this country
    all_dfs_ref[country] = matched_rows, non_matched_rows

# Plot the results
# Extract values for plotting
countries = list(percentage_overlapping.keys())
percentages = [percentage_overlapping[country]['percentage'] * 100 for country in countries]  # Convert to percentage

# Create the bar plot
plt.figure(figsize=(10, 6))
plt.bar(countries, percentages, color='c', alpha=0.7, edgecolor='black')

# Add labels above bars
for i, perc in enumerate(percentages):
    plt.text(i, perc + 1, f"{perc:.1f}%", ha='center', fontsize=10, fontweight='bold')

# Labels and title
plt.xlabel("Country")
plt.ylabel("Percentage of TEs Overlapping Genes (%)")
plt.title("Percentage of Transposable Elements (TEs) Overlapping Exons by Country")
plt.xticks(rotation=45)
plt.ylim(0, max(percentages) + 5)

# Show the plot
plt.savefig('percentage_TEs_overlapping_genes.svg')

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

count_genes = {}  # Store list of top genes for each country
gene_frequencies = {}  # Store frequency of top genes
TE_frequencies_no_overlap = {}
TE_frequencies_overlap = {}
gene_ids = {}  # Store gene IDs for each country

for country, dfs in all_dfs_ref.items():
    df = dfs[0]
    not_overlapping = dfs[1]

    attribute_columns = df['attributes'].str.split(';', expand=True)

    # Get all genes and their counts
    all_gene_counts = attribute_columns[2].value_counts()

    # get gene id for the ones that have insertions
    unique_genes = attribute_columns[2].unique()
    # store gene ids
    gene_ids[country] = (unique_genes, all_gene_counts)

# which genes don't have overlapping TEs
path_annotation = '../../../../data_annotation/VectorBase-54_AaegyptiLVP_AGWG.gff'
df_annotation = read_data_annotation(path_annotation)
# retrieve the rows that have 'exon' in the 'type' column
exons = df_annotation[df_annotation['type'] == 'exon']
# get gene ids from the attributes column
exon_ids = exons['attributes'].str.split(';', expand=True)[2].unique()
# convert to list
exon_ids = exon_ids.tolist()

# number of genes in the annotation file
print(f'Number of genes in the annotation file: {len(exon_ids)}')

for country, gene_ids_info in gene_ids.items():
    gene_ids = gene_ids_info[0]
    gene_counts = gene_ids_info[1]

    # count how many genes have each number of insertions
    gene_counts = gene_counts.value_counts().sort_index()
    # plot this distribution
    plt.figure(figsize=(10, 6))
    gene_counts.plot(kind='bar', color='orange', edgecolor='black')
    plt.xlabel("Number of Insertions")
    plt.ylabel("Number of Genes")
    plt.title(f"Distribution of Genes by Number of Insertions in {country}")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f'distribution_genes_by_number_insertions_{country}.svg')

    results = analyze_go_terms(gene_ids, country)

    outfile = f'gene_ids_{country}.txt'
    with open(outfile, 'w') as f:
        for gene_id in gene_ids:
            f.write(f"{gene_id}\n")
            # number of genes that have insertions
            f.write(f'Number of genes with insertions in {country}: {len(gene_ids)}\n')
            # number of genes that don't have insertions
            f.write(f'Number of genes without insertions in {country}: {len(exon_ids) - len(gene_ids)}\n')
            # percentage of genes without insertions
            f.write(f'Percentage of genes without insertions in {country}: {(len(exon_ids) - len(gene_ids)) / len(exon_ids) * 100:.2f}%\n')
            # percentage of genes with insertions
            f.write(f'Percentage of genes with insertions in {country}: {len(gene_ids) / len(exon_ids) * 100:.2f}%\n')
    

