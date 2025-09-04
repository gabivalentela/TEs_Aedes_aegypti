import pandas as pd

reference_TEs = '/work/users/g/a/gabivla/lab/SV_mosquitos/VCF_genotypes_TEs/separate_chromosomes_windows_reference/merged_reference/VCF_reference_TEs_post_processing_600_absence_evidence.csv'
nonreference_TEs = '/work/users/g/a/gabivla/lab/SV_mosquitos/VCF_genotypes_TEs/separate_chromosomes_windows_reference/merged_nonreference/VCF_nonreference_post_processing.csv'
snps_file = '/proj/dschridelab/tvkent/aedes_data/VCFs/AaegL5_full_filtered_snps_110122_rep_map_100kbgenes_centromere_filtered_sites_mask.vcf.gz'

# counting number of TEs
ref_df = pd.read_csv(reference_TEs)
nonref_df = pd.read_csv(nonreference_TEs)
# remove column - SRR6768005 from both df
ref_df = ref_df.loc[:, ~ref_df.columns.str.contains('SRR6768005')]
nonref_df = nonref_df.loc[:, ~nonref_df.columns.str.contains('SRR6768005')]
print(ref_df.columns.to_list())
# get genotype columns
base_columns = ['chromosome', 'position_start', 'position_end', 'TE_name', 'strand']
genotype_columns_ref = ref_df.columns.difference(base_columns)
genotypes_df_ref = ref_df[genotype_columns_ref]
# remove rows that only contain 0/0 or only contain 1/1
genotypes_df_ref = genotypes_df_ref.loc[~(genotypes_df_ref == '0/0').all(axis=1)]
genotypes_df_ref = genotypes_df_ref.loc[~(genotypes_df_ref == '1/1').all(axis=1)]

# nonref
genotype_columns_nonref = nonref_df.columns.difference(base_columns)
genotypes_df_nonref = nonref_df[genotype_columns_nonref]
# remove rows that only contain 0/0 or only contain 1/1
genotypes_df_nonref = genotypes_df_nonref.loc[~(genotypes_df_nonref == '0/0').all(axis=1)]
genotypes_df_nonref = genotypes_df_nonref.loc[~(genotypes_df_nonref == '1/1').all(axis=1)]

print(f'the number of reference TEs is: {genotypes_df_ref.shape[0]}')
print(f'the number of non-reference TEs is: {genotypes_df_nonref.shape[0]}')

'''# counting number of SNPs
import pandas as pd
import gzip

def read_vcf_with_header(file_path):
    with gzip.open(file_path, 'rt') as f:  # Use gzip.open with 'rt' mode
        lines = f.readlines()
    
    # Find the header line (starts with single #, not ##)
    header_idx = None
    for i, line in enumerate(lines):
        if line.startswith('#') and not line.startswith('##'):
            header_idx = i
            break
    
    if header_idx is None:
        raise ValueError("No header line found")
    
    # Read from header line onwards
    return pd.read_csv(file_path, compression='gzip', sep='\t', 
                       skiprows=header_idx)

snps_df = read_vcf_with_header(snps_file)
# remove columns that start with JB_
snps_df = snps_df.loc[:, ~snps_df.columns.str.startswith('JB_')]
# get genotypes from columns
base_columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
# make a df with only genotype columns
genotype_columns = snps_df.columns.difference(base_columns)
genotypes_df = snps_df[genotype_columns]
# genotypes from columns - first information before the first ;
genotypes_df = genotypes_df.applymap(lambda x: x.split(':')[0] if isinstance(x, str) else x)
# remove any monomorphic SNPs - rows were all info are 0/0 or 1/1
# Remove rows that are monomorphic when ignoring missing data
def filter_polymorphic(df):
    return df[df.apply(lambda row: len(row[row != './.'].unique()) > 1, axis=1)]

polymorphic_df = filter_polymorphic(genotypes_df)
print(f'the number of polymorphic SNPs is: {polymorphic_df.shape[0]}')'''

'''# difference between two lists
unique_ref = set(list_ref) - set(list_snp)
unique_snp = set(list_snp) - set(list_ref)

print(f'Unique reference TEs: {unique_ref}')
print(f'Unique SNPs: {unique_snp}')'''