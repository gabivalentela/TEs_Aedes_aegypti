import sys
from utils import FilesInterpretation
import pandas as pd
from utils import LD_analysis

# read VCF
vcf_ref = '/work/users/g/a/gabivla/lab/SV_mosquitos/VCF_genotypes_TEs/separate_chromosomes_windows_reference/merged_reference/VCF_reference_TEs_post_processing_600_absence_evidence.csv'
vcf_nonref = '/work/users/g/a/gabivla/lab/SV_mosquitos/VCF_genotypes_TEs/separate_chromosomes_windows_reference/merged_nonreference/VCF_nonreference_post_processing.csv'
info_loc = '/nas/longleaf/home/gabivla/lab/SV_mosquitos/location_information/mcclintock_output_locations.csv'

df_ref = pd.read_csv(vcf_ref)
df_nonref = pd.read_csv(vcf_nonref)

# remove first column
#df_ref = df_ref.iloc[:, 1:]
df_nonref = df_nonref.iloc[:, 1:]

# Read only necessary columns from location file
location_df = pd.read_csv(info_loc, usecols=['Genome_code', 'country'])
location_dict = dict(zip(location_df['Genome_code'], location_df['country']))

ref_df_colombia = FilesInterpretation.divide_by_countries(df_ref, location_dict)
nonref_df_colombia = FilesInterpretation.divide_by_countries(df_nonref, location_dict)

genes_with_regulatory_region = {
    'Chr': 'AaegL5_2',
    'position_start': 41628484 - 10000,  # initial position scheme NCBI
    'position_end': 41861946 + 10000,    # final postion + 10KB
    'gene_name': 'rdl'
}

#get SNPs info from VCF
variant_info = [{'variant_name': 'AaegL5_2_41847790_A301S', 'Chromosome': 'AaegL5_2', 'Position_start': 41847790}]
        
vcf_path = '/nas/longleaf/home/gabivla/lab/SV_mosquitos/breakpoint_analysis/ld_analysis/rdl_gene/A301S_RDL/filtered_variants.recode.vcf'
df_SNPs_variants = LD_analysis.find_variant_information(vcf_path, variant_info)

#make table in the same format as the TE one and removing genotypes that are not homozygous for the SNP one before calculating the LD
df_SNPs_variants_country = LD_analysis.add_country_name_SNPs(df_SNPs_variants, info_loc)
#complete_SNP_only_homozygous = LD_analysis.remove_heterozygous_for_SNPs(df_SNPs_variants_country)

# Reference
filtered_df_reference = LD_analysis.get_gene_region(ref_df_colombia, genes_with_regulatory_region)

# Merge chromosome, positions, and strand     
filtered_df_reference['TE_identifier'] = (
    filtered_df_reference['chromosome'] + '_' + 
    filtered_df_reference['position_start'].astype(str) + '_' + 
    filtered_df_reference['position_end'].astype(str) + '_' + 
    filtered_df_reference['TE_name'] + 
    filtered_df_reference['strand']
)

complete_df_ref = LD_analysis.merge_with_te_table(df_SNPs_variants_country, filtered_df_reference)
complete_df_ref = LD_analysis.remove_columns_with_none_pairs(complete_df_ref)

complete_df_ref = LD_analysis.get_encoding_SNPs(complete_df_ref)

#Calculating LD
D_r_values_ref, genotype_data_ref = LD_analysis.calculate_LD_with_estimator(complete_df_ref, significance_threshold=0.2)
ld_table_ref = LD_analysis.make_interpret_table(D_r_values_ref)

# Save the table to a CSV file
ld_table_ref.to_csv(f'Reference_ld_table_flycross_SNPs.csv')

# Retrieve significant values - r squared
significant_rsquared_ref, list_significant_tes_ref = LD_analysis.retrieve_significant_rsquared(
    ld_table_ref, significant_threshold=0.2, output_filename='significant_r_squared_ref.csv'
)

LD_analysis.plot_ld_heatmap_r2(ld_table_ref, list_significant_tes_ref, output_file='ld_heatmap_r2_ref.svg')

genotype_info_ref = LD_analysis.sign_rsquared_retrieve_genotype(genotype_data_ref, significant_rsquared_ref)
# Save to file
genotype_info_ref.to_csv(f'genotype_info_Reference.csv')

# Retrieve only the ones that have the same genotype as the variant
same_genotype_as_variant_ref, list_genotype_te_ref = LD_analysis.same_genotype(genotype_info_ref)
LD_analysis.plot_ld_heatmap_r2(ld_table_ref, list_genotype_te_ref, output_file='same_genotypes_ld_heatmap_r2_ref.svg')

# Save to file
same_genotype_as_variant_ref.to_csv(f'same_genotype_as_variant_Reference.csv')

# Nonreference
print('Making filtered df')
filtered_df_nonreference = LD_analysis.get_gene_region(nonref_df_colombia, genes_with_regulatory_region)

# Merge chromosome, positions, and strand     
filtered_df_nonreference['TE_identifier'] = (
    filtered_df_nonreference['chromosome'] + '_' + 
    filtered_df_nonreference['position_start'].astype(str) + '_' + 
    filtered_df_nonreference['position_end'].astype(str) + '_' + 
    filtered_df_nonreference['TE_name'] + 
    filtered_df_nonreference['strand']
)

complete_df_nonref = LD_analysis.merge_with_te_table(df_SNPs_variants_country, filtered_df_nonreference)
complete_df_nonref = LD_analysis.remove_columns_with_none_pairs(complete_df_nonref)

complete_df_nonref = LD_analysis.get_encoding_SNPs(complete_df_nonref)

# Calculating LD
D_r_values_nonref, genotype_data_nonref = LD_analysis.calculate_LD_with_estimator(complete_df_nonref, significance_threshold=0.2)
ld_table_nonref = LD_analysis.make_interpret_table(D_r_values_nonref)

# Save the table to a CSV file
ld_table_nonref.to_csv(f'Non-Reference_ld_table_flycross_SNPs.csv')

# Retrieve significant values - r squared
significant_rsquared_nonref, list_significant_tes_nonref = LD_analysis.retrieve_significant_rsquared(
    ld_table_nonref, significant_threshold=0.2, output_filename='significant_r_squared_nonref.csv'
)
LD_analysis.plot_ld_heatmap_r2(ld_table_nonref, list_significant_tes_nonref, output_file='ld_heatmap_r2_nonref.svg')

genotype_info_nonref = LD_analysis.sign_rsquared_retrieve_genotype(genotype_data_nonref, significant_rsquared_nonref)

# Save to file
genotype_info_nonref.to_csv(f'genotype_info_Non-Reference.csv')

# Retrieve only the ones that have the same genotype as the variant
same_genotype_as_variant_nonref, list_genotype_te_nonref = LD_analysis.same_genotype(genotype_info_nonref)
LD_analysis.plot_ld_heatmap_r2(ld_table_nonref, list_genotype_te_nonref, output_file='same_genotypes_ld_heatmap_r2_nonref.svg')

# Save to file
same_genotype_as_variant_nonref.to_csv(f'same_genotype_as_variant_Non-Reference.csv')