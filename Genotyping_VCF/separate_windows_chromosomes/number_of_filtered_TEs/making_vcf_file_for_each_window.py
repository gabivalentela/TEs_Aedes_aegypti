import sys
import os
import pandas as pd
import polars as pl
from utils import Making_vcf
import re

def main():
    """
    Main function to execute the workflow: finding BED files, processing TE data, and generating reports.
    """
    list_of_path_ref = ['/work/users/g/a/gabivla/lab/SV_mosquitos/running_mcclintock/general_srr', '/work/users/g/a/gabivla/lab/SV_mosquitos/running_mcclintock/US_genomes', '/work/users/g/a/gabivla/lab/SV_mosquitos/running_mcclintock/flycross', '/work/users/g/a/gabivla/lab/SV_mosquitos/running_mcclintock/extra_genomes']
    list_of_path_nonref = ['/work/users/g/a/gabivla/lab/SV_mosquitos/running_mcclintock_trimmed/general_srr', '/users/g/a/gabivla/running_mcclintock_trimmed_reads/general_srr', '/users/g/a/gabivla/running_mcclintock_trimmed_reads/US_genomes', '/users/g/a/gabivla/running_mcclintock_trimmed_reads/flycross', '/users/g/a/gabivla/running_mcclintock_trimmed_reads/extra_genomes']    
    #list_of_path_nonref = ['/users/g/a/gabivla/running_mcclintock_trimmed_reads/general_srr', '/users/g/a/gabivla/running_mcclintock_trimmed_reads/US_genomes', '/users/g/a/gabivla/running_mcclintock_trimmed_reads/flycross', '/users/g/a/gabivla/running_mcclintock_trimmed_reads/extra_genomes']    
    
    neutral_data_path = '/nas/longleaf/home/gabivla/lab/SV_mosquitos/neutral_data/making_SFS/my_script/homozygous_counts_by_country/SFS_files'
    
    # Value defined for variation in grouping of the transposable elements
    # How much can the positions vary to be considered the same
    variation = sys.argv[1]
    variation = int(variation)
    
    loc_info = Making_vcf.info_location('/nas/longleaf/home/gabivla/lab/SV_mosquitos/location_information/mcclintock_output_locations.csv')
    print(f'making bed list')
    
    list_of_bed = Making_vcf.find_bed_file_temp2(list_of_path_ref)
     
    print(f'this is the bed list: {list_of_bed}')

    print(f'the amount of genomes is: {len(list_of_bed)}')
    
    print(f'making the dictionaries with the ref and nonref data')
    ref_dict_dd, nonref_dict_dd = Making_vcf.find_ref_nonref_df(list_of_bed)
            
    # Get only the data of the most common TE among each genome in a dataframe (only lines that contain the top 1 - 5 TEs) separately
    # Do this for ref and non-ref
    ref_country = Making_vcf.make_df_per_country(ref_dict_dd, loc_info)
    # getting bam files
    list_bam_files = Making_vcf.find_bam_files(list_of_path_ref)
    # find cutoff values
    cutoff_values = Making_vcf.find_insert_size_cutoff('/work/users/g/a/gabivla/lab/SV_mosquitos/insert_sizes_bam/')
    # list of TE families that have reference and nonreference TEs
    path_table_te_families = '/nas/longleaf/home/gabivla/lab/SV_mosquitos/TE_families_reference_and_nonreference_TEs/te_families_ref_and_nonref_df_both_info.csv'
    list_nonref_ref_df = pl.read_csv(path_table_te_families)

    vcf_files_filtering = {}
    vcf_files_no_filtering = {}

    # Loop through each country and process the data
    for i, (country_ref, df_country_ref) in enumerate(ref_country.items()):

        # Get the number of genomes for the current country
        number_of_genomes = df_country_ref['Genome_code'].n_unique()
        
        # Apply the necessary transformations and calculations
        te_content_all_genomes_ref_window = Making_vcf.making_windows_for_comparison(df_country_ref)
        te_content_all_genomes_ref_SNPcheck = Making_vcf.checking_if_TEs_are_the_same_updated(te_content_all_genomes_ref_window, variation)
        
        # making VCF file
        vcf_file = Making_vcf.make_vcf(te_content_all_genomes_ref_SNPcheck)
        
        vcf_files_no_filtering[country_ref] = vcf_file
    
        # remove satellites and simple repeats from file
        vcf_file = Making_vcf.remove_non_tes(vcf_file)

        # remove TE families that are not present in the nonrefence dataset
        vcf_file = Making_vcf.remove_specific_te_families(vcf_file, list_nonref_ref_df)
        
        # add vcf_files to dict
        vcf_files_filtering[country_ref] = vcf_file

    # merge all VCFs
    vcf_all_countries_filtering = Making_vcf.merge_vcf_files(vcf_files_filtering)
    vcf_all_countries_no_filtering = Making_vcf.merge_vcf_files(vcf_files_no_filtering)

    outfile = open('number_of_TEs_filtered_reference.txt', 'a')
    outfile.write(f'number of TEs in the filtered VCF file: {len(vcf_all_countries_filtering)}\n')
    outfile.write(f'number of TEs in the VCF file without filtering: {len(vcf_all_countries_no_filtering)}\n')
    outfile.close()
    print(vcf_all_countries_filtering)
    print(vcf_all_countries_no_filtering)
    
if __name__ == "__main__":
    main()