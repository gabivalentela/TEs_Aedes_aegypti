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
    chromosome = sys.argv[2]
    variation = int(variation)
    
    loc_info = Making_vcf.info_location('/nas/longleaf/home/gabivla/lab/SV_mosquitos/location_information/mcclintock_output_locations.csv')
    print(f'making bed list')
    
    list_of_bed = Making_vcf.find_bed_file_temp2(list_of_path_ref)
    list_of_bed = list_of_bed[:5]
     
    print(f'this is the bed list: {list_of_bed}')

    print(f'the amount of genomes is: {len(list_of_bed)}')
    
    print(f'making the dictionaries with the ref and nonref data')
    ref_dict_dd, nonref_dict_dd = Making_vcf.find_ref_nonref_df(list_of_bed, chromosome)
            
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

    vcf_files = {}

    # Loop through each country and process the data
    for i, (country_ref, df_country_ref) in enumerate(ref_country.items()):
        # keep only TEs that are this TE type - rnd_1_family_823_LINE_RTE_BovB
        df_country_ref = df_country_ref.filter(pl.col('TE_name').str.contains('rnd_1_family_823_LINE_RTE_BovB'))

        # Get the number of genomes for the current country
        number_of_genomes = df_country_ref['Genome_code'].n_unique()
        
        # Apply the necessary transformations and calculations
        te_content_all_genomes_ref_window = Making_vcf.making_windows_for_comparison(df_country_ref)
        print(te_content_all_genomes_ref_window)
        '''te_content_all_genomes_ref_SNPcheck = Making_vcf.checking_if_TEs_are_the_same_updated(te_content_all_genomes_ref_window, variation)
        print(te_content_all_genomes_ref_SNPcheck)
        # making VCF file
        vcf_file = Making_vcf.make_vcf(te_content_all_genomes_ref_SNPcheck, chromosome)
        
        # remove satellites and simple repeats from file
        vcf_file = Making_vcf.remove_non_tes(vcf_file)

        # remove TE families that are not present in the nonrefence dataset
        vcf_file = Making_vcf.remove_specific_te_families(vcf_file, list_nonref_ref_df)
        
        # add vcf_files to dict
        vcf_files[country_ref] = vcf_file

    # merge all VCFs
    vcf_all_countries = Making_vcf.merge_vcf_files(vcf_files)
    # split to windows
    vcf_all_countries = Making_vcf.split_vcf_into_windows(vcf_all_countries, 1000000)
    
    # make one VCF file for each window and make a mkdir for each chromosome
    #output_dir = f'/work/users/g/a/gabivla/lab/SV_mosquitos/VCF_genotypes_TEs/separate_chromosomes_windows/reference_initial_files/{chromosome}'
    output_dir = f'/work/users/g/a/gabivla/lab/SV_mosquitos/VCF_genotypes_TEs/separate_chromosomes_windows/testeeeeeee/{chromosome}'
    
    # Group by chromosome and window_id
    for (chromosome, window_id), group in vcf_all_countries.group_by(["chromosome", "window_id"]):
        # Use output_dir directly (it already has the chromosome name)
        os.makedirs(output_dir, exist_ok=True)

        # Drop the window_id column
        group = group.drop("window_id")
        
        # Define the output file name
        output_file = os.path.join(output_dir, f"VCF_TEs_reference_{chromosome}_{window_id}_initial_new.csv")
        
        # Save the group to CSV
        group.write_csv(output_file)'''
    
if __name__ == "__main__":
    main()