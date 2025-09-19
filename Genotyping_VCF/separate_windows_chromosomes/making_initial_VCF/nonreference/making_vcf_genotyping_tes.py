import sys
import pandas as pd
import polars as pl
from utils import Making_vcf

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
    
    non_ref_list_of_bed = Making_vcf.find_nonref_retroseq(list_of_path_nonref)
    
    print(f'this is the bed list: {non_ref_list_of_bed}')

    print(f'the amount of genomes is: {len(non_ref_list_of_bed)}')
    
    print(f'making the dictionaries with the ref and nonref data')
    non_ref_dict_dd = Making_vcf.make_nonref_retroseq_df(non_ref_list_of_bed, chromosome)
    #checking breakpoints 
    non_ref_retroseq = Making_vcf.check_breakpoints_retroseq(non_ref_dict_dd)
    #grouping individuals
    non_ref_dict_dd = Making_vcf.check_for_tes_between_individuals(non_ref_retroseq, variation)
            
    # Get only the data of the most common TE among each genome in a dataframe (only lines that contain the top 1 - 5 TEs) separately
    # Do this for ref and non-ref
    nonref_country = Making_vcf.make_df_per_country(non_ref_dict_dd, loc_info)

    # list of TE families that have reference and nonreference TEs
    path_table_te_families = '/nas/longleaf/home/gabivla/lab/SV_mosquitos/TE_families_reference_and_nonreference_TEs/te_families_ref_and_nonref_df_both_info.csv'
    list_nonref_ref_df = pl.read_csv(path_table_te_families)

    #dictiionary to store different countries
    vcf_files = {}

    # Loop through each country and process the data
    for i, (country_nonref, df_country_nonref) in enumerate(nonref_country.items()):

        # Get the number of genomes for the current country
        number_of_genomes = df_country_nonref['Genome_code'].n_unique()
            
        # Apply the necessary transformations and calculations
        te_content_all_genomes_nonref_window = Making_vcf.making_windows_for_comparison(df_country_nonref)
        te_content_all_genomes_nonref_SNPcheck = Making_vcf.checking_if_TEs_are_the_same_updated(te_content_all_genomes_nonref_window, variation)
    
        # making VCF file
        vcf_file = Making_vcf.make_vcf(te_content_all_genomes_nonref_SNPcheck, chromosome)

        # remove satellites and simple repeats from file
        vcf_file = Making_vcf.remove_non_tes(vcf_file)
        
        # remove TE families that are not present in the nonrefence dataset
        vcf_file = Making_vcf.remove_specific_te_families(vcf_file, list_nonref_ref_df)

        # add vcf_files to dict
        vcf_files[country_nonref] = vcf_file

    # merge all VCFs
    vcf_all_countries = Making_vcf.merge_vcf_files(vcf_files)
    
    # remove singletons in this final dataset
    vcf_file = Making_vcf.remove_singletons(vcf_all_countries)
    
    # keeping order of TEs
    vcf_file = vcf_file.sort(['chromosome', 'position_start', 'position_end'])
    # save VCF
    vcf_file.write_csv(f'/work/users/g/a/gabivla/lab/SV_mosquitos/VCF_genotypes_TEs/separate_chromosomes_windows_reference/nonreference_files/VCF_nonreference_{chromosome}.csv')
    
if __name__ == "__main__":
    main()