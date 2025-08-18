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
    list_of_path_ref = [
        '/work/users/g/a/gabivla/lab/SV_mosquitos/running_mcclintock/general_srr',
        '/work/users/g/a/gabivla/lab/SV_mosquitos/running_mcclintock/US_genomes',
        '/work/users/g/a/gabivla/lab/SV_mosquitos/running_mcclintock/flycross',
        '/work/users/g/a/gabivla/lab/SV_mosquitos/running_mcclintock/extra_genomes'
    ]

    # getting bam files
    list_bam_files = Making_vcf.find_bam_files(list_of_path_ref)

    # find cutoff values
    cutoff_values = Making_vcf.find_insert_size_cutoff('/work/users/g/a/gabivla/lab/SV_mosquitos/insert_sizes_bam/')

    # list of TE families that have reference and nonreference TEs
    path_table_te_families = '/nas/longleaf/home/gabivla/lab/SV_mosquitos/TE_families_reference_and_nonreference_TEs/te_families_ref_and_nonref_df_both_info.csv'
    list_nonref_ref_df = pl.read_csv(path_table_te_families)

    vcf_all_countries = sys.argv[1]
    chr_window = sys.argv[2]

    # read vcf_all_country as pl df
    vcf_all_countries = pl.read_csv(vcf_all_countries, n_rows=20)

    # get genome codes were the TE is absent
    genome_codes_no_TE = Making_vcf.get_genome_codes_no_TE(vcf_all_countries)
    print('checking deletion reads')

    # add bam files
    updated_genotypes = Making_vcf.check_bam_files(genome_codes_no_TE, list_bam_files, cutoff_values)
    print('BAM files checked')

    # update VCF
    updated_vcf = Making_vcf.update_vcf(vcf_all_countries, updated_genotypes)
    print('VCF updated - updated absent TEs')

    # remove TEs that are fixed in all individuals or missing in only one individual
    # updated_vcf = Making_vcf.remove_fixed_almostfixed_TEs(updated_vcf)
    print('VCF updated - removed fixed TEs')

    # save VCF
    updated_vcf.to_csv(f'/work/users/g/a/gabivla/lab/SV_mosquitos/VCF_genotypes_TEs/separate_chromosomes_windows_reference/reference_final/VCF_TEs_reference_{chr_window}.csv')

if __name__ == "__main__":
    main()