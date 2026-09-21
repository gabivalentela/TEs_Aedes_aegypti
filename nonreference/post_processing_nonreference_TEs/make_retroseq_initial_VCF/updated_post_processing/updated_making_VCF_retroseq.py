import sys
import re
import pandas as pd
import polars as pl
from utils import Making_vcf

def updating_TE_names(df_variation_1000):
    # make dict from original names to new names
    path_dict = '../../../consensus_TEs_dict/consensus_file_dictionary.txt'
    with open(path_dict, 'r') as f:
        lines = f.readlines()
        name_dict = {}
        for line in lines:
            pattern = r'Original:(.+)\s\|\sUpdated:(.+)'
            match = re.search(pattern, line)
            if match:
                original_name = match.group(1)
                new_name = match.group(2)
                # remove whitespace from original_name and new_name
                original_name = original_name.strip()
                new_name = new_name.strip()
                name_dict[new_name] = original_name

    # make new column TE_name_original removing everything after _span* if that is there
    def extract_original_name(te_name):
        if isinstance(te_name, str) and '_span' in te_name:
            return te_name.split('_span')[0]
        return te_name
        
    df_variation_1000['TE_name_original'] = df_variation_1000['TE_name'].apply(extract_original_name)

    # update names in df_variation_1000
    df_variation_1000['TE_name_original'] = df_variation_1000['TE_name_original'].map(name_dict).fillna(df_variation_1000['TE_name'])

    # add the span information back to the original name if it was there
    def add_span_info(te_name, original_name):
        # Catch bad values
        if not isinstance(te_name, str):
            print(f"Non-string te_name encountered: {te_name} (type={type(te_name)})")
            return original_name

        if '_span' in te_name:
            span_info = te_name.split('_span')[1]
            return f"{original_name}_span{span_info}"
        else:
            return original_name

    df_variation_1000['TE_name'] = df_variation_1000.apply(lambda row: add_span_info(row['TE_name'], row['TE_name_original']), axis=1)

    # remove the TE_name_original column
    df_variation_1000.drop(columns=['TE_name_original'], inplace=True)

    return df_variation_1000

def main():
    """
    Main function to execute the workflow: finding BED files, processing TE data, and generating reports.
    """
    list_of_path_nonref = ['/work/users/g/a/gabivla/lab/SV_mosquitos/new_library_pipeline_nonreference/outputs', '/users/g/a/gabivla/new_library_pipeline_nonreference/outputs']    
    
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
    
    #dictionary to store different countries
    vcf_files = {}

    # Loop through each country and process the data
    for i, (country_nonref, df_country_nonref) in enumerate(nonref_country.items()):
        
        df_country_nonref = df_country_nonref.to_pandas()
        # remove the duplicates to make sure we don't have duplicated insertions
        te_content_all_genomes_nonref_SNPcheck = df_country_nonref.drop_duplicates()
        print(f'number of TEs before removing duplicates - {len(df_country_nonref)}')
        print(f'number of TEs after removing duplicates - {len(te_content_all_genomes_nonref_SNPcheck)}')

        te_content_all_genomes_nonref_SNPcheck = pl.from_pandas(te_content_all_genomes_nonref_SNPcheck)

        # making VCF file
        vcf_file = Making_vcf.make_vcf(te_content_all_genomes_nonref_SNPcheck, chromosome)

        # remove satellites and simple repeats from file
        vcf_file = Making_vcf.remove_non_tes(vcf_file)
        
        # add vcf_files to dict
        vcf_files[country_nonref] = vcf_file

    # merge all VCFs
    vcf_all_countries = Making_vcf.merge_vcf_files(vcf_files)
    
    # remove singletons in this final dataset
    #vcf_file = Making_vcf.remove_singletons(vcf_all_countries)
    
    # keeping order of TEs
    vcf_file = vcf_all_countries.sort(['chromosome', 'position_start', 'position_end'])

    # make vcf_file into pandas
    vcf_file = vcf_file.to_pandas()

    # update the names to the original library
    df_updated_names = updating_TE_names(vcf_file)
    
    # save VCF
    df_updated_names.to_csv(f'./nonreference_files/VCF_nonreference_{chromosome}_{variation}.csv')    
if __name__ == "__main__":
    main()