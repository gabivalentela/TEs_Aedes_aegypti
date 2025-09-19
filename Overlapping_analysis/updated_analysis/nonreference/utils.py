import sys
import pandas as pd
import polars as pl
import os
import re
import matplotlib.pyplot as plt
import numpy as np
import glob

@staticmethod
def divide_by_countries(df, location_dict):
    base_cols = ['chromosome', 'position_start', 'position_end', 'TE_name', 'strand']

    genome_cols = [col for col in df.columns if col not in base_cols]
    
    rename_dict = {col: f"{col}_{location_dict.get(col, '')}" for col in genome_cols}
    df = df.rename(columns=rename_dict)

    countries = set()
    for original_col in genome_cols:
        new_col = rename_dict[original_col]
        country = new_col.split('_')[-1]
        countries.add(country)
    
    dfs = {}

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

        # convert to polars
        country_df = pl.DataFrame(country_df)

        # Store the DataFrame in the dictionary
        dfs[country] = country_df

    return dfs

def read_data_annotation(path_annotation):
    data = pd.read_csv(
        path_annotation,
        sep='\t',
        comment='#',  # Skip lines starting with '#'
        header=None,  # GFF files don't usually have column headers
    )
    # Assign column names based on GFF3 specification
    data.columns = [
        'seqid', 'source', 'type', 'start', 'end', 'score',
        'strand', 'phase', 'attributes'
    ]
    return data

def info_location(loc_info):
    """
    Reads the location information from a CSV file into a Polars DataFrame.

    Args:
        loc_info (str): Path to the CSV file containing location information.

    Returns:
        pl.DataFrame: Polars DataFrame with location information.
    """
    df = pl.read_csv(loc_info)
    return df

def find_bed_file_temp2(list_of_initial_path):
    """
    Searches for BED files with a specific naming pattern in a list of initial paths.

    Args:
        list_of_initial_path (list of str): List of directory paths to search.

    Returns:
        list of str: List of paths to BED files ending with 'nonredundant.bed'.
    """
    list_of_bed = []
    for initial_path in list_of_initial_path:
        directory = os.listdir(initial_path)
        for code in directory:
            if 'output_' in code:
                path = f'{initial_path}/{code}'
                file_list = os.listdir(path)
                results = os.listdir(f'{path}/{file_list[0]}')
                if 'results' in results:
                    result_files = os.listdir(f'{path}/{file_list[0]}/results/')
                    if 'temp2' in result_files:
                        files = os.listdir(f'{path}/{file_list[0]}/results/temp2')
                        for file in files:
                            if file.endswith('nonredundant.bed'):
                                bed_file = f'{path}/{file_list[0]}/results/temp2/{file}'
                                list_of_bed.append(bed_file)
    return list_of_bed

def find_nonref_retroseq(list_of_initial_path):
    """
    Finds and returns a list of paths to non-redundant BED files produced by RetroSeq.
    If a specific file or directory is not found, it prints a message.

    Args:
        list_of_initial_path (list of str): A list of paths to the initial directories to search for RetroSeq output.

    Returns:
        list of str: A list of paths to the 'nonredundant.bed' files found within the specified directories.
    """
    list_of_bed = []
    for initial_path in list_of_initial_path:
        try:
            directory = os.listdir(initial_path)
        except FileNotFoundError:
            print(f"Initial path not found: {initial_path}")
            continue
        
        for code in directory:
            if 'output_' in code:
                path = os.path.join(initial_path, code)
                try:
                    file_list = os.listdir(path)
                except FileNotFoundError:
                    print(f"Directory not found: {path}")
                    continue
                
                for file in file_list:
                    if '_75' in file:
                        sub_path = os.path.join(path, file)
                        try:
                            results = os.listdir(sub_path)
                        except FileNotFoundError:
                            print(f"Results directory not found: {sub_path}")
                            continue
                        
                        if 'results' in results:
                            results_path = os.path.join(sub_path, 'results')
                            try:
                                result_files = os.listdir(results_path)
                            except FileNotFoundError:
                                print(f"Results subdirectory not found: {results_path}")
                                continue
                            
                            if 'retroseq' in result_files:
                                retroseq_path = os.path.join(results_path, 'retroseq', 'edited_results', 'breakpoint_6')
                                try:
                                    files = os.listdir(retroseq_path)
                                except FileNotFoundError:
                                    print(f"RetroSeq results path not found: {retroseq_path}")
                                    continue
                                
                                for file_1 in files:
                                    if file_1.endswith('nonredundant.bed'):
                                        bed_file = os.path.join(retroseq_path, file_1)
                                        list_of_bed.append(bed_file)
    
    return list_of_bed

def make_nonref_retroseq_df(list_of_bed):
    """
    Creates a dictionary of Polars DataFrames from a list of BED file paths.

    This function reads each BED file from the provided list of file paths, converts it into a Polars DataFrame, 
    and processes the DataFrame by renaming columns and extracting Transposable Element (TE) names. 
    The processed DataFrames are stored in a dictionary with the file paths as keys.

    Args:
        list_of_bed (list of str): A list of file paths to BED files containing TE insertion data.

    Returns:
        dict: A dictionary where each key is a BED file path and each value is a Polars DataFrame 
              containing the processed TE insertion data.
    """
    non_ref_dict_polars = {}
    
    for file in list_of_bed:
        # Read each BED file into a Polars DataFrame, skipping the first line
        general_df = pl.read_csv(file, separator='\t', skip_rows=1, has_header=False)
        general_df = general_df.rename({"column_1": "chromosome", 
                                        "column_2": "position_start", 
                                        "column_3": "position_end", 
                                        "column_4": "TE_name", 
                                        "column_5": "score", 
                                        "column_6": "strand"}) 
        
        non_result = general_df.with_columns([
            pl.col("TE_name").str.split('|').list.get(0).alias("TE_name")
        ])
        
        non_ref_dict_polars[file] = non_result
    
    return non_ref_dict_polars

def find_ref_nonref_df(list_of_files):
    """
    Reads BED files and filters them into reference and non-reference categories.

    Args:
        list_of_files (list of str): List of paths to BED files.

    Returns:
        tuple: Two dictionaries containing Polars DataFrames for reference and non-reference data.
            - ref_dict_polars (dict): Dictionary with reference data DataFrames.
            - non_ref_dict_polars (dict): Dictionary with non-reference data DataFrames.
    """
    non_ref_dict_polars = {}
    ref_dict_polars = {}
    
    for file in list_of_files:
        # Read each BED file into a Polars DataFrame, skipping the first line
        general_df = pl.read_csv(file, separator='\t', skip_rows=1, has_header=False)
        general_df = general_df.rename({"column_1": "chromosome", "column_2": "position_start", "column_3": "position_end", "column_4": "TE_name", "column_5": "score", "column_6": "strand"}) 
        
        ref_result = general_df.filter(pl.col("TE_name").str.contains("|reference|", literal=True))
        non_result = general_df.filter(pl.col("TE_name").str.contains("|non-reference|", literal=True))

        ref_result = ref_result.with_columns([pl.col("TE_name").str.split('|').list.get(0).alias("TE_name")])
        non_result = non_result.with_columns([pl.col("TE_name").str.split('|').list.get(0).alias("TE_name")])

        # Store Polars DataFrames in dictionaries
        non_ref_dict_polars[file] = non_result
        ref_dict_polars[file] = ref_result

    return ref_dict_polars, non_ref_dict_polars

def check_breakpoints_retroseq(non_ref_retroseq):
    # Make a general df containing all TE insertions in all genomes
    all_unique_combinations = []

    for df in non_ref_retroseq.values():
        unique_combinations = df.select(
            pl.col("chromosome"),
            pl.col("position_start"),
            pl.col("position_end"),
            pl.col("TE_name"), 
            pl.col("score"),
            pl.col("strand")
        ).unique()
        all_unique_combinations.append(unique_combinations)
    
    # Concatenate all unique combinations into one DataFrame
    merged_df = pl.concat(all_unique_combinations).unique()

    updated_dfs = {}

    for genome_name, df in non_ref_retroseq.items():
        genome_combinations = df.select(
            pl.col("chromosome"),
            pl.col("position_start"),
            pl.col("position_end"),
            pl.col("TE_name")
        ).unique()
        
        # Find the insertions in merged_df that are not in the genome
        missing_insertions = merged_df.join(genome_combinations, on=["chromosome", "position_start", "position_end", "TE_name"], how="anti")
        
        match = re.search(r'(.+/results/retroseq/edited_results)/.+', genome_name)
        if match:
            general_path = match.group(1)
        
            # Use glob to find the breakpoint_1 file that ends with 'nonredundant.bed'
            breakpoint_1_files = glob.glob(f'{general_path}/breakpoint_1/*nonredundant.bed')
            
            if breakpoint_1_files:
                breakpoint_1 = breakpoint_1_files[0]
                print(f"Found breakpoint_1 file: {breakpoint_1}")
                general_df = pl.read_csv(breakpoint_1, separator='\t', skip_rows=1, has_header=False)
                
                # Rename the columns appropriately
                general_df = general_df.rename({"column_1": "chromosome", "column_2": "position_start", 
                                                "column_3": "position_end", "column_4": "TE_name", 
                                                "column_5": "score", "column_6": "strand"}) 
                general_df = general_df.drop(["score", "strand"])
                
                # Split TE_name if it contains "|" and take the first part
                non_result = general_df.with_columns([pl.col("TE_name").str.split('|').list.get(0).alias("TE_name")])
                
                matching_insertions = missing_insertions.join(non_result, on=["chromosome", "position_start", "position_end", "TE_name"], how="inner")
                
                updated_df = df.vstack(matching_insertions)
                updated_dfs[genome_name] = updated_df

            else:
                print(f'No "nonredundant.bed" file found in the breakpoint_1 directory. - {genome_name}')
        else:
            print(f"Unable to extract general path from genome name: {genome_name}")

    return updated_dfs

def check_for_tes_between_individuals(dict_individuals, variation):
    # Initialize an empty list to store DataFrames
    df_list = []
    
    # Iterate through the dictionary and process each DataFrame
    for genome_name, df in dict_individuals.items():
        # Extract genome code using regex
        #genome_code = re.search('.+/output_.+/(.+)/results.+', genome_name).group(1)
        genome_code = genome_name
        # Add genome_code as a new column
        df = df.with_columns([pl.lit(genome_code).alias('genome_code')])
        
        # Append the modified DataFrame to the list
        df_list.append(df)
    
    # Concatenate all DataFrames into one big DataFrame
    big_df = pl.concat(df_list, how="vertical")
    
    # Sort by TE_name and position_start
    big_df = big_df.sort(['TE_name', 'position_start'])
    
    # Initialize an empty list to store TE_information for each row
    te_information_list = []
    
    # Iterate over the DataFrame rows
    num_rows = big_df.height
    i = 0  # Initialize the index for iteration
    
    while i < num_rows:
        info_i = big_df.row(i, named=True)
        chromosome_i = info_i['chromosome']
        position_start_i = info_i['position_start']
        position_end_i = info_i['position_end']
        te_name_i = info_i['TE_name']
        
        # Set the window start and end for the current group
        min_start_position = position_start_i
        max_end_position = position_end_i
        chosen_te_information = f"{chromosome_i}_{min_start_position}_{max_end_position}_{te_name_i}"
        
        # Assign initial row the TE information
        te_information_list.append(chosen_te_information)
        
        # Check subsequent rows within the same TE_name for variation
        j = i + 1
        while j < num_rows:
            info_j = big_df.row(j, named=True)
            chromosome_j = info_j['chromosome']
            position_start_j = info_j['position_start']
            position_end_j = info_j['position_end']
            te_name_j = info_j['TE_name']
            
            # Check if TE names match and positions are within the allowed variation
            if te_name_i == te_name_j and chromosome_i == chromosome_j and (
                abs(position_start_j - position_start_i) <= variation or
                abs(position_end_j - position_end_i) <= variation):
                
                # Update the TE_information to be consistent for all matching rows
                te_information_list.append(chosen_te_information)
                
                # Increment the inner index
                j += 1
            else:
                # No further matches for the current group
                break
        
        # Move the outer index to the next row after the last matched row
        i = j
    
    # Add the new TE_information column to the original DataFrame
    big_df = big_df.with_columns([pl.Series(te_information_list).alias("TE_information")])
    
    # Now split the DataFrame back into individual genome DataFrames
    result_dict = {}
    
    # Group by genome_code (passing 'genome_code' as a list)
    for group_tuple in big_df.group_by(['genome_code']):
        # Unpack the group identifier (genome_code) and the DataFrame (df_group)
        genome_code, df_group = group_tuple
        genome_code = genome_code[0]
        
        # Keep only the relevant columns ('TE_information' and original columns)
        #removing te_information since I won't use it anymore
        df_group = df_group.select(['chromosome', 'position_start', 'position_end', 'TE_name', 'score', 'strand', 'TE_information'])
        #df_group = df_group.select(['chromosome', 'position_start', 'position_end', 'TE_name', 'score', 'strand'])
        
        # Store each group in the dictionary with genome_code as the key
        result_dict[genome_code] = df_group
    
    return result_dict

def make_df_per_country(dict_TEs, location):
    """
    Creates separate DataFrames for each country by merging TE data with location information.

    Args:
        dict_TEs (dict): Dictionary where keys are file paths and values are Polars DataFrames.
        location (pl.DataFrame): Polars DataFrame containing location information.

    Returns:
        dict: Dictionary where keys are country names and values are merged DataFrames for each country.
    """
    df = location.clone()  # Create a copy of the location DataFrame to avoid modifying the original

    country_dfs = {}  # Dictionary to hold DataFrames for each country

    for file, bed_file in dict_TEs.items():
        name = re.search(r'.+/.+/output_.+/(.+)/results/.+', file).group(1)
        bed_df = bed_file
        bed_df = bed_df.with_columns([pl.lit(f'{name}').alias("genome_ID")])
        if '_L' in name:
            bed_df = bed_df.with_columns([pl.col("genome_ID").str.split('_L').list.get(0).alias("Genome_code")])      
        elif 'JB_' in name:
            bed_df = bed_df.with_columns([pl.col("genome_ID").replace('R1', 'R2').alias("Genome_code")]) 
        elif '_1' in name:
            bed_df = bed_df.with_columns([pl.col("genome_ID").str.split('_').list.get(0).alias("Genome_code")])  
        else:
            bed_df = bed_df.with_columns(pl.lit(None).alias("Genome_code"))

        merged_df = bed_df.join(df, on='Genome_code', how='left', coalesce=True)
        country_name = merged_df['country'][0]

        if country_name not in country_dfs:
            country_dfs[country_name] = merged_df
        else:
            country_dfs[country_name] = pl.concat([country_dfs[country_name], merged_df])

    return country_dfs