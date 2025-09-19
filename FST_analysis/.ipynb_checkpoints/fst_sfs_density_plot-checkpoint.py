import sys
import pandas as pd
import polars as pl
import os
import re
import matplotlib.pyplot as plt
import numpy as np
import glob
from itertools import combinations

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

def make_one_big_dataframe(dict_TEs, condition):
    """
    Combines multiple Polars DataFrames into a single DataFrame with an additional column.

    Args:
        dict_TEs (dict): Dictionary where keys are file paths and values are Polars DataFrames.
        condition (str): Condition string to append to the 'Ind_code' column.

    Returns:
        pl.DataFrame: Combined Polars DataFrame with the additional 'Ind_code' column.
    """
    dfs = []
    
    for file, df in dict_TEs.items():
        # Make a copy of the DataFrame to avoid modifying the original
        df_copy = df.clone()
        
        # Extract the name using the regex pattern
        name = re.search(r'.+/.+/output_.+/(.+).fastq.+/results/.+', file).group(1)
        
        # Add the 'Ind_code' column with the specified name and condition
        df_copy = df_copy.with_columns([
            pl.lit(f'{name}{condition}').alias("Ind_code")
        ])
        
        # Append the modified DataFrame copy to the list
        dfs.append(df_copy)
    
    # Concatenate all DataFrames in the list
    general_df_complete = pl.concat(dfs, how="vertical")
    
    return general_df_complete

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
    
    # Get top 5 more present TEs
    df_first_5 = sorted_grouped.head(5)
    
    return df_first_5

def retrive_list_top_5_TEs(df):
    """
    Retrieves a list of the top 5 TEs from a DataFrame.

    Args:
        df (pl.DataFrame): DataFrame containing the top 5 TEs.

    Returns:
        list: List of top 5 TEs.
    """
    top_5_list = df['TE_name'].to_list()
    return top_5_list

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

def split_by_te_type(top_te, country_df):
    """
    Filters a DataFrame to retain only rows containing the specified top TE across all genomes.

    Args:
        top_te (str): The TE name to filter for. The function will look for this string in the 'TE_name' column.
        country_df (pl.DataFrame): Polars DataFrame containing TE data, with columns including 'TE_name'.

    Returns:
        pl.DataFrame: A DataFrame containing only the rows where the 'TE_name' matches the specified top TE.
    """
    # Drop the 'Ind_code' column if it exists
    country_df = country_df.drop('Ind_code')
    
    # Create an empty DataFrame with a defined schema for general use
    general_df = pl.DataFrame(schema={"chromosome": pl.Utf8, "position_start": pl.Utf8, "position_end": pl.Utf8, "TE_name": pl.Utf8, "score": pl.Utf8, "strand": pl.Utf8, "genome_ID": pl.Utf8})
    
    # Convert the 'TE_name' column to uppercase and split by '|' to standardize the TE names
    country_df = country_df.with_columns(pl.col('TE_name').str.split('|').list.get(0).str.to_uppercase())
    
    # Create a boolean mask to filter rows where 'TE_name' contains the specified top TE
    mask = country_df[:, 3].str.contains(top_te)
    
    # Use the boolean mask to create a new DataFrame containing only the matching rows
    selected_df = country_df.filter(mask).clone()
    
    # Rename columns to ensure consistency with expected schema
    selected_df.columns = ["chromosome", "position_start", "position_end", "TE_name", "score", "strand", "genome_ID", "Genome_code", "data_analysis", "country", "location"]
    
    return selected_df


def split_by_te_type_nonref(top_te, country_df):
    """
    Filters a DataFrame to retain only rows containing the specified top TE across all genomes.

    Args:
        top_te (str): The TE name to filter for. The function will look for this string in the 'TE_name' column.
        country_df (pl.DataFrame): Polars DataFrame containing TE data, with columns including 'TE_name'.

    Returns:
        pl.DataFrame: A DataFrame containing only the rows where the 'TE_name' matches the specified top TE.
    """
    # Drop the 'Ind_code' column if it exists
    country_df = country_df.drop('Ind_code')
    
    # Create an empty DataFrame with a defined schema for general use
    general_df = pl.DataFrame(schema={"chromosome": pl.Utf8, "position_start": pl.Utf8, "position_end": pl.Utf8, "TE_name": pl.Utf8, "score": pl.Utf8, "strand": pl.Utf8, "genome_ID": pl.Utf8})
    
    # Convert the 'TE_name' column to uppercase and split by '|' to standardize the TE names
    country_df = country_df.with_columns(pl.col('TE_name').str.split('|').list.get(0).str.to_uppercase())
    
    # Create a boolean mask to filter rows where 'TE_name' contains the specified top TE
    mask = country_df[:, 3].str.contains(top_te)
    
    # Use the boolean mask to create a new DataFrame containing only the matching rows
    selected_df = country_df.filter(mask).clone()
    
    # Rename columns to ensure consistency with expected schema
    selected_df.columns = ["chromosome", "position_start", "position_end", "TE_name", "score", "strand", "TE_information","genome_ID", "Genome_code", "data_analysis", "country", "location"]
    
    return selected_df

def making_windows_for_comparison(df_per_te):
    """
    Adds a 'Window' column to the DataFrame that indicates which genomic window each row falls into.

    The function divides each chromosome into fixed-size windows and assigns the start position of the window 
    to each row that overlaps with that window. The window size is set to 100,000 base pairs.

    Args:
        df_per_te (pl.DataFrame): Polars DataFrame containing genomic data with columns 'chromosome', 'position_start', and 'position_end'.

    Returns:
        pl.DataFrame: The original DataFrame with an additional 'Window' column indicating the start position of the genomic window each row falls into.
    """
    # Get unique chromosome values as a list
    chrm_list = df_per_te.select(pl.col("chromosome").unique()).to_series().to_list()
    
    # Initialize "Window" column
    df_per_te = df_per_te.with_columns(pl.lit(0).alias("Window"))
    
    window_size = 100000  # Window size in base pairs (100 Kb)
    for chrm in chrm_list:
        # Get the maximum position for the current chromosome
        max_position = df_per_te.filter(pl.col("chromosome") == chrm)["position_end"].max()
        
        for start_position in range(0, max_position + 1, window_size):
            end_position = start_position + window_size - 1
            
            # Update the 'Window' column based on overlap with the current window
            df_per_te = df_per_te.with_columns(
                pl.when((pl.col("chromosome") == chrm) &
                        (
                            (pl.col("position_start") >= start_position) &
                            (pl.col("position_start") <= end_position) |  # TE starts within the window
                            (pl.col("position_end") >= start_position) &
                            (pl.col("position_end") <= end_position) |  # TE ends within the window
                            (pl.col("position_start") <= start_position) &
                            (pl.col("position_end") >= end_position)  # TE spans the entire window
                        ))
                  .then(start_position)
                  .otherwise(pl.col("Window"))
                  .alias("Window")
            )
    
    return df_per_te

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

    # For each genome, check if any of the insertions are not present in the genome
    for genome_name, df in non_ref_retroseq.items():
        # Select the same unique combination columns from the current genome DataFrame
        genome_combinations = df.select(
            pl.col("chromosome"),
            pl.col("position_start"),
            pl.col("position_end"),
            pl.col("TE_name")
        ).unique()
        
        # Find the insertions in merged_df that are not in the genome
        missing_insertions = merged_df.join(genome_combinations, on=["chromosome", "position_start", "position_end", "TE_name"], how="anti")
        
        # If the insertion is not present, check the breakpoint 1 file
        match = re.search(r'(.+/results/retroseq/edited_results)/.+', genome_name)
        if match:
            general_path = match.group(1)
        
            # Use glob to find the breakpoint_1 file that ends with 'nonredundant.bed'
            breakpoint_1_files = glob.glob(f'{general_path}/breakpoint_1/*nonredundant.bed')
            
            if breakpoint_1_files:
                # Print or use the first file found (assuming there may be more than one match)
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
                
                # Now, check if any TE from missing_insertions are in non_result
                matching_insertions = missing_insertions.join(non_result, on=["chromosome", "position_start", "position_end", "TE_name"], how="inner")
                
                # Add matching_insertions to the original DataFrame
                updated_df = df.vstack(matching_insertions)
                # Update the dictionary with the new DataFrame
                updated_dfs[genome_name] = updated_df

            else:
                print(f'No "nonredundant.bed" file found in the breakpoint_1 directory. - {genome_name}')
        else:
            print(f"Unable to extract general path from genome name: {genome_name}")

    return updated_dfs

def checking_if_TEs_are_the_same(df, variation):
    """
    Filters a DataFrame to keep only one representative row for each unique Transposable Element (TE) 
    within the same window and chromosome, based on a specified variation in start and end positions.
    """

    # Initialize a list to store the final DataFrames
    final_dfs = []

    # Split the DataFrame into separate DataFrames based on Genome_code
    genome_code_dfs = {genome_code: group_df for genome_code, group_df in df.group_by(["Genome_code"])}
    
    for genome_code, group_df in genome_code_dfs.items():
            
        # Define columns to group by for comparison
        group_cols = ['chromosome', 'Window', 'TE_name']

        # Initialize an empty list to store the final rows
        final_rows = []
        num_removed_lines = 0  # Counter for removed lines

        # Group by these columns and iterate over groups
        for list_info, group in group_df.group_by(group_cols):
            num_rows = group.height

            # Initialize a set to track rows that should be removed
            rows_to_remove = set()

            # Iterate over the rows in the group
            for i in range(num_rows):
                info_i = group.row(i, named=True)
                position_start_i = info_i['position_start']
                position_end_i = info_i['position_end']
                te_name_i = info_i['TE_name']

                # Check against subsequent rows j in the group
                for j in range(i + 1, num_rows):
                    info_j = group.row(j, named=True)
                    position_start_j = info_j['position_start']
                    position_end_j = info_j['position_end']
                    te_name_j = info_j['TE_name']

                    # Check if the start or end positions differ by at most `variation` and TE names match
                    if (abs(position_start_j - position_start_i) <= variation or
                        abs(position_end_i - position_end_j) <= variation) and te_name_i == te_name_j:

                        # Mark the row info_j for removal
                        rows_to_remove.add(j)  # Mark row j for removal

            # After checking all rows, add only those not marked for removal
            for k in range(num_rows):
                if k not in rows_to_remove:
                    final_rows.append(group.row(k, named=True))

        # Create a new Polars DataFrame from final_rows and add it to the final_dfs list
        final_df = pl.DataFrame(final_rows)
        final_dfs.append(final_df)

    # Concatenate all final DataFrames into one
    merged_final_df = pl.concat(final_dfs)

    return merged_final_df

def calculate_te_frequency_removed_fixed(df, number_of_genomes):
    """
    Calculate the frequency of Transposable Elements (TEs) that are present in a specified number of genomes, 
    and filter out rows corresponding to these TEs from the original DataFrame.

    This function performs the following steps:
    1. Aggregates 'Genome_code' values for each TE by chromosome, position, and TE name.
    2. Counts the number of occurrences of each TE.
    3. Filters TEs that appear in the specified number of genomes.
    4. Filters out rows from the original DataFrame that correspond to these TEs.

    Parameters:
    df (pd.DataFrame): The input DataFrame containing columns 'chromosome', 'position_start', 
                       'position_end', 'TE_name', and 'Genome_code'.
    number_of_genomes (int): The number of genomes a TE must be present in to be considered.

    Returns:
    pd.DataFrame: A DataFrame with rows corresponding to TEs present in exactly the specified number of genomes
                  removed from the original DataFrame.
    """
    df = df.to_pandas()
    # Group by specified columns and aggregate the 'Genome_code' column
    grouped_df = df.groupby(['chromosome', 'position_start', 'position_end', 'TE_name'])['Genome_code'].agg(lambda x: '; '.join(set(x))).reset_index()
    # Count the frequency of occurrences
    frequency_df = df.groupby(['chromosome', 'position_start', 'position_end', 'TE_name']).size().reset_index(name='Frequency_of_TE')
    
    # Merge the two dataframes on the grouping columns
    result = grouped_df.merge(frequency_df, on=['chromosome', 'position_start', 'position_end', 'TE_name'])
    
    # Sort the result by 'Frequency_of_TE' in descending order
    result = result.sort_values(by='Frequency_of_TE', ascending=False)
    
    # Filter to keep only TEs with the specified frequency
    fixed_TEs = result[result['Frequency_of_TE'] == number_of_genomes]
    
    # Convert fixed_TEs to a list of tuples for filtering
    fixed_TE_tuples = set(fixed_TEs[['position_start', 'position_end', 'TE_name']].itertuples(index=False, name=None))
    
    # Filter the initial DataFrame to exclude rows with these values
    filtered_df = df[~df[['position_start', 'position_end', 'TE_name']].apply(tuple, axis=1).isin(fixed_TE_tuples)]

    # Count the frequency of occurrences
    new_frequency_df = filtered_df.groupby(['chromosome', 'position_start', 'position_end', 'TE_name']).size().reset_index(name='Frequency_of_TE')

    return new_frequency_df

def calculate_te_frequency_removed_fixed_nonref(df, number_of_genomes):
    """
    Calculate the frequency of Transposable Elements (TEs) that are present in a specified number of genomes, 
    and filter out rows corresponding to these TEs from the original DataFrame.

    This function performs the following steps:
    1. Aggregates 'Genome_code' values for each TE by chromosome, position, and TE name.
    2. Counts the number of occurrences of each TE.
    3. Filters TEs that appear in the specified number of genomes.
    4. Filters out rows from the original DataFrame that correspond to these TEs.

    Parameters:
    df (pd.DataFrame): The input DataFrame containing columns 'chromosome', 'position_start', 
                       'position_end', 'TE_name', and 'Genome_code'.
    number_of_genomes (int): The number of genomes a TE must be present in to be considered.

    Returns:
    pd.DataFrame: A DataFrame with rows corresponding to TEs present in exactly the specified number of genomes
                  removed from the original DataFrame.
    """
    df = df.to_pandas()
    
    # Step 1: Group by 'TE_information' and aggregate the 'Genome_code' values for each TE
    grouped_df = df.groupby(['TE_information'])['Genome_code'].agg(lambda x: '; '.join(set(x))).reset_index()
    
    # Step 2: Count the frequency of occurrences of each TE
    frequency_df = df.groupby(['TE_information']).size().reset_index(name='Frequency_of_TE')
    
    # Step 3: Merge the grouped DataFrame with the frequency DataFrame on 'TE_information'
    result = grouped_df.merge(frequency_df, on=['TE_information'])
    
    # Step 4: Sort the result by 'Frequency_of_TE' in descending order
    result = result.sort_values(by='Frequency_of_TE', ascending=False)
    
    # Step 5: Filter TEs that appear in exactly the specified number of genomes
    fixed_TEs = result[result['Frequency_of_TE'] == number_of_genomes]
    
    # Step 6: Convert fixed_TEs to a set of tuples for filtering
    fixed_TE_tuples = set(fixed_TEs[['TE_information']].itertuples(index=False, name=None))
    
    # Step 7: Filter the initial DataFrame to exclude rows with these values
    filtered_df = df[~df[['TE_information']].apply(tuple, axis=1).isin(fixed_TE_tuples)]
    
    # Step 8: Recount the frequency of occurrences in the filtered DataFrame
    new_frequency_df = filtered_df.groupby(['TE_information']).size().reset_index(name='Frequency_of_TE')

    return new_frequency_df


def sfs_data(frequency_of_tes, top_te):
    """
    Generate a DataFrame summarizing the occurrences of TEs based on their frequencies.

    This function calculates the frequency of each Transposable Element (TE) occurrence, 
    counts how many TEs have each frequency, and creates a DataFrame summarizing these counts.

    Parameters:
    frequency_of_tes (pd.DataFrame): DataFrame containing a 'Frequency_of_TE' column 
                                     with the frequencies of TEs.
    top_te (str): The name of the column to use as the TE name in the result DataFrame.

    Returns:
    pd.DataFrame: A DataFrame with the following columns:
                  - 'TE_name': The name of the TE (as provided by `top_te`).
                  - 'occurence': The unique frequencies of TEs.
                  - `top_te`: The count of TEs for each frequency.
    """
    frequency = frequency_of_tes['Frequency_of_TE']
    grouped = frequency.value_counts().sort_values(ascending=False)
    # Create a DataFrame with renamed columns
    result_df = pd.DataFrame({'TE_name': top_te, 'occurence': grouped.index, top_te: grouped.values})
    return result_df


def sfs_general_data(dataframes):
    """
    Create a DataFrame containing the frequency counts and proportions of all Transposable Elements (TEs) 
    analyzed across multiple input DataFrames.

    Parameters:
    - dataframes (list of pd.DataFrame): A list of DataFrames where each DataFrame contains frequency counts 
      for different TEs. Each DataFrame should have columns 'TE_name' and 'occurence', with the TE names 
      in the first row and occurrence counts in the second column.

    Returns:
    - pd.DataFrame: A DataFrame where each row represents a unique occurrence count, and each column represents 
      the proportion of that occurrence count for each TE.
    """
    
    # Find the maximum occurrence across all dataframes
    max_occurrence = max(df['occurence'].max() for df in dataframes)
    
    # Create a master DataFrame with all possible occurrences
    master_df = pd.DataFrame({'occurence': range(1, max_occurrence + 1)})
    
    # Create an empty dictionary to store the proportions for each TE
    te_proportions = {}

    # Process each DataFrame
    for df in dataframes:
        name = df.iloc[0]['TE_name']
        
        # Check for missing occurrences in the DataFrame
        missing_occurrences = set(range(1, max_occurrence + 1)) - set(df['occurence'])
        
        # Add missing occurrences with a count of 0
        for missing_occurrence in missing_occurrences:
            new_row = pd.DataFrame({'TE_name': [name], 'occurence': [missing_occurrence], name: [0]})
            df = pd.concat([df, new_row], ignore_index=True)

        df = df.sort_values(by=['occurence'])
        sel_df = df[['occurence', name]]
        total_count = sel_df[name].sum()
        proportions = sel_df[name] / total_count
        master_df = pd.merge(master_df, sel_df, on='occurence', how='left', suffixes=('', '_df'))
        te_proportions[name] = proportions

    # Add the proportions to the master_df
    for te, proportions in te_proportions.items():
        master_df[te] = proportions

    print(master_df)
    return master_df

def add_other_information(master_df, neutral_data_path, country):
    """
    Adds neutral information to the master DataFrame based on the specified country.

    Parameters:
    master_df (pd.DataFrame): The DataFrame containing the master SFS data.
    neutral_data_path (str): The path to the directory containing neutral data files.
    country (str): The country to match against the neutral data files.

    Returns:
    pd.DataFrame: The master DataFrame with added neutral information as proportions.
    """
    # Dictionary mapping country names to their corresponding names in the neutral data files
    country_map = {
        "USA": "USA",
        "Colombia": "Rio_Claro",
        "Kenya": "Kenya",
        "Senegal": "Senegal",
        "Brazil": "Brazil",
        "Gabon": "Gabon"
    }

    # Get the corresponding country_name from the dictionary
    country_name_to_match = country_map.get(country, country)
    
    list_of_neutral = os.listdir(neutral_data_path)
    print("Country to match:", country)
    print("Files in neutral data path:", list_of_neutral)
    
    for sfs in list_of_neutral:
        # Use a raw string for the regex pattern
        match = re.search(r'homozygous_counts_(\w+)_rep.+mask_folded.sfs', sfs)  
        if match:
            country_name = match.group(1)
            
            if country_name_to_match == country_name:
                print(f'Matching file: {sfs} -- to this country: {country}')
                neutral = os.path.join(neutral_data_path, sfs)
                # Reading neutral info from sfs file
                neutral_df = pd.read_csv(neutral, sep=r'\s+', header=None, engine='python')
                # removes first column
                neutral_df = neutral_df.iloc[:, 1:]
                neutral_new_df = pd.melt(neutral_df)
                neutral_new_df = neutral_new_df.rename(columns={"variable": "occurence", "value": "neutral"})
                total_neutral = neutral_new_df['neutral'].sum()
                proportions_neutral = neutral_new_df['neutral'] / total_neutral
                merged_df = pd.concat([master_df, proportions_neutral], axis=1)
                merged_df = merged_df.drop(columns=['occurence'])
                print(merged_df)
                return merged_df

def fold_my_sfs(master_df):
    """
    Fold the Site Frequency Spectrum (SFS) data for each column in the input DataFrame 
    except for the 'occurence' column. The folding process involves flipping the SFS 
    and adding it to itself, then retaining only the first half of the resulting array.
    NaN values are removed before the folding process, and all resulting folded SFS arrays 
    are padded with NaNs to ensure they have the same length.

    Parameters:
    - master_df (pd.DataFrame): A DataFrame where each column represents a different 
      TE's frequency counts, and the 'occurence' column represents the occurrence counts.

    Returns:
    - pd.DataFrame: A DataFrame where each column corresponds to the folded SFS of 
      the original columns from the input DataFrame. If the resulting DataFrame would be empty, 
      returns the original DataFrame with the 'occurence' column removed.
    """
    # Create a dictionary to store the folded SFS for each column
    folded_sfs_dict = {}
    max_length = 0
    
    # Iterate through columns, excluding 'occurence'
    for col in master_df.columns:
        if col != 'occurence':
            SFS = master_df[col].dropna().tolist()
            n = len(SFS)
            
            # Folded SFS array
            folded_sfs = np.zeros((n + 1) // 2)
            
            # Sum the extremes
            for i in range((n + 1) // 2):
                if i != n - i - 1:
                    folded_sfs[i] = SFS[i] + SFS[n - i - 1]
                else:
                    folded_sfs[i] = SFS[i]  # Middle element in case of odd length
            
            folded_sfs_dict[col] = folded_sfs
            max_length = max(max_length, len(folded_sfs))

    # Pad arrays with NaN to ensure all arrays have the same size
    for col in folded_sfs_dict:
        if len(folded_sfs_dict[col]) < max_length:
            folded_sfs_dict[col] = np.pad(folded_sfs_dict[col], (0, max_length - len(folded_sfs_dict[col])), constant_values=np.nan)

    # Create a DataFrame from the folded SFS dictionary
    folded_sfs_df = pd.DataFrame(folded_sfs_dict)
    
    # Check if folded_sfs_df is empty, and return the original DataFrame if true
    if folded_sfs_df.empty:
        return master_df.drop(columns=['occurence'])

    return folded_sfs_df

def check_for_tes_between_individuals(dict_individuals, variation=100):
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

            
def plot_sfs_master_df(ax, df, suf, country):
    """
    Plots the SFS master DataFrame as a grouped bar plot.

    Parameters:
    ax (matplotlib.axes.Axes): The matplotlib Axes object to plot on.
    df (pd.DataFrame): The DataFrame containing the SFS data.
    suf (str): The suffix for the plot title.
    country (str): The country for the plot title.

    Returns:
    None
    """
    import numpy as np
    import matplotlib.pyplot as plt
    
    # Identify the variable column names dynamically
    variable_columns = [col for col in df.columns]

    # Set the width of the bars
    bar_width = 0.1

    # Create an array of indices for the x-axis positions
    x = np.arange(len(df.index)) + 1  # Adding 1 to index to start from 1

    # Create the grouped bar plot dynamically for variable columns
    colormap = plt.get_cmap('Blues', len(variable_columns))  # Use 'Blues' colormap
    
    for i, col in enumerate(variable_columns):
        ax.bar(x + (i - len(variable_columns) / 2) * bar_width, df[col], bar_width, label=col, color=colormap(i), edgecolor='black')

    # Set x-axis labels and title
    ax.set_xlabel('Number of genomes with TE insertion', fontsize=16)  # Adjust the font size
    ax.set_ylabel('Proportion', fontsize=16)  # Adjust the font size
    ax.set_xticks(x)  # Set x-axis labels to the modified 'x' array
    ax.set_title(f"{country} - {suf}", fontsize=18)

    # Adjust the legend font size and location
    legend = ax.legend(loc='upper right')
    for item in legend.get_texts():
        item.set_fontsize(14)  # You can adjust the legend font size here

def make_vcf_file(df, output_file):
    # Convert to pandas if not already
    df = df.to_pandas() if hasattr(df, 'to_pandas') else df

    # Create info column more efficiently
    df['info'] = df['chromosome'].astype(str) + '_' + \
                 df['position_start'].astype(str) + '_' + \
                 df['position_end'].astype(str) + '_' + \
                 df['TE_name'].astype(str)
    
    df = df[['info', 'Genome_code']]
    
    # Create a pivot table - much faster than iterating
    vcf_df = df.pivot_table(
        index='info', 
        columns='Genome_code', 
        aggfunc=lambda x: '1/1' if len(x) > 0 else '0/0', 
        fill_value='0/0'
    )
    
    # Reset index to make 'info' a column
    vcf_df = vcf_df.reset_index()
    
    # Rename and split info column
    vcf_df = vcf_df.rename(columns={'info': 'CHROM_POS_END'})
    
    # Save to CSV
    #vcf_df.to_csv(output_file, index=False)
    
    return vcf_df

def calculate_al_freq(df):
    # get the number of samples
    number_of_samples = df.shape[1] - 1
    # calculate frequency of first allele - 1/1
    df['1/1'] = df.iloc[:, 1:].apply(lambda x: (x == '1/1').sum() / number_of_samples, axis=1)
    # calculate frequency of second allele - 0/1
    df['0/0'] = df.iloc[:, 1:].apply(lambda x: (x == '0/0').sum() / number_of_samples, axis=1)

    return df

def calculate_al_freq(df):
    # Retain CHROM_POS_END column and get the number of samples
    number_of_samples = df.shape[1] - 1
    # Calculate frequency of first allele - 1/1
    df['1/1'] = df.iloc[:, 1:].apply(lambda x: (x == '1/1').sum() / number_of_samples, axis=1)
    # Calculate frequency of second allele - 0/0
    df['0/0'] = df.iloc[:, 1:].apply(lambda x: (x == '0/0').sum() / number_of_samples, axis=1)
    
    return df

def calculate_FST(freq_p1, freq_p2):
    # Align DataFrames by CHROM_POS_END
    merged = pd.merge(freq_p1, freq_p2, on="CHROM_POS_END", suffixes=('_p1', '_p2'))
    
    # Calculate allele frequencies
    p1 = merged['1/1_p1']
    p2 = merged['1/1_p2']

    # Calculate q1 and q2
    q1 = 1 - p1
    q2 = 1 - p2

    # Calculate total allele frequency
    p_t = (p1 + p2) / 2
    q_t = 1 - p_t

    # Calculate expected heterozygosity
    hs_1 = 2 * p1 * q1
    hs_2 = 2 * p2 * q2
    hs = (hs_1 + hs_2) / 2

    ht = 2 * p_t * q_t

    # Calculate FST
    fst = (ht - hs) / ht

    # Return CHROM_POS_END and FST as a DataFrame
    fst_df = pd.DataFrame({
        'CHROM_POS_END': merged['CHROM_POS_END'],
        'FST': fst
    })

    return fst_df

def main():
    """
    Main function to execute the workflow: finding BED files, processing TE data, and generating reports.
    """
    list_of_path_ref = ['/work/users/g/a/gabivla/lab/SV_mosquitos/running_mcclintock/general_srr', '/work/users/g/a/gabivla/lab/SV_mosquitos/running_mcclintock/US_genomes', '/work/users/g/a/gabivla/lab/SV_mosquitos/running_mcclintock/flycross', '/work/users/g/a/gabivla/lab/SV_mosquitos/running_mcclintock/extra_genomes']
    list_of_path_nonref = ['/work/users/g/a/gabivla/lab/SV_mosquitos/running_mcclintock_trimmed/general_srr', '/users/g/a/gabivla/running_mcclintock_trimmed_reads/general_srr', '/users/g/a/gabivla/running_mcclintock_trimmed_reads/US_genomes', '/users/g/a/gabivla/running_mcclintock_trimmed_reads/flycross', '/users/g/a/gabivla/running_mcclintock_trimmed_reads/extra_genomes']
    
    neutral_data_path = '/nas/longleaf/home/gabivla/lab/SV_mosquitos/neutral_data/making_SFS/my_script/homozygous_counts_by_country/SFS_files'
    
    # Value defined for variation in grouping of the transposable elements
    # How much can the positions vary to be considered the same
    variation = sys.argv[1]
    variation = int(variation)

    loc_info = info_location('/nas/longleaf/home/gabivla/lab/SV_mosquitos/location_information/mcclintock_output_locations.csv')
    print(f'making bed list')
    
    list_of_bed = find_bed_file_temp2(list_of_path_ref)
    non_ref_list_of_bed = find_nonref_retroseq(list_of_path_nonref)
    print(f'this is the bed list: {list_of_bed}')
    
    print(f'the amount of genomes is: {len(list_of_bed)}')
    
    print(f'making the dictionaries with the ref and nonref data')
    ref_dict_dd, non_ref_dict_dd_temp2 = find_ref_nonref_df(list_of_bed)
    non_ref_dict_dd = make_nonref_retroseq_df(non_ref_list_of_bed)
    #checking breakpoints 
    non_ref_retroseq = check_breakpoints_retroseq(non_ref_dict_dd)
    #grouping individuals
    non_ref_dict_dd = check_for_tes_between_individuals(non_ref_retroseq, variation)
        
    print(f'making one big df')
    general_df_nonref = make_one_big_dataframe(non_ref_dict_dd, '.nonr')
    
    # Find the top 5 TEs present in the genomes on the nonref results
    top_te_df = count_TE_occurence(general_df_nonref)
    print(f'NonReference - {top_te_df}')

    # Retrieve a list with the top 5 TEs found
    top_5_list = retrive_list_top_5_TEs(top_te_df)

    # Get only the data of the most common TE among each genome in a dataframe (only lines that contain the top 1 - 5 TEs) separately
    # Do this for ref and non-ref
    ref_country = make_df_per_country(ref_dict_dd, loc_info)
    nonref_country = make_df_per_country(non_ref_dict_dd, loc_info)
    
    # everything for ref
    # Initialize the dictionary to store the VCF filenames for country pairs
    country_vcf_dict_ref = {}

    # Generate all combinations of two countries
    country_pairs_ref = list(combinations(ref_country.keys(), 2))

    # Iterate over each country and its DataFrame
    vcf_files_ref = {}
    
    outfile_ref = open('FST_values_ref_percentage.txt', 'a')

    for i, (country_ref, df_country_ref) in enumerate(ref_country.items()):

        # Create VCF file for the current country
        vcf_file = make_vcf_file(df_country_ref, f'output_VCF_file_{country_ref}.vcf')
        vcf_files_ref[country_ref] = vcf_file  # Store it in a dictionary for reuse

    # Now, create the country pair dictionary
    for country1, country2 in country_pairs_ref:
        country_vcf_dict_ref[(country1, country2)] = (vcf_files_ref[country1], vcf_files_ref[country2])

    fig, axes = plt.subplots(4, 4, figsize=(40, 30))    
    
    # Iterate over country pairs and their corresponding VCF data
    for i, (pair, vcf_tuple) in enumerate(country_vcf_dict_ref.items()):
        country1, country2 = pair
        vcf_country_1, vcf_country_2 = vcf_tuple

        # Calculate frequency of alleles
        counts_1 = calculate_al_freq(vcf_country_1)
        counts_2 = calculate_al_freq(vcf_country_2)

        # Calculate FST values
        fst_pop1_pop2 = calculate_FST(counts_1, counts_2)
        fst_pop1_pop2.to_csv(f'/work/users/g/a/gabivla/lab/SV_mosquitos/FST_data/TEs/FST_values_{country1}_{country2}_ref.csv', index=False)

        # Get the FST values as a list or numpy array
        fst_values = fst_pop1_pop2['FST'].values

        # Compute the histogram counts and bin edges
        counts, bins = np.histogram(fst_values, bins=20, range=(0, 1))

        # Convert counts to percentages
        total_count = sum(counts)
        percentages = (counts / total_count) * 100

        # Plot the histogram as percentages
        axes[i // 4, i % 4].bar(bins[:-1], percentages, width=np.diff(bins), edgecolor='black', alpha=0.7)

        # Label the plot
        axes[i // 4, i % 4].set_title(f'FST Histogram: {country1} vs {country2}')
        axes[i // 4, i % 4].set_xlabel('FST value')
        axes[i // 4, i % 4].set_ylabel('Percentage')
        axes[i // 4, i % 4].set_xlim(0, 1)
        axes[i // 4, i % 4].set_ylim(0, 100)

        # Optionally, write the mean FST value for each pair
        fst_value = np.mean(fst_pop1_pop2['FST'])
        outfile_ref.write(f'FST value between {country1} and {country2}: {fst_value}\n')

    # Adjust layout
    plt.tight_layout()

    # Save the figure
    plt.savefig('fst_histograms_ref.png')

    #---------------------------------------------------------------------------------------------------
    #everything for nonref
    # Initialize the dictionary to store the VCF filenames for country pairs
    country_vcf_dict_nonref = {}

    # Generate all combinations of two countries
    country_pairs_nonref = list(combinations(nonref_country.keys(), 2))

    # Iterate over each country and its DataFrame
    vcf_files_nonref = {}

    outfile_nonref = open('FST_values_nonref_percentage.txt', 'a')

    for country_nonref, df_country_nonref in nonref_country.items():
        # Create VCF file for the current country
        vcf_file = make_vcf_file(df_country_nonref, f'output_VCF_file_{country_nonref}.vcf')
        vcf_files_nonref[country_nonref] = vcf_file  # Store it in a dictionary for reuse

    # Now, create the country pair dictionary
    for country1, country2 in country_pairs_nonref:
        country_vcf_dict_nonref[(country1, country2)] = (vcf_files_nonref[country1], vcf_files_nonref[country2])

    # Create a figure with subplots
    fig, axes = plt.subplots(4, 4, figsize=(40, 30))

    # Iterate over country pairs and their corresponding VCF data
    for i, (pair, vcf_tuple) in enumerate(country_vcf_dict_nonref.items()):
        country1, country2 = pair
        vcf_country_1, vcf_country_2 = vcf_tuple

        # Calculate frequency of alleles
        counts_1 = calculate_al_freq(vcf_country_1)
        counts_2 = calculate_al_freq(vcf_country_2)

        # Calculate FST values
        fst_pop1_pop2 = calculate_FST(counts_1, counts_2)
        fst_pop1_pop2.to_csv(f'/work/users/g/a/gabivla/lab/SV_mosquitos/FST_data/TEs/FST_values_{country1}_{country2}_nonref.csv', index=False)

        # Get the FST values as a list or numpy array
        fst_values = fst_pop1_pop2['FST'].values

        # Compute the histogram counts and bin edges
        counts, bins = np.histogram(fst_values, bins=20, range=(0, 1))

        # Convert counts to percentages
        total_count = sum(counts)
        percentages = (counts / total_count) * 100

        # Plot the histogram as percentages
        axes[i // 4, i % 4].bar(bins[:-1], percentages, width=np.diff(bins), edgecolor='black', alpha=0.7)

        # Label the plot
        axes[i // 4, i % 4].set_title(f'FST Histogram: {country1} vs {country2}')
        axes[i // 4, i % 4].set_xlabel('FST value')
        axes[i // 4, i % 4].set_ylabel('Percentage')
        axes[i // 4, i % 4].set_xlim(0, 1)
        axes[i // 4, i % 4].set_ylim(0, 100)

        # Optionally, write the mean FST value for each pair
        fst_value = np.mean(fst_pop1_pop2['FST'])
        outfile_nonref.write(f'FST value between {country1} and {country2}: {fst_value}\n')

    # Adjust layout
    plt.tight_layout()

    # Save the figure
    plt.savefig('fst_histograms_nonref.png')
        
if __name__ == "__main__":
    main()
