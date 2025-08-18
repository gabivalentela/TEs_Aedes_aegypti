import sys
import pandas as pd
import polars as pl
import os
import re
import matplotlib.pyplot as plt
import numpy as np
import glob
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import seaborn as sns
import pysam
import re
import multiprocessing as mp
from functools import partial
from collections import defaultdict

# make class with functions
class Making_vcf:
    @staticmethod
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
    
    @staticmethod
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
    
    @staticmethod
    def find_ref_nonref_df(list_of_files, chromosome):
        """
        Reads BED files and filters them into reference and non-reference categories, 
        while also filtering rows based on the specified chromosome.

        Args:
            list_of_files (list of str): List of paths to BED files.
            chromosome (str): The chromosome to filter the rows by.

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

            # Filter by the chromosome
            filtered_df = general_df.filter(pl.col("chromosome") == chromosome)

            # Filter reference and non-reference rows
            ref_result = filtered_df.filter(pl.col("TE_name").str.contains("|reference|", literal=True))
            non_result = filtered_df.filter(pl.col("TE_name").str.contains("|non-reference|", literal=True))

            # Clean up TE_name column
            ref_result = ref_result.with_columns([pl.col("TE_name").str.split('|').list.get(0).alias("TE_name")])
            non_result = non_result.with_columns([pl.col("TE_name").str.split('|').list.get(0).alias("TE_name")])

            # Store Polars DataFrames in dictionaries
            non_ref_dict_polars[file] = non_result
            ref_dict_polars[file] = ref_result

        return ref_dict_polars, non_ref_dict_polars
    
    @staticmethod
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

            # remove unnecessary columns
            merged_df = merged_df.drop(["data-analysis-ngs-mapper-2", "score", "genome_ID"])

            if country_name not in country_dfs:
                country_dfs[country_name] = merged_df
            else:
                country_dfs[country_name] = pl.concat([country_dfs[country_name], merged_df])

        return country_dfs
    
    @staticmethod
    def find_bam_files(general_path):
        list_bam_files = {}
        for path in general_path:
            list_files = os.listdir(path)
            for file in list_files:
                new_path = os.path.join(path, file)
                files = os.listdir(new_path)
                # Check if the first two characters are uppercase
                for genome in files:
                    if len(genome) >= 2 and genome[:2].isupper():
                        if '_1' in genome and '_L001_L002_R1_001' not in genome:
                            genome_name = genome.replace('_1', '')
                        if '_L001_L002_R1_001' in genome:
                            genome_name = genome.replace('_L001_L002_R1_001', '')
                        new_new_path = os.path.join(new_path, genome, 'intermediate', 'mapped_reads')
                        new_list = os.listdir(new_new_path)
                        for bam in new_list:
                            if bam.endswith('.bam'):
                                list_bam_files[genome_name] = os.path.join(new_new_path, bam)
                                
        return list_bam_files

    @staticmethod
    def making_windows_for_comparison(df_per_te: pl.DataFrame) -> pl.DataFrame:
        """
        Adds a 'Window' column to indicate which genomic window(s) each TE overlaps.

        For TEs that overlap multiple windows, a separate row is created for each window.

        Args:
            df_per_te (pl.DataFrame): DataFrame with 'chromosome', 'position_start', 'position_end'.

        Returns:
            pl.DataFrame: Expanded DataFrame with 'Window' column indicating the window start position.
        """
        window_size = 100000  # size of each genomic window
        
        chrm_list = df_per_te.select(pl.col("chromosome").unique()).to_series().to_list()

        windowed_rows = []  # will collect dictionaries for all overlapping rows

        for chrm in chrm_list:
            df_chr = df_per_te.filter(pl.col("chromosome") == chrm)
            max_position = df_chr["position_end"].max()

            for start_position in range(0, max_position + 1, window_size):
                end_position = start_position + window_size - 1

                # Filter TEs overlapping this window
                overlapping = df_chr.filter(
                    (pl.col("position_end") >= start_position) & (pl.col("position_start") <= end_position)
                )

                for row in overlapping.iter_rows(named=True):
                    row_copy = row.copy()
                    row_copy["Window"] = start_position
                    windowed_rows.append(row_copy)

        # Convert the list of dicts back to a Polars DataFrame
        return pl.DataFrame(windowed_rows)
    
    @staticmethod
    def checking_if_TEs_are_the_same_updated(df, variation):
        """
        Merges similar Transposable Elements (TEs) within the same window and chromosome based on 
        a specified variation in start/end positions. Creates new names reflecting the merged spans.
        """
        final_dfs = []
        genome_code_dfs = {genome_code: group_df for genome_code, group_df in df.group_by(["Genome_code"])}
        
        for genome_code, group_df in genome_code_dfs.items():
            group_cols = ['chromosome', 'Window', 'TE_name']
            final_rows = []
            
            # Group by these columns and iterate over groups
            for list_info, group in group_df.group_by(group_cols):
                num_rows = group.height
                
                # Track which rows are part of which merged groups
                merged_groups = []
                processed_indices = set()
                
                # Iterate over the rows in the group
                for i in range(num_rows):
                    if i in processed_indices:
                        continue
                        
                    current_group = []
                    info_i = group.row(i, named=True)
                    current_group.append(info_i)
                    processed_indices.add(i)
                    
                    # Check against subsequent rows
                    for j in range(i + 1, num_rows):
                        if j in processed_indices:
                            continue
                            
                        info_j = group.row(j, named=True)
                        
                        # Check if positions are within variation
                        if (abs(info_j['position_start'] - info_i['position_start']) <= variation or
                            abs(info_i['position_end'] - info_j['position_end']) <= variation) and \
                            info_i['TE_name'] == info_j['TE_name']:
                            current_group.append(info_j)
                            processed_indices.add(j)
                    
                    if current_group:
                        merged_groups.append(current_group)
                
                # Process each merged group
                for group in merged_groups:
                    if len(group) == 1:
                        # If no merging occurred, add the original row
                        final_rows.append(group[0])
                    else:
                        # Create a merged entry
                        base_entry = group[0].copy()
                        # Find earliest start and latest end positions
                        min_start = min(entry['position_start'] for entry in group)
                        max_end = max(entry['position_end'] for entry in group)
                        # Create new name incorporating the span
                        base_entry['TE_name'] = f"{base_entry['TE_name']}_span_{min_start}_{max_end}"
                        # Update positions
                        base_entry['position_start'] = min_start
                        base_entry['position_end'] = max_end
                        final_rows.append(base_entry)
            
            # Create a new Polars DataFrame from final_rows and add it to final_dfs
            final_df = pl.DataFrame(final_rows)
            final_dfs.append(final_df)
        
        # Concatenate all final DataFrames
        merged_final_df = pl.concat(final_dfs)
        return merged_final_df
    
    @staticmethod
    def remove_duplicates(df):
        # drop window column
        df = df.drop("Window")
        # remove duplicates
        df = df.unique(subset=["chromosome", "position_start", "position_end", "TE_name", "Genome_code"])
        
        return df
    
    @staticmethod
    def remove_duplicates_final(df):
        # remove duplicates
        df = df.unique(subset=["chromosome", "position_start", "position_end", "TE_name"])
        # keep sorted by chromosome and position start and end
        df = df.sort(["chromosome", "position_start", "position_end"])
        return df


    @staticmethod
    def calculate_te_frequency_ref(df, number_of_genomes):
        """
        Calculate the frequency of Transposable Elements (TEs) that are present in a specified number of genomes.

        Parameters:
        df (pd.DataFrame): The input DataFrame containing columns 'chromosome', 'position_start', 
                            'position_end', 'TE_name', and 'Genome_code'.
        number_of_genomes (int): Unused parameter (kept for consistency with function signature).

        Returns:
        pd.DataFrame: A DataFrame with 'chromosome', 'position_start', 'position_end', 'TE_name', 
                    'Genome_code', and 'Frequency_of_TE'.
        """
        df = df.to_pandas()
        
        # Step 1: Group by 'chromosome', 'position_start', 'position_end', and 'TE_name' and aggregate 'Genome_code'
        grouped_df = df.groupby(['chromosome', 'position_start', 'position_end', 'TE_name'])['Genome_code'].agg(lambda x: '; '.join(set(x))).reset_index()
        
        # Step 2: Count the frequency of occurrences of each TE
        frequency_df = df.groupby(['chromosome', 'position_start', 'position_end', 'TE_name']).size().reset_index(name='Frequency_of_TE')
        
        # Step 3: Merge the grouped DataFrame with the frequency DataFrame on the grouping columns
        result = grouped_df.merge(frequency_df, on=['chromosome', 'position_start', 'position_end', 'TE_name'])
        
        # Step 4: Sort the result by 'Frequency_of_TE' in descending order
        result = result.sort_values(by='Frequency_of_TE', ascending=False)
        
        return result
    
    @staticmethod
    def make_vcf(df, chromosome):
        # Drop unnecessary columns
        df = df.drop(["score", "data-analysis-ngs-mapper-2", "Window", "country", "location"])

        # Create a genotype column
        df = df.with_columns(pl.lit("1/1").alias("genotype"))

        # Pivot to make genome_IDs into columns
        vcf_df = df.pivot(
            index=["chromosome", "position_start", "position_end", "TE_name", "strand"],
            columns="Genome_code",
            values="genotype",
            aggregate_function="first"
        ).fill_null("0/0")  # Fill missing values with "0/0"

        # get only rows where the chromosome information is the same as the input
        vcf_df = vcf_df.filter(pl.col("chromosome") == chromosome)
        
        return vcf_df

    @staticmethod
    def remove_non_tes(df):
        # remove satellites and simple repeats
        ## Satellite
        satellite_name = 'Satellite'
        ## Simple repeat
        simple_repeat_name = 'Simple_repeat'
        ## buffer
        buffer_name = 'buffer'

        unique_te_names = df['TE_name'].unique().to_list()

        list_of_tes_to_remove = []
        for TE in unique_te_names:
            if satellite_name in TE or simple_repeat_name in TE or buffer_name in TE:
                list_of_tes_to_remove.append(TE)

        # the ~ is to remove the rows that contain the TE names in the list
        df = df.filter(~pl.col("TE_name").is_in(list_of_tes_to_remove))
        return df
    
    @staticmethod
    def remove_specific_te_families(df, list_TEs):
        list_of_te_names = list_TEs['TE'].unique().to_list()

        # first make two dfs - one with the _span_ values and one without
        merged_TEs = df.filter(pl.col('TE_name').str.contains('_span_'))
        nonmerged_TEs = df.filter(~pl.col('TE_name').str.contains('_span_'))

        # remove TEs that are not in the list from the nonmerged TEs
        # without the ~ it keeps the TEs in the list
        df_nonmerged = nonmerged_TEs.filter(pl.col("TE_name").is_in(list_of_te_names))

        # remove TEs that are not in the list from the merged TEs
        ## split the TE_name column
        merged_TEs = merged_TEs.with_columns([pl.col("TE_name").str.split('_span_').list.get(0).alias("TE_name")])
        ## remove the TEs that are not in the list
        df_merged = merged_TEs.filter(pl.col("TE_name").is_in(list_of_te_names))
        ## re add the _span_ to the TE_name column
        df_merged = df_merged.with_columns([pl.concat_str([pl.col("TE_name"), pl.lit("_span_"), pl.col("position_start").cast(pl.Utf8), pl.lit("_"), pl.col("position_end").cast(pl.Utf8)], separator='').alias("TE_name")])
        ## since the position information is in the start and end position I don't need to re add the span option?
        # merge the two dfs
        df = pl.concat([df_nonmerged, df_merged])
        
        return df
    
    @staticmethod
    def merge_vcf_files(vcf_files):
        # Start with None for the merged dataframe
        merged_df = None
        
        # Define the identifier columns
        id_columns = ['chromosome', 'position_start', 'position_end', 'TE_name', 'strand']
        
        for country, vcf in vcf_files.items():
            print(f"Processing {country} VCF")
            
            # If this is the first dataframe, use it as the base
            if merged_df is None:
                merged_df = vcf.clone()
                continue
            
            # Get sample columns from the current vcf
            sample_columns = [col for col in vcf.columns if col not in id_columns]
            
            # Create temporary composite keys for comparison
            merged_df = merged_df.with_columns(
                pl.concat_str([
                    pl.col('chromosome'),
                    pl.col('position_start').cast(pl.Utf8),
                    pl.col('position_end').cast(pl.Utf8),
                    pl.col('TE_name'),
                    pl.col('strand')
                ], separator='-').alias('_temp_key')
            )
            
            vcf = vcf.with_columns(
                pl.concat_str([
                    pl.col('chromosome'),
                    pl.col('position_start').cast(pl.Utf8),
                    pl.col('position_end').cast(pl.Utf8),
                    pl.col('TE_name'), 
                    pl.col('strand')
                ], separator='-').alias('_temp_key')
            )
            
            # Find keys in vcf that don't exist in merged_df
            merged_keys = set(merged_df.get_column('_temp_key').to_list())
            vcf_keys = set(vcf.get_column('_temp_key').to_list())
            new_keys = vcf_keys - merged_keys
            
            # Get new rows that need to be added
            new_rows = vcf.filter(pl.col('_temp_key').is_in(new_keys))
            
            # Add 0/0 values for all existing sample columns in merged_df
            existing_sample_cols = [c for c in merged_df.columns if c not in id_columns and c != '_temp_key']
            if existing_sample_cols:
                default_values = {col: '0/0' for col in existing_sample_cols}
                new_rows = new_rows.with_columns(**{col: pl.lit('0/0') for col in existing_sample_cols})
            
            # Add new rows to merged_df
            merged_df = pl.concat([merged_df, new_rows.select(merged_df.columns)])
            
            # For existing rows, update with new sample data
            # First, create a mapping for quick lookup
            vcf_dict = {row['_temp_key']: row for row in vcf.to_dicts()}
            
            # For each new sample column, add it to merged_df
            for sample_col in sample_columns:
                if sample_col not in merged_df.columns:
                    # Add new column with default 0/0 values
                    merged_df = merged_df.with_columns(pl.lit('0/0').alias(sample_col))
                
                # Update values for existing TEs
                updates = []
                
                for key in merged_keys.intersection(vcf_keys):
                    vcf_value = vcf_dict[key][sample_col]
                    updates.append((key, vcf_value))
                
                # Apply updates for this sample column
                if updates:
                    update_dict = dict(updates)
                    merged_df = merged_df.with_columns(
                        pl.when(pl.col('_temp_key').is_in(list(update_dict.keys())))
                        .then(pl.col('_temp_key').map_dict(update_dict))
                        .otherwise(pl.col(sample_col))
                        .alias(sample_col)
                    )
            
            # Remove temporary key
            merged_df = merged_df.drop('_temp_key')
        
        # sort the df by chromosome and position start and position end
        merged_df = merged_df.sort(['chromosome', 'position_start', 'position_end'])
        return merged_df

    @staticmethod
    def get_genome_codes_no_TE(df):
        """
        Get a dictionary mapping TE names to genome codes where the TE is absent (0/0).

        Args:
            df (pl.DataFrame): The input DataFrame containing TE presence/absence data.

        Returns:
            dict: Dictionary with TE names as keys and lists of genome codes with 0/0 as values.
        """

        # Extract sample columns (all except the first four columns)
        sample_columns = df.columns[4:]

        # Dictionary to store results
        te_absent_dict = {}

        for row in df.iter_rows(named=True):
            te_name = row["TE_name"] + "_" + row["chromosome"] + "_" + str(row["position_start"]) + "_" + str(row["position_end"]) + "_" + row["strand"]
            absent_samples = [sample for sample in sample_columns if row[sample] == "0/0"]

            if absent_samples:
                te_absent_dict[te_name] = absent_samples

        return te_absent_dict
    
    @staticmethod
    def extract_large_insert_reads(bam_file, chromosome, start, end, insert_threshold):
        """Same implementation as before but with early exit optimizations"""

        # how much variation in the distance between the flanking reads should be allowed
        te_size = end - start
        # use insertion size compared to the TE size to check if there is too much variation - tolerance for that could be the insert threshold
        tolerance = insert_threshold

        if bam_file == 'Not Found':
            return './.'

        try:
            with pysam.AlignmentFile(bam_file, "rb") as bam:
                # Check upstream region
                for read in bam.fetch(chromosome, max(0, start - 1000), start):
                    read_end = read.reference_start + read.query_length
                    mate_start = read.next_reference_start
                    mate_end = read.next_reference_start + read.query_length
                    insert_size = abs(mate_end - read.reference_start)
                    x = insert_size - 2 * read.query_length

                    if read.flag in {99, 163} and read.template_length > insert_threshold and mate_start > end and mate_end < (end + 1000) and abs(x - te_size) <= tolerance:
                        return '0/0'  # Early exit once we find a qualifying read
                
                # Check downstream region
                for read in bam.fetch(chromosome, end, end + 1000):
                    read_end = read.reference_start + read.query_length
                    mate_start = read.next_reference_start
                    mate_end = read.next_reference_start + read.query_length
                    insert_size = abs(mate_end - read.reference_start)
                    x = insert_size - 2 * read.query_length

                    if read.flag in {99, 163} and read.template_length > insert_threshold and mate_end < start and mate_start > (start - 1000) and abs(x - te_size) <= tolerance:
                        return '0/0'  # Early exit once we find a qualifying read
                        
        except Exception as e:
            return './.'  # Error case
            
        return './.'  # No qualifying reads found
    
    @staticmethod
    def process_te_batch(te_batch, bam_files, cutoff_values):
        """Process a batch of TEs for a specific BAM file"""
        results = {}
        
        for te_info in te_batch:
            TE_name, genome, chromosome, start, end = te_info
            bam_file = bam_files.get(genome, 'Not Found')
            insert_threshold = int(cutoff_values.get(genome, 1000))
            
            te_name_genome = f'{TE_name}_{genome}'
            results[te_name_genome] = Making_vcf.extract_large_insert_reads(
                bam_file, chromosome, start, end,insert_threshold
            )
            
        return results

    @staticmethod
    def check_bam_files(dict_absent_tes, bam_files, cutoff_values, num_processes=None):
        """
        Parallelized version of check_bam_files
        
        Args:
            dict_absent_tes (dict): {TE_name: [list of genome codes where TE is absent]}
            bam_files (dict): {genome_code: bam_file_path}
            cutoff_values (dict): {genome_code: insert_size_threshold}
            num_processes (int): Number of processes to use. Default is CPU count.
            
        Returns:
            dict: {genome_code: '0/0' or './.'}
        """
        if num_processes is None:
            num_processes = mp.cpu_count()
        
        # Group TEs by genome to reduce BAM file opening/closing
        genome_te_batches = defaultdict(list)
        
        for TE_name, absent_genomes in dict_absent_tes.items():
            match = re.search(r'(.+)_([A-Z].+)_(\d+)_(\d+)_(.+)_(.+)', TE_name)
            if not match:
                print(f'No match found for {TE_name}')
                continue
                
            chromosome, start, end = match.group(2), int(match.group(3)), int(match.group(4))
            
            for genome in absent_genomes:
                genome_te_batches[genome].append((TE_name, genome, chromosome, start, end))
        
        # Create a flat list of all TE batches
        all_batches = []
        for genome, te_list in genome_te_batches.items():
            # Further batch by chromosome to minimize BAM file operations
            chrom_batches = defaultdict(list)
            for te_info in te_list:
                chrom_batches[te_info[2]].append(te_info)
                
            for chrom_batch in chrom_batches.values():
                all_batches.append(chrom_batch)
        
        # Use multiprocessing to process the batches
        results = {}
        with mp.Pool(processes=num_processes) as pool:
            batch_results = pool.map(
                partial(Making_vcf.process_te_batch, bam_files=bam_files, cutoff_values=cutoff_values),
                all_batches
            )
            
            # Merge all batch results
            for batch_result in batch_results:
                results.update(batch_result)
                
        return results
                
    @staticmethod
    def find_insert_size_cutoff(path):
        list_tables = os.listdir(path)
        dict_cutoff_values = {}
        for csv in list_tables:
            genome_code = re.search('insert_sizes_(.+).csv', csv).group(1)
            if csv.endswith('.csv'):
                df = pd.read_csv(os.path.join(path, csv))
                df = df[df['Cumulative Percentage'] >= 95]
                insert_size_cutoff = df.iloc[0]['Insert Size']
                dict_cutoff_values[genome_code] = insert_size_cutoff

        return dict_cutoff_values
    
    @staticmethod
    def update_vcf(vcf_file, updated_genotypes):
        """
        Update the VCF DataFrame with new genotypes.

        Args:
            vcf_file (DataFrame): VCF-style DataFrame with TE information.
            updated_genotypes (dict): {TE_name_genome: genotype} where TE_name_genome = 'TE_chr_start_end_genome'

        Returns:
            DataFrame: Updated VCF DataFrame.
        """

        # convert to pandas
        vcf_file = vcf_file.to_pandas()

        for te_info, genotype in updated_genotypes.items():
            # Match pattern like 'TE_chr_start_end_genome'
            match = re.search(r'(.+)_([A-Z].+)_(\d+)_(\d+)_([-+])_(.+)', te_info)
            if not match:
                print(f"Skipping invalid TE info format: {te_info}")
                continue

            te_name, chromosome, start, end, strand, genome_code = match.groups()
            start, end = int(start), int(end)  # Convert to integers

            # Find the row index
            row_match = vcf_file[
                (vcf_file['chromosome'] == chromosome) &
                (vcf_file['position_start'] == start) &
                (vcf_file['position_end'] == end) &
                (vcf_file['TE_name'] == te_name) &
                (vcf_file['strand'] == strand)
            ]

            if row_match.empty:
                print(f"No match found for {te_info}")
                continue

            # Get the index of the matching row
            row_index = row_match.index[0]

            # Ensure the genome code column exists
            if genome_code not in vcf_file.columns:
                print(f"Genome code column {genome_code} not found in VCF file.")
                continue

            # Update the genotype
            vcf_file.at[row_index, genome_code] = genotype

        return vcf_file  # Return updated DataFrame
    
    @staticmethod
    def remove_fixed_almostfixed_TEs(df):
        # Select only genome columns (exclude metadata columns)
        df_genomes = df.drop(['chromosome', 'position_start', 'position_end', 'TE_name', 'strand'], axis='columns')
        
        # Count occurrences of '0/0' per row - if 0/0 is not present the TE is fixed or almost fixed
        count_0_0 = (df_genomes == '0/0').sum(axis=1)
        
        # Keep rows where at least two '0/0' are present
        df = df[count_0_0 >= 2]
        
        return df
    
    @staticmethod
    def split_vcf_into_windows(df, window_size):
        """
        Split the VCF DataFrame into windows of a specified size, assigning each TE to
        only the first window it overlaps.

        Args:
            df (pl.DataFrame): The input VCF DataFrame (Polars).
            window_size (int): The size of each window.

        Returns:
            pl.DataFrame: The VCF DataFrame split into windows with a window_id column.
        """
        chromosomes = df['chromosome'].unique()
        window_records = []

        for chrom in chromosomes:
            chrom_df = df.filter(pl.col('chromosome') == chrom)
            max_pos = chrom_df['position_end'].max()
            
            for row in chrom_df.iter_rows(named=True):
                te_start = row['position_start']
                te_end = row['position_end']

                # Iterate over windows and stop at the first overlapping one
                for window_id, start_pos in enumerate(range(0, max_pos, window_size)):
                    end_pos = start_pos + window_size

                    if te_start < end_pos and te_end > start_pos:
                        # TE overlaps this window; assign it and break
                        row_with_window = dict(row)
                        row_with_window['window_id'] = window_id
                        window_records.append(row_with_window)
                        break  # Only assign to first overlapping window

        return pl.DataFrame(window_records)