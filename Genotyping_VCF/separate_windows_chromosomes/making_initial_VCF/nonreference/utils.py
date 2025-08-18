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
    
    @staticmethod
    def make_nonref_retroseq_df(list_of_bed, chromosome):
        """
        Creates a dictionary of Polars DataFrames from a list of BED file paths.

        This function reads each BED file from the provided list of file paths, converts it into a Polars DataFrame, 
        and processes the DataFrame by renaming columns, extracting Transposable Element (TE) names, 
        and filtering by the specified chromosome. The processed DataFrames are stored in a dictionary with the file paths as keys.

        Args:
            list_of_bed (list of str): A list of file paths to BED files containing TE insertion data.
            chromosome (str): The chromosome to filter the rows by.

        Returns:
            dict: A dictionary where each key is a BED file path and each value is a Polars DataFrame 
                containing the processed TE insertion data filtered by chromosome.
        """
        non_ref_dict_polars = {}

        for file in list_of_bed:
            general_df = pl.read_csv(file, separator='\t', skip_rows=1, has_header=False)
            general_df = general_df.rename({"column_1": "chromosome", 
                                            "column_2": "position_start", 
                                            "column_3": "position_end", 
                                            "column_4": "TE_name", 
                                            "column_5": "score", 
                                            "column_6": "strand"}) 

            filtered_df = general_df.filter(pl.col("chromosome") == chromosome)

            non_result = filtered_df.with_columns([
                pl.col("TE_name").str.split('|').list.get(0).alias("TE_name")
            ])

            non_ref_dict_polars[file] = non_result

        return non_ref_dict_polars

    @staticmethod
    def check_breakpoints_retroseq(non_ref_retroseq):
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
        
        merged_df = pl.concat(all_unique_combinations).unique()

        updated_dfs = {}

        for genome_name, df in non_ref_retroseq.items():
            genome_combinations = df.select(
                pl.col("chromosome"),
                pl.col("position_start"),
                pl.col("position_end"),
                pl.col("TE_name")
            ).unique()
            
            missing_insertions = merged_df.join(genome_combinations, on=["chromosome", "position_start", "position_end", "TE_name"], how="anti")
            
            match = re.search(r'(.+/results/retroseq/edited_results)/.+', genome_name)
            if match:
                general_path = match.group(1)
            
                breakpoint_1_files = glob.glob(f'{general_path}/breakpoint_1/*nonredundant.bed')
                
                if breakpoint_1_files:
                    breakpoint_1 = breakpoint_1_files[0]
                    print(f"Found breakpoint_1 file: {breakpoint_1}")
                    general_df = pl.read_csv(breakpoint_1, separator='\t', skip_rows=1, has_header=False)
                    
                    general_df = general_df.rename({"column_1": "chromosome", "column_2": "position_start", 
                                                    "column_3": "position_end", "column_4": "TE_name", 
                                                    "column_5": "score", "column_6": "strand"}) 
                    general_df = general_df.drop(["score", "strand"])
                    
                    non_result = general_df.with_columns([pl.col("TE_name").str.split('|').list.get(0).alias("TE_name")])
                    
                    matching_insertions = missing_insertions.join(non_result, on=["chromosome", "position_start", "position_end", "TE_name"], how="inner")
                    
                    updated_df = df.vstack(matching_insertions)
                    updated_dfs[genome_name] = updated_df

                else:
                    print(f'No "nonredundant.bed" file found in the breakpoint_1 directory. - {genome_name}')
            else:
                print(f"Unable to extract general path from genome name: {genome_name}")

        return updated_dfs
    
    @staticmethod
    def check_for_tes_between_individuals(dict_individuals, variation):
        df_list = []
        
        for genome_name, df in dict_individuals.items():
            genome_code = genome_name 
            df = df.with_columns([pl.lit(genome_code).alias('genome_code')])
            df_list.append(df)
        
        big_df = pl.concat(df_list, how="vertical").sort(['TE_name', 'chromosome', 'position_start'])
        
        te_information_list = []
        updated_positions = []
        num_rows = big_df.height
        i = 0
        
        while i < num_rows:
            info_i = big_df.row(i, named=True)
            chromosome_i = info_i['chromosome']
            te_name_i = info_i['TE_name']
            
            min_start_position = info_i['position_start']
            max_end_position = info_i['position_end']
            
            j = i + 1
            while j < num_rows:
                info_j = big_df.row(j, named=True)
                chromosome_j = info_j['chromosome']
                te_name_j = info_j['TE_name']
                
                if te_name_i == te_name_j and chromosome_i == chromosome_j and (
                    abs(info_j['position_start'] - min_start_position) <= variation or
                    abs(info_j['position_end'] - max_end_position) <= variation):
                    
                    min_start_position = min(min_start_position, info_j['position_start'])
                    max_end_position = max(max_end_position, info_j['position_end'])
                    j += 1
                else:
                    break
            
            chosen_te_information = f"{chromosome_i}_{min_start_position}_{max_end_position}_{te_name_i}"
            
            for k in range(i, j):
                te_information_list.append(chosen_te_information)
                updated_positions.append((min_start_position, max_end_position))
            
            i = j
        
        big_df = big_df.with_columns([
            pl.Series([x[0] for x in updated_positions]).alias("position_start"),
            pl.Series([x[1] for x in updated_positions]).alias("position_end"),
            pl.Series(te_information_list).alias("TE_information")
        ])
        
        result_dict = {}
        
        for group_tuple in big_df.group_by(['genome_code']):
            genome_code, df_group = group_tuple
            genome_code = genome_code[0]
            df_group = df_group.select(['chromosome', 'position_start', 'position_end', 'TE_name', 'score', 'strand', 'TE_information'])
            result_dict[genome_code] = df_group
        
        return result_dict
    
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
        df = location.clone()

        country_dfs = {}

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
    def making_windows_for_comparison(df_per_te: pl.DataFrame) -> pl.DataFrame:
        """
        Adds a 'Window' column to indicate which genomic window(s) each TE overlaps.

        For TEs that overlap multiple windows, a separate row is created for each window.

        Args:
            df_per_te (pl.DataFrame): DataFrame with 'chromosome', 'position_start', 'position_end'.

        Returns:
            pl.DataFrame: Expanded DataFrame with 'Window' column indicating the window start position.
        """
        window_size = 100000
        
        chrm_list = df_per_te.select(pl.col("chromosome").unique()).to_series().to_list()

        windowed_rows = []

        for chrm in chrm_list:
            df_chr = df_per_te.filter(pl.col("chromosome") == chrm)
            max_position = df_chr["position_end"].max()

            for start_position in range(0, max_position + 1, window_size):
                end_position = start_position + window_size - 1

                overlapping = df_chr.filter(
                    (pl.col("position_end") >= start_position) & (pl.col("position_start") <= end_position)
                )

                for row in overlapping.iter_rows(named=True):
                    row_copy = row.copy()
                    row_copy["Window"] = start_position
                    windowed_rows.append(row_copy)

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
            
            for list_info, group in group_df.group_by(group_cols):
                num_rows = group.height
                
                merged_groups = []
                processed_indices = set()
                
                for i in range(num_rows):
                    if i in processed_indices:
                        continue
                        
                    current_group = []
                    info_i = group.row(i, named=True)
                    current_group.append(info_i)
                    processed_indices.add(i)
                    
                    for j in range(i + 1, num_rows):
                        if j in processed_indices:
                            continue
                            
                        info_j = group.row(j, named=True)
                        
                        if (abs(info_j['position_start'] - info_i['position_start']) <= variation or
                            abs(info_i['position_end'] - info_j['position_end']) <= variation) and \
                            info_i['TE_name'] == info_j['TE_name']:
                            current_group.append(info_j)
                            processed_indices.add(j)
                    
                    if current_group:
                        merged_groups.append(current_group)
                
                for group in merged_groups:
                    if len(group) == 1:
                        final_rows.append(group[0])
                    else:
                        base_entry = group[0].copy()
                        min_start = min(entry['position_start'] for entry in group)
                        max_end = max(entry['position_end'] for entry in group)
                        base_entry['TE_name'] = f"{base_entry['TE_name']}_span_{min_start}_{max_end}"
                        base_entry['position_start'] = min_start
                        base_entry['position_end'] = max_end
                        final_rows.append(base_entry)
            
            final_df = pl.DataFrame(final_rows)
            final_dfs.append(final_df)
        
        merged_final_df = pl.concat(final_dfs)
        return merged_final_df
    
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
        
        grouped_df = df.groupby(['chromosome', 'position_start', 'position_end', 'TE_name'])['Genome_code'].agg(lambda x: '; '.join(set(x))).reset_index()
        
        frequency_df = df.groupby(['chromosome', 'position_start', 'position_end', 'TE_name']).size().reset_index(name='Frequency_of_TE')
        
        result = grouped_df.merge(frequency_df, on=['chromosome', 'position_start', 'position_end', 'TE_name'])
        
        result = result.sort_values(by='Frequency_of_TE', ascending=False)
        
        return result
    
    @staticmethod
    def make_vcf(df, chromosome):
        df = df.drop(["score", "data-analysis-ngs-mapper-2", "Window", "country", "location"])

        df = df.with_columns(pl.lit("1/1").alias("genotype"))

        vcf_df = df.pivot(
            index=["chromosome", "position_start", "position_end", "TE_name", "strand"],
            columns="Genome_code",
            values="genotype",
            aggregate_function="first"
        ).fill_null("0/0")  # Fill missing values with "0/0"

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

        # without the ~ it keeps the TEs in the list
        df = df.filter(pl.col("TE_name").is_in(list_of_te_names))
        return df
    
    @staticmethod
    def remove_singletons(df):
        if df is None:
            raise ValueError("No non-reference TEs in this chromsome")
        
        meta_cols = ['chromosome', 'position_start', 'position_end', 'TE_name', 'strand']
        
        genome_cols = [col for col in df.columns if col not in meta_cols]
        
        if not genome_cols:
            return df

        count_expr = sum(pl.col(col) == '1/1' for col in genome_cols).alias('count_1_1')
        df_with_count = df.with_columns(count_expr)
        
        filtered_df = df_with_count.filter(pl.col('count_1_1') != 1)
        
        return filtered_df.drop('count_1_1')
    
    @staticmethod
    def merge_vcf_files(vcf_files):
        merged_df = None
        
        id_columns = ['chromosome', 'position_start', 'position_end', 'TE_name', 'strand']
        
        for country, vcf in vcf_files.items():
            print(f"Processing {country} VCF")
            
            if merged_df is None:
                merged_df = vcf.clone()
                continue
            
            sample_columns = [col for col in vcf.columns if col not in id_columns]
            
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
            
            merged_keys = set(merged_df.get_column('_temp_key').to_list())
            vcf_keys = set(vcf.get_column('_temp_key').to_list())
            new_keys = vcf_keys - merged_keys
            
            new_rows = vcf.filter(pl.col('_temp_key').is_in(new_keys))
            
            existing_sample_cols = [c for c in merged_df.columns if c not in id_columns and c != '_temp_key']
            if existing_sample_cols:
                default_values = {col: '0/0' for col in existing_sample_cols}
                new_rows = new_rows.with_columns(**{col: pl.lit('0/0') for col in existing_sample_cols})
            
            merged_df = pl.concat([merged_df, new_rows.select(merged_df.columns)])
            
            vcf_dict = {row['_temp_key']: row for row in vcf.to_dicts()}
            
            for sample_col in sample_columns:
                if sample_col not in merged_df.columns:
                    merged_df = merged_df.with_columns(pl.lit('0/0').alias(sample_col))
                
                updates = []
                
                for key in merged_keys.intersection(vcf_keys):
                    vcf_value = vcf_dict[key][sample_col]
                    updates.append((key, vcf_value))
                
                if updates:
                    update_dict = dict(updates)
                    merged_df = merged_df.with_columns(
                        pl.when(pl.col('_temp_key').is_in(list(update_dict.keys())))
                        .then(pl.col('_temp_key').map_dict(update_dict))
                        .otherwise(pl.col(sample_col))
                        .alias(sample_col)
                    )
            
            merged_df = merged_df.drop('_temp_key')
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
        sample_columns = df.columns[4:]

        te_absent_dict = {}

        for row in df.iter_rows(named=True):
            te_name = row["TE_name"] + "_" + row["chromosome"] + "_" + str(row["position_start"]) + "_" + str(row["position_end"])
            absent_samples = [sample for sample in sample_columns if row[sample] == "0/0"]

            if absent_samples:
                te_absent_dict[te_name] = absent_samples

        return te_absent_dict
    
    @staticmethod
    def extract_large_insert_reads(bam_file, chromosome, start, end, insert_threshold):
        """
        Check for large insert reads in a BAM file.

        Args:
            bam_file (str): Path to the BAM file.
            chromosome (str): Chromosome name.
            start (int): Start position of TE.
            end (int): End position of TE.
            insert_threshold (int): Minimum insert size to consider.

        Returns:
            str: '0/0' if large insert reads are found, './.' otherwise.
        """
        if bam_file == 'Not Found':
            return './.'

        try:
            with pysam.AlignmentFile(bam_file, "rb") as bam:
                for read in bam.fetch(chromosome, max(0, start - 1000), start):
                    if read.flag in {99, 163} and read.template_length > insert_threshold:
                        return '0/0'  # Large insert read found

                for read in bam.fetch(chromosome, end, end + 1000):
                    if read.flag in {99, 163} and read.template_length > insert_threshold:
                        return '0/0'  # Large insert read found

        except Exception as e:
            return './.'  # Error case, treat as missing data

        return './.'  # No large insert reads found

    @staticmethod
    def check_bam_files(dict_absent_tes, bam_files, cutoff_values):
        """
        Check BAM files and determine TE genotype ('0/0' if large insert reads are found, './.' if not).

        Args:
            dict_absent_tes (dict): {TE_name: [list of genome codes where TE is absent]}
            bam_files (dict): {genome_code: bam_file_path}
            cutoff_values (dict): {genome_code: insert_size_threshold}

        Returns:
            dict: {genome_code: '0/0' or './.'}
        """
        genotype_results = {}

        for TE_name, absent_genomes in dict_absent_tes.items():
            match = re.search(r'(.+)_([A-Z].+)_(\d+)_(\d+)', TE_name)
            if not match:
                print('no match found')
                continue  # Skip if TE name format is incorrect

            chromosome, start, end = match.group(2), int(match.group(3)), int(match.group(4))

            for genome in absent_genomes:
                bam_file = bam_files.get(genome, 'Not Found')
                insert_threshold = int(cutoff_values.get(genome, 1000))
                
                te_name_genome = f'{TE_name}_{genome}'

                # Get genotype from BAM file
                genotype_results[te_name_genome] = Making_vcf.extract_large_insert_reads(
                    bam_file, chromosome, start, end, insert_threshold
                )

        return genotype_results
                
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
            match = re.search(r'(.+)_([A-Z].+)_(\d+)_(\d+)_([A-Z]\w+)', te_info)
            if not match:
                print(f"Skipping invalid TE info format: {te_info}")
                continue

            te_name, chromosome, start, end, genome_code = match.groups()
            start, end = int(start), int(end)  # Convert to integers

            # Find the row index
            row_match = vcf_file[
                (vcf_file['chromosome'] == chromosome) &
                (vcf_file['position_start'] == start) &
                (vcf_file['position_end'] == end) &
                (vcf_file['TE_name'] == te_name)
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
