import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from itertools import chain
import pysam
from ld_estimator import pairwise_ld
import polars as pl
import sys
import glob

class FilesInterpretation:
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

        dfs = []

        for i, country in enumerate(countries):

            country_cols = [col for col in df.columns 
                        if col not in base_cols and col.endswith('_' + country)]
            
            if not country_cols:
                continue

            country_df = df[base_cols + country_cols]
            
            country_df = country_df.dropna(subset=country_cols, how='all')

            if country == 'Colombia':
                return country_df

class LD_analysis:
    @staticmethod
    def find_variant_information(file_path, variant_info):
        variant_info_list = []

        vcf_file = pysam.VariantFile(file_path)

        for variant in variant_info:
            name = variant['variant_name']
            chromosome = variant['Chromosome']
            position_start = variant['Position_start']

            for record in vcf_file:
                if record.pos == position_start:
                    genotype_info = {'Variant': name}
                    for sample in record.samples:
                        genotype = record.samples[sample]['GT']
                        genotype_info[sample] = "/".join(map(str, genotype)) if genotype is not None else "./."

                    variant_info_list.append(genotype_info)

                    break

        vcf_file.close()

        df = pd.DataFrame(variant_info_list)

        df = df.set_index('Variant')

        return df

    @staticmethod
    def add_country_name_SNPs(df, file_loc): 
        info_location_country = pd.read_csv(file_loc)

        info_location_country['space'] = info_location_country['country']
        
        match_list = df.columns

        genome_space_dict = pd.Series(info_location_country.space.values, index=info_location_country.Genome_code).to_dict()

        new_columns = []
        for col in match_list:
            if col in genome_space_dict:
                space_info = genome_space_dict[col]
                new_columns.append(f"{col}_{space_info}")
            else:
                new_columns.append(col)
        
        df.columns = new_columns
        return df

    @staticmethod
    def remove_heterozygous_for_SNPs(general_df):
        columns_to_keep = [
            col for col in general_df.columns 
            if general_df[col].isin(['0/0', '1/1']).all()
        ]
        
        cleaned_df = general_df[columns_to_keep]
        
        return cleaned_df

    def get_gene_region(df, gene_info):
        gene_name = gene_info['gene_name']
        chromosome = gene_info['Chr']
        start = gene_info['position_start']
        end = gene_info['position_end']

        # Filter the dataframe
        filtered_df = df[(df['chromosome'] == chromosome) & 
                    (df['position_start'] >= start) & 
                    (df['position_end'] <= end)]
                                
        return filtered_df

    @staticmethod
    def remove_fixed_TEs_country(df):
        df['Country'] = df['Genome_code'].str.split('_').str[-1]
        
        unique_counts = df.groupby('Country')['Genome_code'].nunique()
        
        print("Number of genomes for each country:")
        print(unique_counts)
        
        df['number_of_genomes'] = df['Country'].map(unique_counts)
        
        grouped_df = df.groupby(['chromosome', 'position_start', 'position_end', 'TE_name', 'Country']).agg({
            'Genome_code': lambda x: '; '.join(set(x)),
            'number_of_genomes': 'first'
        }).reset_index()
        
        grouped_df['Frequency_of_TE'] = df.groupby(['chromosome', 'position_start', 'position_end', 'TE_name', 'Country']).size().values
        fixed_TEs = grouped_df[grouped_df['Frequency_of_TE'] == grouped_df['number_of_genomes']]
        
        fixed_TE_tuples = set(fixed_TEs[['position_start', 'position_end', 'TE_name']].itertuples(index=False, name=None))
        
        filtered_df = df[~df[['position_start', 'position_end', 'TE_name']].apply(tuple, axis=1).isin(fixed_TE_tuples)]
        
        return filtered_df

    @staticmethod
    def make_list_of_te_insertions(gene_df):

        selected_columns_df = gene_df[['chromosome', 'position_start', 'position_end', 'TE_name', 'Genome_code']]

        selected_columns_df.loc[:, 'TE_name'] = selected_columns_df['TE_name'].str.split('|').str[0]
        
        grouped_df = selected_columns_df.groupby(['chromosome', 'position_start', 'position_end', 'TE_name'])['Genome_code'].apply(lambda x: ';'.join(x)).reset_index()
        
        return grouped_df
    
    @staticmethod
    def make_table_insertion(unique_df):
        unique_df['TE_identifier'] = unique_df['chromosome'] + '_' + unique_df['position_start'].astype(str) + '_' + unique_df['position_end'].astype(str) + '_' + unique_df['TE_name']

        unique_df['Genome_code'] = unique_df['Genome_code'].str.split(';')

        all_genome_codes = sorted(set(code for codes in unique_df['Genome_code'] for code in codes))

        presence_dict = {}
        for idx, row in unique_df.iterrows():
            te_id = row['TE_identifier']
            for code in all_genome_codes:
                if code in row['Genome_code']:
                    presence_dict.setdefault(te_id, {})[code] = (1, 1)
                else:
                    presence_dict.setdefault(te_id, {})[code] = (0, 0)
        
        presence_df = pd.DataFrame(presence_dict).T.reset_index()
        presence_df.rename(columns={'index': 'TE_identifier'}, inplace=True)

        return presence_df

    @staticmethod
    def merge_with_te_table(df_snp, df_te):
        df_snp = df_snp.reset_index()
        df_snp = df_snp.rename(columns={'Variant': 'Variant_name'})  # Rename if needed
        
        df_te = df_te.rename(columns={'TE_identifier': 'Variant_name'})
        
        common_columns = df_snp.columns.intersection(df_te.columns)
        print(df_te.columns)
        print(common_columns)
        
        df_snp = df_snp[common_columns]
        df_te = df_te[common_columns]
        
        merged_df = pd.concat([df_snp, df_te], ignore_index=True)
        
        return merged_df

    @staticmethod
    def merge_presence_with_variants(presence_table_nonref, all_variants_df):
        df = all_variants_df.reset_index(names=['Variant_name'])
        df.columns = [df.columns[0]] + [col + '_Colombia' for col in df.columns[1:]]
        
        matching_columns = [col for col in presence_table_nonref.columns if col in df.columns]
        
        merged_df = pd.concat(
            [presence_table_nonref[matching_columns], df[matching_columns]],
            ignore_index=True
        )

        return merged_df
    
    @staticmethod
    def remove_columns_with_none_pairs(df):
        columns_to_remove = []

        for col in df.columns:
            if any(df[col] == 'None/None'):
                columns_to_remove.append(col)

        df = df.drop(columns=columns_to_remove)

        return df
    
    @staticmethod
    def get_encoding_SNPs(df):
        variants = [
            'AaegL5_3_316080722_V410L',
            'AaegL5_3_315983763_V1016I',
            'AaegL5_3_315999297_I915K',
            'AaegL5_3_316014588_S723T']
        
        variant_rows = df[df['Variant_name'].isin(variants)].copy()

        other_rows = df[~df['Variant_name'].isin(variants)].copy()
        df_encoding_1 = variant_rows.copy()
        df_encoding_2 = variant_rows.copy()

        df_encoding_1['Variant_name'] = df_encoding_1['Variant_name'] + '_encoding_1'
        df_encoding_2['Variant_name'] = df_encoding_2['Variant_name'] + '_encoding_2'

        df_encoded = pd.concat([df_encoding_1, df_encoding_2], ignore_index=True)

        sample_cols = [col for col in df_encoded.columns if col != 'Variant_name']

        df_encoded.loc[df_encoded['Variant_name'].str.contains('_encoding_1'), sample_cols] = \
            df_encoded.loc[df_encoded['Variant_name'].str.contains('_encoding_1'), sample_cols].replace({'0/1': '1/1', '1/0': '1/1'})

        df_encoded.loc[df_encoded['Variant_name'].str.contains('_encoding_2'), sample_cols] = \
            df_encoded.loc[df_encoded['Variant_name'].str.contains('_encoding_2'), sample_cols].replace({'0/1': '0/0', '1/0': '0/0'})

        final_df = pd.concat([other_rows, df_encoded], ignore_index=True)
        return final_df

    @staticmethod
    def calcule_LD_equation(combination_1, combination_2):
        total_count_combination1 = len(combination_1)
        total_count_combination2 = len(combination_2)
        
        countA = combination_1.count("1/1")
        pA = countA / total_count_combination1

        countB = combination_2.count("1/1")
        pB = countB / total_count_combination2

        match_count = sum(1 for c1, c2 in zip(combination_1, combination_2) if c1 == "1/1" and c2 == "1/1")
        pAB = match_count / total_count_combination1

        # Calculating D
        D = pAB - (pA * pB)

        # Calculating D_max
        if pA * pB == 0 or (1 - pA) * (1 - pB) == 0:
            Dmax = 0
        elif D > 0:
            Dmax = min(pA * (1 - pB), (1 - pA) * pB)
        else:
            Dmax = min(pA * pB, (1 - pA) * (1 - pB))

        # Calculating D'
        D_prime = None if Dmax == 0 else D / Dmax

        # Calculating R_squared
        denominator = pA * (1 - pA) * pB * (1 - pB)
        r_squared = None if denominator == 0 else (D ** 2) / denominator

        return D_prime, r_squared

    @staticmethod
    def calculate_LD_with_estimator(presence_table, significance_threshold):
        # Extract TE_identifier and variant columns
        te_identifiers = presence_table['Variant_name'].values
        variants = presence_table.drop(columns='Variant_name').values
        
        num_columns_minus_one = len(presence_table.columns) - 1
        is_haploid = [False] * num_columns_minus_one
        
        results = {}
        genotypes = {}
        
        # Iterate over each pair of rows to calculate LD
        for i in range(len(te_identifiers)):
            for j in range(i + 1, len(te_identifiers)):
                te_id_1 = te_identifiers[i]
                te_id_2 = te_identifiers[j]
                variants_1 = variants[i]
                variants_2 = variants[j]
        
                # Create the combinations
                combination_1 = list(variants_1)
                combination_2 = list(variants_2)

                # Define TE pair
                te_pair = (te_id_1, te_id_2)
                
                # Calculating LD
                #ld = pairwise_ld(combination_1, combination_2, is_haploid)
                D_prime, r_square = LD_analysis.calcule_LD_equation(combination_1, combination_2)
                outfile = open('r_squaredata_d_prime.txt', 'a')
                outfile.write(f'{te_pair} - D_prime: {D_prime} - r_square: {r_square} - combination_1: {combination_1} - combination_2: {combination_2}\n')
                
                # Ensure r_square is not None before comparing
                if r_square is not None and r_square > significance_threshold and ('V410L' in te_id_1 or 'V410L' in te_id_2):
                    print(f'{te_pair} - D_prime: {D_prime} - r_square: {r_square} - combination_1: {combination_1} - combination_2: {combination_2}')
                    genotypes[te_pair] = {'D': D_prime, 'r': r_square, 'combination_1': combination_1, 'combination_2': combination_2}

                results[te_pair] = {'D': D_prime, 'r': r_square}
                
        return results, genotypes

    @staticmethod
    def make_interpret_table(LD_values):
        te_identifiers = set()
        for te_pair in LD_values.keys():
            te_identifiers.update(te_pair)
        
        te_identifiers = sorted(te_identifiers) 
        
        ld_table = pd.DataFrame(index=te_identifiers, columns=te_identifiers)
        
        for te_pair, values in LD_values.items():
            te1, te2 = te_pair
            ld_table.loc[te1, te2] = values['r']
            ld_table.loc[te2, te1] = values['r']
        
        np.fill_diagonal(ld_table.values, 1.0)
        
        return ld_table
    
    @staticmethod  
    def retrieve_significant_rsquared(df, significant_threshold, output_filename):
        # Get the row names (indexes)
        row_names = df.index

        # Extract both encoding columns
        col1 = 'AaegL5_3_316080722_V410L_encoding_1'
        col2 = 'AaegL5_3_316080722_V410L_encoding_2'

        if col1 not in df.columns or col2 not in df.columns:
            raise ValueError(f"One or both columns {col1} and {col2} not found in the DataFrame.")

        # Create a new DataFrame with TEs and both encodings
        selected_columns = pd.DataFrame({
            'TEs': row_names,
            col1: df[col1],
            col2: df[col2]
        })

        filtered_columns = selected_columns[
            (selected_columns[col1] > significant_threshold) |
            (selected_columns[col2] > significant_threshold)
        ]

        print(f'Filtered rows with values in {col1} or {col2} > {significant_threshold}:')
        print(filtered_columns)

        filtered_columns.to_csv(output_filename, index=False)

        te_list = filtered_columns['TEs'].tolist()
        print(f'List of TEs with significant values: {te_list}')

        return filtered_columns, te_list

    @staticmethod 
    def plot_ld_heatmap_r2(ld_table, list_TEs, output_file, highlight_variant='AaegL5_3_316080722_V410L'):
        import pandas as pd
        import numpy as np
        import matplotlib.pyplot as plt
        import seaborn as sns
        """
        Plot a heatmap for r^2 values from the lower triangle of the LD table.
        Ensures all data is plotted while only showing labels for variants in `list_TEs`.
        Highlights a specific variant name in red.
        """
        numeric_table = ld_table.apply(pd.to_numeric, errors='coerce')
        
        mask = np.triu(np.ones_like(numeric_table, dtype=bool), k=1)
        
        lower_triangle_table = numeric_table.where(mask == False)
        
        plt.figure(figsize=(20, 10), dpi=300)
        
        sns.heatmap(lower_triangle_table, 
                    cmap='Reds',  # Blue-green color palette
                    annot=False,    # Don't show individual values
                    cbar_kws={'label': r'$r^2$ Values'},
                    square=True,    # Make cells square
                    vmin=0,         # Minimum value for color scaling
                    vmax=1,         # Maximum value for color scaling
                    mask=lower_triangle_table.isna())  # Mask NaN values
        
        plt.title(r'$r^2$ Linkage Disequilibrium Heatmap', fontsize=60)
        plt.xlabel('Variants', fontsize=30)
        plt.ylabel('Variants', fontsize=30)
        
        ax = plt.gca()
        
        x_ticks = list(range(len(ld_table.columns)))
        x_ticklabels = ld_table.columns
        x_filtered_ticks = [i for i, label in enumerate(x_ticklabels) if label in list_TEs]
        
        y_ticks = list(range(len(ld_table.index)))
        y_ticklabels = ld_table.index
        y_filtered_ticks = [i for i, label in enumerate(y_ticklabels) if label in list_TEs]
        
        plt.xticks(x_filtered_ticks, 
                [x_ticklabels[i] for i in x_filtered_ticks], 
                rotation=90, ha='right', fontsize=30)
        
        plt.yticks(y_filtered_ticks, 
                [y_ticklabels[i] for i in y_filtered_ticks], 
                rotation=0, fontsize=30)
        
        for ax_labels, filtered_ticks in [
            (plt.gca().get_xticklabels(), x_filtered_ticks),
            (plt.gca().get_yticklabels(), y_filtered_ticks)
        ]:
            for label, tick in zip(ax_labels, filtered_ticks):
                if label.get_text() == highlight_variant:
                    label.set_color('red')
                    label.set_fontweight('bold')
        
        plt.tight_layout()
        
        plt.savefig(output_file)
        plt.close()
        
        print(f"Heatmap for r^2 saved to {output_file}")

    @staticmethod 
    def sign_rsquared_retrieve_genotype(genotype_dict, significant_rsquared):
        info_genotypes = {}
        
        for key, value in genotype_dict.items():
            # key is a tuple of TEs, value contains combination_1 and combination_2
            te_list = list(key)  # Convert tuple to list for indexing
            
            for i, item in enumerate(te_list):
                # Use combination based on position in the tuple
                if i == 0:  # First TE in tuple gets combination_1
                    genotype = value['combination_1']
                else:  # Second TE in tuple gets combination_2
                    genotype = value['combination_2']
                info_genotypes[item] = genotype
        
        df = pd.DataFrame(list(info_genotypes.items()), columns=['TE', 'Genotypes'])
        return df
    
    @staticmethod
    def same_genotype(genotype_df):
        if genotype_df.empty:
            return genotype_df, [] 

        genotype_df['Genotypes_str'] = genotype_df['Genotypes'].astype(str)

        first_row_genotype = genotype_df.iloc[0]['Genotypes_str']

        matching_rows = genotype_df[genotype_df['Genotypes_str'] == first_row_genotype]

        if matching_rows.empty:
            return matching_rows, []

        matching_identifiers = matching_rows.iloc[:, 0].tolist()

        return matching_rows, matching_identifiers
    
    @staticmethod 
    def combining_TEs(complete_df, variation):
        first_row = complete_df.iloc[0:1].reset_index(drop=True)
        complete_df = complete_df.iloc[1:].reset_index(drop=True)
        
        split_columns = complete_df['Variant_name'].str.extract(r'([^_]+_[^_]+)_(\d+)_(\d+)_(.*)')
        
        split_columns.columns = ['Chr', 'position_start', 'position_end', 'TE_name']
        
        split_columns['position_start'] = split_columns['position_start'].astype(int)
        split_columns['position_end'] = split_columns['position_end'].astype(int)
        
        complete_df[['Chr', 'position_start', 'position_end', 'TE_name']] = split_columns
        
        group_cols = ['TE_name']

        final_rows = []
        num_removed_lines = 0 
        number_of_chosen_rows = 0
        combined_groups = []  

        for name, group in complete_df.groupby(group_cols):
            num_rows = len(group)
            min_start_position = float('inf')
            chosen_row = None

            combined_row = group.iloc[0].copy()  # Start with the first row of the group
            genotype_cols = complete_df.columns.difference(['Variant_name', 'Chr', 'position_start', 'position_end', 'TE_name'])
            
            # Iterate over rows in the group
            for i in range(1, num_rows):
                info_i = group.iloc[i]
                position_start_i = info_i['position_start']
                position_end_i = info_i['position_end']

                similar_found = False

                # Check against subsequent rows j in the group
                for j in range(i, num_rows):
                    info_j = group.iloc[j]
                    position_start_j = info_j['position_start']
                    position_end_j = info_j['position_end']
                    
                    # Check if the start positions differ by at least `variation`
                    if (abs(position_start_j - position_start_i) <= variation or
                            abs(position_end_i - position_end_j) <= variation):
                        similar_found = True
                        num_removed_lines += 1

                        # Combine genotype information from row `info_j`
                        for col in genotype_cols:
                            combined_value = combined_row[col]
                            new_value = info_j[col]

                            # Choose the non-zero value, or custom logic to merge
                            if combined_value == (0, 0) and new_value != (0, 0):
                                combined_row[col] = new_value
                            elif combined_value != (0, 0) and new_value != (0, 0):
                                combined_row[col] = max(combined_value, new_value)
                        break
                
                if similar_found:
                    # Log combined rows for clarity
                    combined_groups.append((name, group.iloc[i:]))
                else:
                    number_of_chosen_rows += 1

            # After combining, store the merged row
            final_rows.append(combined_row)

        # Create a new DataFrame from final_rows
        final_df = pd.DataFrame(final_rows)
        
        '''# Print combined group information
        for group_name, combined in combined_groups:
            print(f"TE group '{group_name}' had the following rows combined:")
            for row in combined:
                print(row)
            print()  # Add a newline for readability'''
        
        print(f"Number of chosen lines: {number_of_chosen_rows}")
        print(f'This is the final df: {final_df.shape}')
        
        final_df = final_df.drop(columns=['Chr', 'position_start', 'position_end', 'TE_name'])
        final_df = pd.concat([first_row, final_df], ignore_index=True)
        
        return final_df