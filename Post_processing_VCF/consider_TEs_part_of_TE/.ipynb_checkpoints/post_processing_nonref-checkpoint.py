import pandas as pd
import polars as pl
import sys

def remove_rows_without_TEs(df):
    # Identify rows that contain at least one '1/1' in columns after the first 5
    mask = (df.iloc[:, 5:] == '1/1').any(axis=1)
    # Keep only the rows that satisfy the condition
    df = df[mask]
    return df

def check_TEs_parts(df, max_distance=1000):
    merged_rows = []
    seen_indices = set()  # Track indices that were merged
    
    # Process all TEs in a consistent order
    for te_name in sorted(df.keys()):
        te_df = df[te_name]
        pl_df = pl.DataFrame(te_df)  # Ensure it's a Polars DataFrame
        
        # First collect all chromosomes and sort them to ensure consistent ordering
        all_chromosomes = sorted(pl_df["chromosome"].unique())
        # Process each chromosome in the sorted order
        for chromosome in all_chromosomes:
            # Filter DataFrame for this chromosome and sort by position
            group = pl_df.filter(pl.col("chromosome") == chromosome).sort("position_start")
            num_rows = group.shape[0]  # Get number of rows
            
            # Debug: Print all positions in this chromosome
            print(f"\nPositions in chromosome {chromosome}:")
            for i in range(num_rows):
                row = group.row(i, named=True)
                print(f"{row['position_start']} {row['position_end']} {row['TE_name']} {row.get('strand', '-')} {row['chromosome']}")
            
            # Check for overlaps or nearby TEs
            print("\nOverlaps or nearby TEs detected:")
            for i in range(num_rows):
                if i in seen_indices:
                    continue  # Skip if already merged
                
                merged = False  # Track if this row gets merged
                info_i = group.row(i, named=True)
                position_start_i = int(info_i["position_start"])
                position_end_i = int(info_i["position_end"])
                te_name_i = info_i["TE_name"]
                strand_i = info_i.get("strand", "-")
                chromosome_i = info_i["chromosome"]

                # Check against subsequent rows
                for j in range(i + 1, num_rows):
                    if j in seen_indices:
                        continue  # Skip if already merged

                    info_j = group.row(j, named=True)
                    position_start_j = int(info_j["position_start"])
                    position_end_j = int(info_j["position_end"])
                    te_name_j = info_j["TE_name"]
                    strand_j = info_j.get("strand", "-")
                    chromosome_j = info_j["chromosome"]
 
                    # Only check positions if the TE names, strands, and chromosomes match
                    if te_name_i == te_name_j and strand_i == strand_j and chromosome_i == chromosome_j:
                        # Check if regions overlap
                        regions_overlap = (
                            (position_start_i <= position_start_j <= position_end_i) or 
                            (position_start_i <= position_end_j <= position_end_i) or
                            (position_start_j <= position_start_i <= position_end_j)
                        )
                        
                        # Calculate distance between regions
                        if position_start_j > position_end_i:
                            distance = position_start_j - position_end_i
                        elif position_start_i > position_end_j:
                            distance = position_start_i - position_end_j
                        else:
                            distance = 0  # They overlap
                        
                        # Check if regions are within max_distance of each other
                        regions_nearby = distance <= max_distance
                        
                        if regions_overlap or regions_nearby:
                            # Merge regions
                            merged_start = min(position_start_i, position_start_j)
                            merged_end = max(position_end_i, position_end_j)
                            
                            if regions_overlap:
                                reason = "overlap"
                            else:
                                reason = f"nearby (distance={distance}bp)"
                            
                            print(f"Merge needed ({reason}): {position_start_i}-{position_end_i} and {position_start_j}-{position_end_j}")
                            print(f"Merged region: {merged_start}-{merged_end}")
                            
                            # Merge genotypes according to the rules
                            genotypes = {col: info_i[col] for col in pl_df.columns if col not in ["chromosome", "position_start", "position_end", "TE_name", "strand"]}

                            for col in genotypes.keys():
                                if (info_j[col] == "1/1") and (info_i[col] == "1/1"):
                                    genotypes[col] = "1/1"
                                elif (info_j[col] == "0/0") and (info_i[col] == "1/1") or (info_j[col] == "1/1") and (info_i[col] == "0/0"):
                                    genotypes[col] = "1/1"
                                elif (info_j[col] == "./.") and (info_i[col] == "1/1") or (info_j[col] == "1/1") and (info_i[col] == "./."):
                                    genotypes[col] = "1/1"
                                elif (info_j[col] == "0/0") and (info_i[col] == "0/0"):
                                    genotypes[col] = "0/0"
                                elif (info_j[col] == "./.") and (info_i[col] == "./."):
                                    genotypes[col] = "./."

                            merged_rows.append({
                                "chromosome": chromosome_i,
                                "strand": strand_i,
                                "TE_name": te_name_i,
                                "position_start": merged_start,
                                "position_end": merged_end,
                                **genotypes
                            })

                            seen_indices.add(i)
                            seen_indices.add(j)
                            merged = True
                            break  # Stop checking once merged

                # If the row was not merged, add it as it is
                if not merged and i not in seen_indices:
                    merged_rows.append(info_i)

    # Create a new DataFrame with merged results and sort by chromosome and position_start
    if merged_rows:
        merged_df = pl.DataFrame(merged_rows).sort(["chromosome", "position_start"])
        print(merged_df)
        return merged_df
    else:
        return None
   
def one_df_for_each_TE(df):
    # Get all unique TE names
    TE_names = df['TE_name'].unique()
    # Create a dictionary with each TE name as a key and the corresponding DataFrame as a value
    te_dfs = {te: df[df['TE_name'] == te] for te in TE_names}
    return te_dfs

def checking_for_almost_fixed_fixed_TEs(df):
    df = df.to_pandas()
    # Get genotype columns (ignoring the first 5 columns)
    genotype_data = df.iloc[:, 5:]

    # Condition: At least two columns must contain "0/0"
    at_least_two_ref = (genotype_data == "0/0").sum(axis=1) >= 2

    # Keep only rows that satisfy this condition
    df = df[at_least_two_ref]

    return df

non_reference_VCF = '/work/users/g/a/gabivla/lab/SV_mosquitos/VCF_genotypes_TEs/separate_chromosomes_windows_reference/merged_nonreference/VCF_nonreference_TEs_pre_processing.csv'
df_nonref = pd.read_csv(non_reference_VCF)
df_nonref = df_nonref.sort_values(by=['chromosome', 'position_start'])
df_nonref_updated = remove_rows_without_TEs(df_nonref)
# check for duplicated rows in this df
duplicated_rows = df_nonref_updated[df_nonref_updated.duplicated()]
print(duplicated_rows)
df_nonref_updated.to_csv('/work/users/g/a/gabivla/lab/SV_mosquitos/VCF_genotypes_TEs/separate_chromosomes_windows_reference/merged_nonreference/VCF_nonreference_post_processing.csv')
