import pandas as pd

def check_for_evidence_at_least_one_genome(df, genome_cols):
    """
    Check if at least one genome has '0/0' genotype for each TE insertion.
    
    Parameters:
    df - DataFrame with TE data
    genome_cols - List of columns representing different genomes
    
    Returns:
    DataFrame with a new column indicating if at least one genome has '0/0'
    """
    # Create a new column to indicate if at least one genome has '0/0'
    df['evidence_at_least_one_genome'] = df[genome_cols].apply(lambda row: '0/0' in row.values, axis=1)

    return df

def update_genotypes(df, genome_cols):
    """
    Update genotype information in the DataFrame.
    If evidence for at least one genome 0/0 - change all ./. to 0/0.
    
    Parameters:
    df - DataFrame with TE data and evidence column
    genome_cols - List of columns representing different genomes
    
    Returns:
    DataFrame with updated genotypes
    """
    # Make a copy to avoid modifying the original
    df_updated = df.copy()
    
    # Get rows where evidence_at_least_one_genome is True
    evidence_mask = df_updated['evidence_at_least_one_genome'] == True
    
    # For rows with evidence, replace './.' with '0/0' in genome columns
    for col in genome_cols:
        df_updated.loc[evidence_mask, col] = df_updated.loc[evidence_mask, col].replace('./.', '0/0')
    
    return df_updated

# Running functions
reference_path = '/work/users/g/a/gabivla/lab/SV_mosquitos/VCF_genotypes_TEs/separate_chromosomes_windows_reference/merged_reference/VCF_reference_TEs_post_processing_600.csv'

# open file as df 
df = pd.read_csv(reference_path)

# fixing df
df = df.drop(columns=['Unnamed: 0'])

# Define base columns that should not be modified
base_cols = ['chromosome', 'position_start', 'position_end', 'TE_name', 'strand']

# Get genome columns at once
genome_cols = [col for col in df.columns if col not in base_cols]

# Process the data
df_evidence = check_for_evidence_at_least_one_genome(df, genome_cols)
df_updated = update_genotypes(df_evidence, genome_cols) 

print("Sample of updated data:")
print(df_updated.head())

# Show summary of what was updated vs what was kept
evidence_true_rows = df_updated[df_updated['evidence_at_least_one_genome'] == True]
evidence_false_rows = df_updated[df_updated['evidence_at_least_one_genome'] == False]

print(f"\nSummary:")
print(f"Rows with evidence=True (./. changed to 0/0): {len(evidence_true_rows)}")
print(f"Rows with evidence=False (./. kept as ./.): {len(evidence_false_rows)}")
print(f"Total rows: {len(df_updated)}")

# remove the evidence column after updating genotypes
df_updated = df_updated.drop(columns=['evidence_at_least_one_genome'])
# make df with updated genotypes
df_updated.to_csv('/work/users/g/a/gabivla/lab/SV_mosquitos/VCF_genotypes_TEs/separate_chromosomes_windows_reference/merged_reference/VCF_reference_TEs_post_processing_600_absence_evidence.csv', index=False)