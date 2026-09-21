import pandas as pd
import os
import re
import sys

# Load TE reference data
ref_df = sys.argv[1] 
df = pd.read_csv(ref_df)

# Your existing location mapping code
location_path = '/nas/longleaf/home/gabivla/lab/SV_mosquitos/location_information/mcclintock_output_locations.csv'
location_df = pd.read_csv(location_path, usecols=['Genome_code', 'country'])
location_dict = dict(zip(location_df['Genome_code'], location_df['country']))

base_cols = ['Unnamed: 0', 'chromosome', 'position_start', 'position_end', 'TE_name', 'strand']
genome_cols = [col for col in df.columns if col not in base_cols]
rename_dict = {col: f"{col}_{location_dict.get(col, '')}" for col in genome_cols}
df = df.rename(columns=rename_dict)

# Function to find TEs that contain sweep locations
def find_TEs_containing_sweeps(te_df, sweep_df, country, chrom):
    """
    Find TEs where the sweep location falls within the TE boundaries
    
    Parameters:
    - te_df: DataFrame with TE information
    - sweep_df: DataFrame with sweep finder results
    - country: country name
    - chrom: chromosome name
    
    Returns:
    - DataFrame with TEs that contain sweep locations and their sweep statistics
    """
    overlapping_results = []
    
    for _, sweep_row in sweep_df.iterrows():
        sweep_chrom = sweep_row['chrom'] 
        location = sweep_row['location']
        lr = sweep_row['LR']
        alpha = sweep_row['alpha']
        pop = sweep_row['pop']

        # Find TEs where the sweep location falls WITHIN the TE boundaries
        matching_tes = te_df[
            (te_df['chromosome'] == sweep_chrom) &
            (te_df['position_start'] <= location) &
            (te_df['position_end'] >= location)
        ].copy()
        
        if not matching_tes.empty:
            matching_tes['sweep_location'] = location
            matching_tes['sweep_LR'] = lr
            matching_tes['sweep_alpha'] = alpha
            matching_tes['country'] = country
            matching_tes['chromosome_file'] = chrom
            matching_tes['population'] = pop
            # Distance is 0 since sweep is within TE, but we can calculate center distance
            te_center = (matching_tes['position_start'] + matching_tes['position_end']) / 2
            matching_tes['distance_from_TE_center'] = abs(te_center - location)
            overlapping_results.append(matching_tes)
    
    if overlapping_results:
        return pd.concat(overlapping_results, ignore_index=True)
    else:
        return pd.DataFrame()

# Get sweep file from command line argument
sweep_file = sys.argv[2]
ref_or_nonref = sys.argv[3]  # 'ref' or 'nonref'

# Extract country and chromosome from filename
country = re.search(r'.+no_relatives_(.+)_raw', sweep_file).group(1)
chrom = re.search(r'(AaegL5_\d+)_no_relatives.+', sweep_file).group(1)

print(f"Processing {country}, Chromosome: {chrom}")
print(f"Sweep file: {sweep_file}")

# Read sweep results
sweep_df = pd.read_csv(sweep_file, sep='\t')
print(f"Loaded {len(sweep_df)} sweep regions")

# Find TEs that contain sweep locations
overlapping = find_TEs_containing_sweeps(df, sweep_df, country, chrom)

# Save results
if not overlapping.empty:
    # Sort by sweep LR (likelihood ratio) to prioritize strongest sweeps
    overlapping = overlapping.sort_values('sweep_LR', ascending=False)
    
    # Create output filename based on input
    output_path = f'./results/TEs_containing_sweeps_{country}_{chrom}_{ref_or_nonref}.csv'
    overlapping.to_csv(output_path, index=False)
    
    print(f"\n{'='*50}")
    print(f"Found {len(overlapping)} TEs containing sweep locations")
    print(f"Results saved to: {output_path}")
    
    # Show summary statistics
    print(f"\nSweep LR statistics:")
    print(f"  Mean: {overlapping['sweep_LR'].mean():.2f}")
    print(f"  Max: {overlapping['sweep_LR'].max():.2f}")
    print(f"  Min: {overlapping['sweep_LR'].min():.2f}")
else:
    print("\nNo TEs contain sweep locations for this file")