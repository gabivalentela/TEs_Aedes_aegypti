import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def pair_data(list_files_tes, list_files_snps, TEs_data, SNPs_data):
    dict_data = {}
    
    for file in list_files_tes:
        # Construct the full path for the TE file
        path_TE = os.path.join(TEs_data, file)
        
        # Extract country1, country2, and pair type using regex
        country_1 = re.search(r'FST_values_(.+)_(.+)_.+.csv', file).group(1)
        country_2 = re.search(r'FST_values_(.+)_(.+)_.+.csv', file).group(2)
        pair_type = re.search(r'FST_values_(.+)_(.+)_(.+).csv', file).group(3)
        
        # Sort the countries alphabetically to ensure consistent key order
        sorted_countries = tuple(sorted([country_1, country_2]))
        print(f"Processing pair: {sorted_countries} with pair type: {pair_type}")
        
        # Find the corresponding SNP file
        for SNPs in list_files_snps:
            if country_1 in SNPs and country_2 in SNPs:
                # Construct the full path for the SNP file
                snp_path = os.path.join(SNPs_data, SNPs)
                
                # Initialize an entry for this country combination if not already present
                if sorted_countries not in dict_data:
                    dict_data[sorted_countries] = {}
                
                # Add the pair_type data to the dictionary
                dict_data[sorted_countries][pair_type] = [path_TE, snp_path]
    
    return dict_data


# Define a function for plotting histograms
def plot_fst_histograms(files_dict):
    num_combinations = 15
    cols = 4  # Number of columns in the subplot grid
    rows = -(-num_combinations // cols)  # Ceiling division to determine rows

    # Create a figure with subplots
    fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 5 * rows))
    axes = axes.flatten()  # Flatten axes for easier indexing

    for idx, (countries, value) in enumerate(files_dict.items()):
        df_snps  = value['SNPs']
        df_tes = value['combined_ref_nonref']
        pair_type = 'combined_ref_nonref'  # Use type for labeling

        # Get FST values for both datasets
        fst_values_tes = df_tes['FST'].values
        fst_values_snps = df_snps.values

        # Compute histogram counts and bin edges for TEs
        counts_tes, bins = np.histogram(fst_values_tes, bins=20, range=(0, 1))
        total_count_tes = sum(counts_tes)
        percentages_tes = (counts_tes / total_count_tes) * 100

        # Compute histogram counts and bin edges for SNPs
        counts_snps, _ = np.histogram(fst_values_snps, bins=bins)  # Use the same bins for consistency
        total_count_snps = sum(counts_snps)
        percentages_snps = (counts_snps / total_count_snps) * 100

        # Compute bin centers for plotting lines
        bin_centers = (bins[:-1] + bins[1:]) / 2

        # Plot on the current axis
        ax = axes[idx]

        # Plot TEs as a blue line
        ax.plot(bin_centers, percentages_tes, color='blue', alpha=0.7, label='TEs', linewidth=2)

        # Plot SNPs as an orange line
        ax.plot(bin_centers, percentages_snps, color='orange', alpha=0.7, label='SNPs', linewidth=2, linestyle='--')

        # Add labels and title
        ax.set_title(f'{countries} ({pair_type})')
        ax.set_xlabel('FST Value')
        ax.set_ylabel('Percentage')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 100)
        ax.legend()

    # Turn off unused subplots
    for idx in range(num_combinations, len(axes)):
        fig.delaxes(axes[idx])

    plot_type = 'combined_ref_nonref'

    # Adjust layout and save the figure
    plt.tight_layout()
    plt.savefig(f'FST_histograms_{plot_type}_1000.png')

def combine_df(files_SNPs_TEs):
    # Collect unique country combinations (excluding the pair type)
    combined_data = {}  # Dictionary to store combined DataFrames for each country
    for country, files in files_SNPs_TEs.items():
        print(f"Processing country: {country}")
        
        nonref = files['nonref'][0]
        ref = files['ref'][0]
        SNPs = files['ref'][1]
        
        print(f"Non-ref file: {nonref}")
        print(f"Ref file: {ref}")
        print(f"SNPs file: {SNPs}")

        # Read in the CSV files
        df_ref = pd.read_csv(ref)
        df_nonref = pd.read_csv(nonref)
        df_SNPs = pd.read_csv(SNPs)

        # Combine ref and non-ref DataFrames
        df_combined = pd.concat([df_ref, df_nonref], ignore_index=True)

        # Store the combined DataFrame in the dictionary
        combined_data[country] = {
            "combined_ref_nonref": df_combined,
            "SNPs": df_SNPs
        }

    return combined_data

import pandas as pd
import matplotlib.pyplot as plt

def bin_and_plot_combined_fst(files_dict, bin_edges=None):
    if bin_edges is None:
        bin_edges = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0]  # Default bin edges

    for countries, value in files_dict.items():
        # Extract data for this country pair
        df_snps = value['SNPs']
        df_tes = value['combined_ref_nonref']

        # Extract FST values
        fst_values_tes = df_tes['FST'].values
        fst_values_snps = df_snps.values.flatten()

        # Binning FST values
        snp_bins = pd.cut(fst_values_snps, bins=bin_edges, labels=False, include_lowest=True)
        te_bins = pd.cut(fst_values_tes, bins=bin_edges, labels=False, include_lowest=True)

        # Count values in bins
        snp_bin_counts = pd.Series(snp_bins).value_counts(sort=False).sort_index()
        te_bin_counts = pd.Series(te_bins).value_counts(sort=False).sort_index()

        # Create a bar plot for this country pair
        plt.figure(figsize=(12, 8))
        plt.yscale('log')  # Set y-axis to logarithmic scale
        bar_width = 0.4
        bin_labels = [f"{bin_edges[i]}-{bin_edges[i+1]}" for i in range(len(bin_edges)-1)]
        x = range(len(bin_labels))

        # Plot SNP counts
        plt.bar(
            [xi - bar_width / 2 for xi in x],
            snp_bin_counts,
            width=bar_width,
            label="SNPs",
            color="darkblue"
        )

        # Plot TE counts
        plt.bar(
            [xi + bar_width / 2 for xi in x],
            te_bin_counts,
            width=bar_width,
            label="TEs",
            color="blue",
            hatch="//"
        )

        # Add labels and legend
        plt.xlabel('FST Bins')
        plt.ylabel('Count')
        plt.title(f'FST Values for TEs and SNPs for {countries[0]}-{countries[1]}')
        plt.xticks(x, bin_labels, rotation=45)
        plt.legend(loc='upper right')
        plt.tight_layout()

        # Save the plot to a file
        plot_filename = f"fst_bins_plot_1000_{countries[0]}_{countries[1]}.png"
        plt.savefig(plot_filename)
        plt.close()  # Close the plot to free up memory


TEs_data = '/work/users/g/a/gabivla/lab/SV_mosquitos/FST_data/TEs_variation_1000'
SNPs_data = '/work/users/g/a/gabivla/lab/SV_mosquitos/FST_data/SNPs'

# Load data
list_files_tes = os.listdir(TEs_data)
list_files_snps = os.listdir(SNPs_data)
files_SNPs_TEs = pair_data(list_files_tes, list_files_snps, TEs_data, SNPs_data)

combined_df = combine_df(files_SNPs_TEs)
#initial plot of the data
plot_fst_histograms(combined_df)

# Run the function
binned_data1 = bin_and_plot_combined_fst(combined_df)
