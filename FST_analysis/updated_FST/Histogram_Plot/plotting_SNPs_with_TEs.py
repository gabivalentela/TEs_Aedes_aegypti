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
    cols = 4  
    rows = -(-num_combinations // cols)  

    # Create a figure with subplots
    fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 5 * rows))
    axes = axes.flatten()  

    for idx, (countries, value) in enumerate(files_dict.items()):
        df_snps  = value['SNPs']
        df_tes = value['combined_ref_nonref']
        pair_type = 'combined_ref_nonref'  

        # Get FST values for both datasets
        fst_values_tes = df_tes['FST'].values
        fst_values_snps = df_snps.values

        # Compute histogram counts and bin edges for TEs
        counts_tes, bins = np.histogram(fst_values_tes, bins=20, range=(0, 1))
        total_count_tes = sum(counts_tes)
        percentages_tes = (counts_tes / total_count_tes) * 100

        # Compute histogram counts and bin edges for SNPs
        counts_snps, _ = np.histogram(fst_values_snps, bins=bins)  
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

    for idx in range(num_combinations, len(axes)):
        fig.delaxes(axes[idx])

    plot_type = 'combined_ref_nonref'

    # Adjust layout and save the figure
    plt.tight_layout()
    plt.savefig(f'FST_histograms_{plot_type}_1000.png')

def combine_df(files_SNPs_TEs):
    # Collect unique country combinations 
    combined_data = {}  
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

def bin_and_plot_combined_fst(group, files_dict, target_countries, bin_edges=None):
    if bin_edges is None:
        bin_edges = np.arange(0, 1.01, 0.1)  

    # Define a color mapping for each country pair
    color_mapping = {
        ("Colombia", "USA"): "#8B0000",  # Dark Red
        ("Brazil", "Colombia"): "#9D00FF",  # Deep Purple
        ("Colombia", "Kenya"): "#005F73",  # Dark Teal Blue
        ("Colombia", "Senegal"): "#007200",  # Deep Forest Green
        ("Colombia", "Gabon"): "#C2185B",  # Dark Pink/Magenta
        ("Brazil", "USA"): "#FF8C00",  # Dark Orange
        ("Kenya", "USA"): "#4682B4",  # Steel Blue
        ("Senegal", "USA"): "#556B2F",  # Dark Olive Green
        ("Gabon", "USA"): "#DC143C",  # Crimson
        ("Brazil", "Kenya"): "#8A2BE2",  # Blue Violet
        ("Brazil", "Senegal"): "#D2691E",  # Chocolate Brown
        ("Brazil", "Gabon"): "#FF4500",  # Orange Red
        ("Kenya", "Senegal"): "#2E8B57",  # Sea Green
        ("Gabon", "Kenya"): "#8B4513",  # Saddle Brown
        ("Gabon", "Senegal"): "#B22222"  # Firebrick
    }

    # Initialize storage for bin proportions
    bin_proportions = {}

    for countries, value in files_dict.items():
        if countries not in target_countries:
            continue 

        df_snps = value['SNPs']
        df_tes = value['combined_ref_nonref']

        # Extract FST values
        fst_values_tes = df_tes['FST'].values
        fst_values_snps = df_snps.values.flatten()

        # Bin counts
        snp_bin_counts = pd.Series(pd.cut(fst_values_snps, bins=bin_edges, labels=False, include_lowest=True)) \
                         .value_counts(sort=False).sort_index() \
                         .reindex(range(len(bin_edges) - 1), fill_value=0)

        te_bin_counts = pd.Series(pd.cut(fst_values_tes, bins=bin_edges, labels=False, include_lowest=True)) \
                        .value_counts(sort=False).sort_index() \
                        .reindex(range(len(bin_edges) - 1), fill_value=0)

        # Convert to proportions
        total_snps = snp_bin_counts.sum()
        total_tes = te_bin_counts.sum()

        snp_proportions = snp_bin_counts / total_snps if total_snps > 0 else snp_bin_counts
        te_proportions = te_bin_counts / total_tes if total_tes > 0 else te_bin_counts

        # Store proportions
        bin_proportions[countries] = {
            'SNPs': snp_proportions,
            'TEs': te_proportions
        }

    # Create a combined bar plot
    plt.figure(figsize=(12, 8))
    bar_width = 0.1
    bin_labels = [f"{bin_edges[i]:.1f}-{bin_edges[i+1]:.1f}" for i in range(len(bin_edges) - 1)]
    x = range(len(bin_labels))
    offset = 0  

    for countries, proportions in bin_proportions.items():
        color = color_mapping.get(tuple(sorted(countries)), "#808080") 

        plt.bar(
            [xi + offset for xi in x],
            proportions['SNPs'],
            width=bar_width,
            label=f"SNPs ({countries[0]}-{countries[1]})",
            color=color
        )

        plt.bar(
            [xi + offset + bar_width for xi in x],
            proportions['TEs'],
            width=bar_width,
            label=f"TEs ({countries[0]}-{countries[1]})",
            color=color,
            hatch="//"
        )

        offset += bar_width * 2

    # Add labels and legend
    plt.xlabel('FST Bins', fontsize=25)
    plt.ylabel('Log₁₀(Proportion)', fontsize=25)
    plt.title('Proportion of FST Values for TEs and SNPs', fontsize=15)
    plt.xticks([xi + bar_width * len(bin_proportions) / 2 for xi in x], bin_labels, rotation=45, fontsize=20)
    plt.yticks(fontsize=14)
    plt.legend(fontsize=16, loc='upper left', bbox_to_anchor=(1, 1))  
    plt.yscale('log')  
    plt.ylim(1e-4, 1)  
    plt.tight_layout()

    # Save or show the plot
    plt.savefig(f"new_log_proportion_fst_bins_plot_{group}.svg")

TEs_data = '/work/users/g/a/gabivla/lab/SV_mosquitos/FST_analysis/TEs_absence_evidence'
SNPs_data = '/work/users/g/a/gabivla/lab/SV_mosquitos/FST_analysis/SNPs'

# Load data
list_files_tes = os.listdir(TEs_data)
list_files_snps = os.listdir(SNPs_data)
files_SNPs_TEs = pair_data(list_files_tes, list_files_snps, TEs_data, SNPs_data)

combined_df = combine_df(files_SNPs_TEs)
#initial plot of the data
plot_fst_histograms(combined_df)
# Target countries plot
target_countries_1 = [
    ('Brazil', 'Colombia'),
    ('Colombia', 'USA'),
    ('Brazil', 'USA')
]

target_countries_2 = [
    ('Gabon', 'Senegal'),
    ('Gabon', 'Kenya'),
    ('Kenya', 'Senegal')
]

target_countries_3 = [
    ('Brazil', 'Kenya'),
    ('Colombia', 'Kenya'),
    ('Kenya', 'USA')]

target_countries_4 = [
    ('Colombia', 'Kenya'),
    ('Colombia', 'Senegal'),
    ('Colombia', 'Gabon')]

target_countries_5 = [
    ('Brazil', 'Gabon'),
    ('Gabon', 'USA'),
    ('Colombia', 'Gabon')]

target_countries_6 = [
    ('Brazil', 'Senegal'),
    ('Colombia', 'Senegal'),
    ('Senegal', 'USA')]

# Run the function
binned_data1 = bin_and_plot_combined_fst('americas', combined_df, target_countries_1)
binned_data2 = bin_and_plot_combined_fst('african', combined_df, target_countries_2)
binned_data2 = bin_and_plot_combined_fst('americas_africa', combined_df, target_countries_3)
binned_data2 = bin_and_plot_combined_fst('africa_america', combined_df, target_countries_4)
binned_data2 = bin_and_plot_combined_fst('africa_america2', combined_df, target_countries_5)
binned_data2 = bin_and_plot_combined_fst('africa_america3', combined_df, target_countries_6)
