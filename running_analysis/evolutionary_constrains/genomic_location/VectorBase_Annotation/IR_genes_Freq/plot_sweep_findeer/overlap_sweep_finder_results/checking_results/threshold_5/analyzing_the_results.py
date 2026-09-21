import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df_ref = pd.read_csv('../../results/reference_grouped.csv')
df_nonref = pd.read_csv('../../results/nonreference_grouped.csv')

def split_sweep_lr_by_population(df):
    """
    Split the sweep_LR column into separate columns for each population.
    Handles rows with different numbers of populations.
    """
    # Create a list to store the expanded rows
    expanded_data = []
    
    for idx, row in df.iterrows():
        # Split populations and sweep_LR values
        populations = row['population'].split(';')
        sweep_lr_values = row['sweep_LR'].split(';')
        
        # Create a dictionary for this row with the new columns
        row_dict = row.to_dict()
        
        # Add sweep_LR columns for each population
        for pop, lr_value in zip(populations, sweep_lr_values):
            col_name = f'sweep_LR_{pop}'
            try:
                row_dict[col_name] = float(lr_value)
            except ValueError:
                row_dict[col_name] = np.nan
        
        expanded_data.append(row_dict)
    
    # Convert to dataframe
    df_expanded = pd.DataFrame(expanded_data)
    
    return df_expanded


def get_sweep_lr_distributions(df):
    """
    Get the distribution of sweep_LR values for each population/country.
    Returns a dictionary with population names as keys and arrays of sweep_LR values.
    """
    # Find all sweep_LR columns
    sweep_lr_cols = [col for col in df.columns if col.startswith('sweep_LR_')]
    
    distributions = {}
    
    for col in sweep_lr_cols:
        # Extract population name
        pop_name = col.replace('sweep_LR_', '')
        
        # Get non-null values
        values = df[col].dropna()
        
        distributions[pop_name] = {
            'values': values.values,
            'count': len(values),
            'mean': values.mean() if len(values) > 0 else np.nan,
            'median': values.median() if len(values) > 0 else np.nan,
            'std': values.std() if len(values) > 0 else np.nan,
            'min': values.min() if len(values) > 0 else np.nan,
            'max': values.max() if len(values) > 0 else np.nan,
            'non_zero_count': (values > 0).sum()
        }
    
    return distributions


def plot_sweep_lr_distributions(distributions, title="Sweep LR Distributions by Population"):
    """
    Plot histograms of sweep_LR distributions for each population.
    """
    n_populations = len(distributions)
    
    # Calculate grid dimensions
    n_cols = 3
    n_rows = (n_populations + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 5*n_rows))
    axes = axes.flatten() if n_populations > 1 else [axes]
    
    for idx, (pop_name, dist_info) in enumerate(sorted(distributions.items())):
        ax = axes[idx]
        values = dist_info['values']
        
        # Filter out zeros for better visualization (optional)
        non_zero_values = values[values > 0]
        
        if len(non_zero_values) > 0:
            ax.hist(non_zero_values, bins=50, edgecolor='black', alpha=0.7)
            ax.set_xlabel('Sweep LR')
            ax.set_ylabel('Frequency')
            ax.set_title(f'{pop_name}\n(n={dist_info["non_zero_count"]} non-zero, mean={dist_info["mean"]:.2f})')
            ax.grid(True, alpha=0.3)
        else:
            ax.text(0.5, 0.5, 'No non-zero values', ha='center', va='center')
            ax.set_title(f'{pop_name}\n(n=0)')
    
    # Hide unused subplots
    for idx in range(n_populations, len(axes)):
        axes[idx].axis('off')
    
    plt.suptitle(title, fontsize=16, y=1.00)
    plt.tight_layout()
    return fig


# Apply the function to both dataframes
print("Processing reference data...")
df_ref_expanded = split_sweep_lr_by_population(df_ref)
print("Processing non-reference data...")
df_nonref_expanded = split_sweep_lr_by_population(df_nonref)

# Get distributions
print("\nGetting distributions for reference data...")
ref_distributions = get_sweep_lr_distributions(df_ref_expanded)
print("\nGetting distributions for non-reference data...")
nonref_distributions = get_sweep_lr_distributions(df_nonref_expanded)

# Print summary statistics
print("\n" + "="*60)
print("REFERENCE DATA - Summary Statistics:")
print("="*60)
for pop_name, dist_info in sorted(ref_distributions.items()):
    print(f"\n{pop_name}:")
    print(f"  Total values: {dist_info['count']}")
    print(f"  Non-zero values: {dist_info['non_zero_count']}")
    print(f"  Mean: {dist_info['mean']:.2f}")
    print(f"  Median: {dist_info['median']:.2f}")
    print(f"  Std: {dist_info['std']:.2f}")
    print(f"  Min: {dist_info['min']:.2f}")
    print(f"  Max: {dist_info['max']:.2f}")

print("\n" + "="*60)
print("NON-REFERENCE DATA - Summary Statistics:")
print("="*60)
for pop_name, dist_info in sorted(nonref_distributions.items()):
    print(f"\n{pop_name}:")
    print(f"  Total values: {dist_info['count']}")
    print(f"  Non-zero values: {dist_info['non_zero_count']}")
    print(f"  Mean: {dist_info['mean']:.2f}")
    print(f"  Median: {dist_info['median']:.2f}")
    print(f"  Std: {dist_info['std']:.2f}")
    print(f"  Min: {dist_info['min']:.2f}")
    print(f"  Max: {dist_info['max']:.2f}")

# Plot distributions
fig1 = plot_sweep_lr_distributions(ref_distributions, "Reference Data - Sweep LR Distributions")
plt.savefig('ref_sweep_lr_distributions.png', dpi=300, bbox_inches='tight')

fig2 = plot_sweep_lr_distributions(nonref_distributions, "Non-Reference Data - Sweep LR Distributions")
plt.savefig('nonref_sweep_lr_distributions.png', dpi=300, bbox_inches='tight')

plt.show()

print("\nPlots saved as 'ref_sweep_lr_distributions.png' and 'nonref_sweep_lr_distributions.png'")

## start identifying outliers

def get_top_outliers(df, percentile=80):
    """
    Get all rows (TEs) that have at least one sweep_LR value 
    above the specified percentile threshold.
    
    Parameters:
    -----------
    df : DataFrame
        The expanded dataframe with sweep_LR columns
    percentile : float
        The percentile threshold (default 99 for top 1%)
    
    Returns:
    --------
    outliers_df : DataFrame
        Rows that have at least one sweep_LR value above threshold
    thresholds : dict
        Dictionary of thresholds for each population
    """
    # Find all sweep_LR columns
    sweep_lr_cols = [col for col in df.columns if col.startswith('sweep_LR_')]
    
    # Calculate thresholds for each population
    thresholds = {}
    for col in sweep_lr_cols:
        pop_name = col.replace('sweep_LR_', '')
        values = df[col].dropna()
        threshold = np.percentile(values, percentile)
        # just use 200 as the threshold for all countries
        threshold = 200
        thresholds[pop_name] = threshold
    
    # Create a boolean mask for rows that exceed threshold in ANY population
    mask = pd.Series([False] * len(df), index=df.index)
    
    for col in sweep_lr_cols:
        pop_name = col.replace('sweep_LR_', '')
        threshold = thresholds[pop_name]
        # Mark rows where this population's sweep_LR exceeds its threshold
        mask |= (df[col] > threshold)
    
    # Get the outlier rows
    outliers_df = df[mask].copy()
    
    # Add a column showing which populations triggered the outlier status
    def get_outlier_populations(row):
        outlier_pops = []
        for col in sweep_lr_cols:
            pop_name = col.replace('sweep_LR_', '')
            if pd.notna(row[col]) and row[col] > thresholds[pop_name]:
                outlier_pops.append(f"{pop_name}({row[col]:.2f})")
        return '; '.join(outlier_pops)
    
    outliers_df['outlier_populations'] = outliers_df.apply(get_outlier_populations, axis=1)
    
    return outliers_df, thresholds


def summarize_outliers(outliers_df, thresholds, data_name="Data"):
    """
    Print summary statistics about the outliers.
    """
    print("\n" + "="*60)
    print(f"{data_name} - TOP 1% OUTLIERS SUMMARY")
    print("="*60)
    
    print(f"\nTotal number of outlier TEs: {len(outliers_df)}")
    
    print("\nThresholds (99th percentile) for each population:")
    for pop_name, threshold in sorted(thresholds.items()):
        print(f"  {pop_name}: {threshold:.2f}")
    
    # Count how many times each population appears as an outlier
    print("\nNumber of TEs that are outliers for each population:")
    sweep_lr_cols = [col for col in outliers_df.columns if col.startswith('sweep_LR_')]
    
    for col in sorted(sweep_lr_cols):
        pop_name = col.replace('sweep_LR_', '')
        threshold = thresholds[pop_name]
        count = (outliers_df[col] > threshold).sum()
        print(f"  {pop_name}: {count}")
    
    # Show distribution of how many populations each TE is an outlier in
    def count_outlier_pops(row):
        count = 0
        for col in sweep_lr_cols:
            pop_name = col.replace('sweep_LR_', '')
            if pd.notna(row[col]) and row[col] > thresholds[pop_name]:
                count += 1
        return count
    
    outlier_counts = outliers_df.apply(count_outlier_pops, axis=1)
    
    print("\nDistribution of TEs by number of populations they're outliers in:")
    for i in sorted(outlier_counts.unique()):
        count = (outlier_counts == i).sum()
        print(f"  Outlier in {i} population(s): {count} TEs")
    
    return outlier_counts


# Get top 1% outliers for both datasets
print("\nIdentifying top 1% outliers in reference data...")
ref_outliers, ref_thresholds = get_top_outliers(df_ref_expanded, percentile=99)

print("\nIdentifying top 1% outliers in non-reference data...")
nonref_outliers, nonref_thresholds = get_top_outliers(df_nonref_expanded, percentile=99)

# Summarize outliers
ref_outlier_counts = summarize_outliers(ref_outliers, ref_thresholds, "REFERENCE DATA")
nonref_outlier_counts = summarize_outliers(nonref_outliers, nonref_thresholds, "NON-REFERENCE DATA")

# Save outliers to CSV files
ref_outliers.to_csv('reference_top1percent_outliers.csv', index=False)
nonref_outliers.to_csv('nonreference_top1percent_outliers.csv', index=False)

print("\n" + "="*60)
print("Outlier files saved:")
print("  - reference_top1percent_outliers.csv")
print("  - nonreference_top1percent_outliers.csv")
print("="*60)

# Display first few outliers
print("\nFirst 10 reference outliers:")
display_cols = ['population', 'sweep_LR', 'outlier_populations']
# Add sweep_LR individual columns if you want to see them
sweep_lr_cols = [col for col in ref_outliers.columns if col.startswith('sweep_LR_')]
display_cols_full = display_cols + sweep_lr_cols[:3]  # Show first 3 sweep_LR columns as example

print(ref_outliers[display_cols].head(10))

print("\nFirst 10 non-reference outliers:")
print(nonref_outliers[display_cols].head(10))