import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter

# Load datasets
def load_data(nonref_path, ref_path):
    return pd.read_csv(nonref_path), pd.read_csv(ref_path)

# Function to process the dataframe
def process_dataframe(df):
    freq_columns = [col for col in df.columns if col.startswith("count_1/1_")]
    
    # Identify countries with nonzero frequency
    df["nonzero_countries"] = df[freq_columns].apply(
        lambda row: [col.replace("count_1/1_", "") for col in freq_columns if row[col] > 0], axis=1
    )
    
    # Check if a gene_symbol is associated with each row
    df["has_gene_symbol"] = df["gene_symbol"].notna()
    
    return df

# Function to count country occurrences in rows with a gene symbol
def count_countries(df):
    country_counter = Counter()
    
    # Filter rows with a gene symbol
    filtered_df = df[df["has_gene_symbol"]]
    
    # Count country occurrences
    for countries in filtered_df["nonzero_countries"]:
        country_counter.update(countries)
    
    return country_counter

# Function to create bar plot
def plot_country_occurrences(nonref_counts, ref_counts, output_file):
    countries = sorted(set(nonref_counts.keys()).union(ref_counts.keys()))
    nonref_values = [nonref_counts[country] for country in countries]
    ref_values = [ref_counts[country] for country in countries]
    
    x = range(len(countries)) 
    
    plt.figure(figsize=(12, 6))
    bar_width = 0.4
    
    # Plot bars
    plt.bar(x, nonref_values, width=bar_width, label="Nonref", color="blue", alpha=0.7)
    plt.bar([i + bar_width for i in x], ref_values, width=bar_width, label="Ref", color="orange", alpha=0.7)
    
    # Formatting
    plt.xlabel("Country")
    plt.ylabel("Count")
    plt.title("Country Occurrences in Rows with Gene Symbols (Nonref vs Ref)")
    plt.xticks([i + bar_width / 2 for i in x], countries, rotation=45, ha="right")
    plt.legend()
    plt.tight_layout()
    
    # Save plot
    plt.savefig(output_file)
    plt.close()

# Main execution
if __name__ == "__main__":
    df_nonref, df_ref = load_data('../genes_overlap_tes_nonref_no_duplicates.csv', '../genes_overlap_tes_ref_no_duplicates.csv')
    
    df_nonref = process_dataframe(df_nonref)
    df_ref = process_dataframe(df_ref)
    
    nonref_counts = count_countries(df_nonref)
    ref_counts = count_countries(df_ref)
    
    plot_country_occurrences(nonref_counts, ref_counts, 'country_occurrences.png')
