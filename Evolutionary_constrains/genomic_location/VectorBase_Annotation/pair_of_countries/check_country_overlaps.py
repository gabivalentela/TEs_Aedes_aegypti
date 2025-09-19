import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import re
from collections import Counter

# Function to load and merge CSV files from a directory
def load_and_merge_te_info(directory, pattern):
    csv_files = [f for f in os.listdir(directory) if f.endswith('.csv')]
    dfs = []
    for file in csv_files:
        match = re.search(pattern, file)
        if match:
            country_combo = match.group(1)
            df = pd.read_csv(os.path.join(directory, file))
            df['country_combo'] = country_combo
            df = df[['country_combo', 'TE_information']]
            dfs.append(df)
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()

# Function to merge TE information with main dataset
def merge_te_data(main_df, te_df):
    return main_df.merge(te_df, on="TE_information", how="left")

# Function to count TEs associated with a gene per country combination
def count_te_with_gene(df):
    return df[df["gene_symbol"].notna()].groupby("country_combo")["TE_information"].count().reset_index()

# Function to plot bar chart
def plot_te_counts(te_counts, title, filename):
    plt.figure(figsize=(30, 25))
    sns.barplot(x="country_combo", y="num_TEs_with_gene", data=te_counts, palette="tab20")
    plt.xticks(rotation=45, ha="right", fontsize=25)
    plt.xlabel("Country Combination", fontsize=25)
    plt.yticks(fontsize=25)
    plt.ylabel("Number of TEs with Associated Gene", fontsize=25)
    plt.title(title)
    plt.savefig(filename)
    plt.close()

# Load main datasets
df_nonref = pd.read_csv('../genes_overlap_tes_nonref_no_duplicates.csv')
df_ref = pd.read_csv('../genes_overlap_tes_ref_no_duplicates.csv')

# Load and merge TE information
tes_nonref = load_and_merge_te_info('../../../nonreference/results/', r'top_10_diff_(.+).csv')
tes_ref = load_and_merge_te_info('../../../reference/results/', r'top_10_diff_(.+).csv')

df_nonref = merge_te_data(df_nonref, tes_nonref)
df_ref = merge_te_data(df_ref, tes_ref)

# Count and plot for nonreference
te_counts_nonref = count_te_with_gene(df_nonref)
te_counts_nonref.columns = ["country_combo", "num_TEs_with_gene"]
plot_te_counts(te_counts_nonref, "Number of TEs with Gene per Country Combination (Nonreference)", 'num_TEs_with_gene_per_country_combo_nonreference.png')

# Count and plot for reference
te_counts_ref = count_te_with_gene(df_ref)
te_counts_ref.columns = ["country_combo", "num_TEs_with_gene"]
plot_te_counts(te_counts_ref, "Number of TEs with Gene per Country Combination (Reference)", 'num_TEs_with_gene_per_country_combo_reference.png')
