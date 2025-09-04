import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

'''# Read your data
df_ref = pd.read_csv('./TEs_overlapping_genes_ref.txt')
print(df_ref)

# Filter columns that start with 'count_1/1_'
count_columns = [col for col in df_ref.columns if col.startswith('count_1/1_')]

# Melt the DataFrame to have 'Country' as a variable
df_melted = df_ref.melt(id_vars=['TE_information'], value_vars=count_columns,
                        var_name='Country', value_name='Count')

# Remove 'count_1/1_' from the country name
df_melted['Country'] = df_melted['Country'].str.replace('count_1/1_', '')

# Plot one plot for each insertion (TE)
unique_TE = df_melted['TE_information'].unique()

for te in unique_TE:
    plt.figure(figsize=(10, 6))
    te_data = df_melted[df_melted['TE_information'] == te]
    sns.barplot(data=te_data, x='Country', y='Count', color='lightblue')
    plt.title(f'{te.replace("_", " ")}', fontsize=20)    
    plt.ylabel('Count')
    plt.xlabel('Country')
    plt.xticks(rotation=45, fontsize=25)
    plt.tight_layout()
    plt.savefig(f'{te}_ref.svg')'''

# total genomes per country
total_number_of_genomes = {
    'Brazil': 18, 
    'Colombia': 24, 
    'USA': 28, 
    'Senegal': 20, 
    'Gabon': 13, 
    'Kenya': 19
}

# Read data
df_nonref = pd.read_csv('./TEs_overlapping_genes_nonref.txt')

# Find all count columns
count_columns = [col for col in df_nonref.columns if col.startswith('count_1/1_')]

# Make proportion columns
for col in count_columns:
    # extract the country name from the column name
    country = col.replace("count_1/1_", "")
    df_nonref[f"Proportion_{country}"] = df_nonref[col] / total_number_of_genomes[country]

# Save results
df_nonref.to_csv("information_plotting_frequency.csv", index=False)

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# total genomes per country
total_number_of_genomes = {
    'Brazil': 18, 
    'Colombia': 24, 
    'USA': 28, 
    'Senegal': 20, 
    'Gabon': 13, 
    'Kenya': 19
}

# Melt the DataFrame
df_melted = df_nonref.melt(
    id_vars=['TE_information'], 
    value_vars=count_columns,
    var_name='Country', 
    value_name='Count'
)

# Clean up country names
df_melted['Country'] = df_melted['Country'].str.replace('count_1/1_', '')

# Add proportion column
df_melted['Proportion'] = df_melted.apply(
    lambda row: row['Count'] / total_number_of_genomes[row['Country']], axis=1
)

# Plot one plot for each TE, but using Proportion
unique_TE = df_melted['TE_information'].unique()

for te in unique_TE:
    plt.figure(figsize=(15, 8))
    te_data = df_melted[df_melted['TE_information'] == te]
    
    # use proportion instead of count
    sns.barplot(data=te_data, x='Country', y='Proportion', color='lightgray')
    
    plt.title(f'{te.replace("_", " ")}', fontsize=25)    
    plt.ylabel('Proportion of genomes with TE', fontsize=25)
    plt.xlabel('Country', fontsize=25)
    plt.yticks(fontsize=25)
    plt.xticks(rotation=45, fontsize=25)
    plt.tight_layout()
    plt.savefig(f'{te}_nonref_proportion.svg')
