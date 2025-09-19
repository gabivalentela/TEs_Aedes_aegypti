import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns

# Load nonref and ref datasets
nonref = '../../../nonreference/results/'
ref = '../../../reference/results/'

# List files
nonref_files = os.listdir(nonref)
ref_files = os.listdir(ref)

# Load nonref data
list_df_nonref = []
for file in nonref_files:
    if file.endswith('.csv'):
        path = os.path.join(nonref, file)
        df = pd.read_csv(path)
        # Keep only relevant columns
        df = df[['TE_information', 'frequency_diff']]
        df['Category'] = 'NonRef'  
        list_df_nonref.append(df)

nonref_df = pd.concat(list_df_nonref)

# Load ref data
list_df_ref = []
for file in ref_files:
    if file.endswith('.csv'):
        path = os.path.join(ref, file)
        df = pd.read_csv(path)
        # Keep only relevant columns
        df = df[['TE_information', 'frequency_diff']]
        df['Category'] = 'Ref'  
        list_df_ref.append(df)

ref_df = pd.concat(list_df_ref)

# Extract TE names
nonref_df['TE_name'] = nonref_df['TE_information'].str.extract(r'AaegL5_\d_\d+_\d+_rnd_\d+_family_\d+_(.+)_.+')
ref_df['TE_name'] = ref_df['TE_information'].str.extract(r'AaegL5_\d_\d+_\d+_rnd_\d+_family_\d+_(.+)_.+')

# Count TE names
nonref_count = nonref_df['TE_name'].value_counts()
ref_count = ref_df['TE_name'].value_counts()

# Plot TE counts separately
fig, axes = plt.subplots(1, 2, figsize=(15, 6))

# NonRef TE count plot
sns.barplot(x=ref_count.index[:10], y=ref_count.values[:10], ax=axes[0], color='lightblue')
axes[0].set_title('Top 10 TE Counts (Ref)')
axes[0].set_ylabel('Count')
axes[0].set_xlabel('TE Name')
axes[0].tick_params(axis='x', rotation=90)

# Ref TE count plot
sns.barplot(x=nonref_count.index[:10], y=nonref_count.values[:10], ax=axes[1], color='lightcoral')
axes[1].set_title('Top 10 TE Counts (NonRef)')
axes[1].set_ylabel('Count')
axes[1].set_xlabel('TE Name')
axes[1].tick_params(axis='x', rotation=90)

# Adjust layout
plt.tight_layout()
plt.savefig('TE_counts_ref_nonref.svg')
