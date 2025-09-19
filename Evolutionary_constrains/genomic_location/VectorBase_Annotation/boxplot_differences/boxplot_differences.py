import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns

# Load nonref and ref datasets
nonref = '../../../nonreference/results'
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
        df['Category'] = 'NonRef'  # Add category column
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
        df['Category'] = 'Ref'  # Add category column
        list_df_ref.append(df)

ref_df = pd.concat(list_df_ref)

# Combine both datasets
combined_df = pd.concat([nonref_df, ref_df])

# Create subplots
fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=True)

# Plot NonRef data
sns.boxplot(y=ref_df['frequency_diff'], ax=axes[0], color='lightblue')
axes[0].set_title('Ref')
axes[0].set_ylabel('Frequency Difference')

# Plot Ref data
sns.boxplot(y=nonref_df['frequency_diff'], ax=axes[1], color='lightcoral')
axes[1].set_title('NonRef')

# Improve layout
plt.tight_layout()
plt.savefig('boxplot_differences.svg')

mean_nonref = nonref_df["frequency_diff"].mean()
mean_ref = ref_df["frequency_diff"].mean()
print(f"Mean Frequency Difference (NonRef): {mean_nonref}")
print(f"Mean Frequency Difference (Ref): {mean_ref}")