import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def calculate_frequencies(df, numer_supporting_reads):
    count_11 = (df == '1/1').sum().sum()
    count_01 = (df == '0/1').sum().sum()
    count_00 = (df == '0/0').sum().sum()
    count_missing = (df == './.').sum().sum()

    counts = [count_11, count_01, count_00, count_missing]
    labels = ['1/1', '0/1', '0/0', './.']

    total_TEs = len(df)

    counts_df = pd.DataFrame({
        f'genotype_{numer_supporting_reads}': labels,
        'count': counts
    })

    counts_df.loc[len(counts_df.index)] = ['total_TEs', total_TEs]

    return counts_df


# Read data
path_2_5_either = '../post_processing/merged_results_5_2_either_filtered.csv'

df_2_5_either = pd.read_csv(path_2_5_either)

# Calculate counts
frequencies_either_supporting_reads = calculate_frequencies(
    df_2_5_either, 'either'
)

# ---------------------------------------
# Prepare data for plotting
# ---------------------------------------

genotypes = ['1/1', '0/1', '0/0', './.']

plot_df = pd.DataFrame({
    'Genotype': genotypes,
    'Either supporting reads': frequencies_either_supporting_reads.loc[:3, 'count'].values
})


# ---------------------------------------
# Calculate proportions
# ---------------------------------------

plot_proportion_df = plot_df.copy()

plot_proportion_df['Either supporting reads'] = (
    plot_proportion_df['Either supporting reads'] /
    len(df_2_5_either)
)

# ---------------------------------------
# Plot 1: Counts
# ---------------------------------------

x = np.arange(len(genotypes))
width = 0.25

fig, ax = plt.subplots(figsize=(10, 6))

ax.bar(
    x,
    plot_df['Either supporting reads'],
    width,
    label='Either supporting reads'
)

ax.set_xlabel('Genotype')
ax.set_ylabel('Count')
ax.set_title('Genotype counts by number of supporting reads')
ax.set_xticks(x)
ax.set_xticklabels(genotypes)
ax.legend()

plt.tight_layout()
plt.savefig('genotype_counts_by_supporting_reads.png')
plt.show()


# ---------------------------------------
# Plot 2: Proportions
# ---------------------------------------

fig, ax = plt.subplots(figsize=(10, 6))

ax.bar(
    x,
    plot_proportion_df['Either supporting reads'],
    width,
    label='Either supporting reads'
)

ax.set_xlabel('Genotype')
ax.set_ylabel('Proportion')
ax.set_title('Genotype proportions by number of supporting reads')
ax.set_xticks(x)
ax.set_xticklabels(genotypes)
ax.legend()

plt.tight_layout()
plt.savefig('genotype_proportions_by_supporting_reads.png')
plt.show()