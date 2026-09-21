import re
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


# ============================================================
# Input SLURM output files
# ============================================================

reference = '../../reference/slurm-458158.out'
non_reference = '../../nonreference/slurm-458321.out'
SNPs_data = '../../SNPs/slurm-477482.out'


# ============================================================
# Desired country order
# ============================================================

country_order = [
    'Brazil',
    'Colombia',
    'USA',
    'Gabon',
    'Kenya',
    'Senegal'
]


# ============================================================
# Parse SLURM output
# ============================================================

def parse_slurm_fst(slurm_file):

    fst_data = {}

    pattern = (
        r'^Genome-wide weighted FST\s*'
        r'\((.+?)\s+vs\s+(.+?)\):\s*([0-9.eE+-]+)'
    )

    with open(slurm_file, 'r') as file:

        for line in file:
            if 'Unknown' in line:
                continue

            match = re.search(pattern, line.strip())

            if match:

                country_1, country_2, fst_value = match.groups()

                country_1 = country_1.strip()
                country_2 = country_2.strip()

                fst_value = float(fst_value)

                # Make pair independent of order
                pair = tuple(sorted([country_1, country_2]))

                fst_data[pair] = fst_value
    print(fst_data)
    return fst_data


# ============================================================
# Make FST dataframe
# ============================================================

def make_fst_dataframe(fst_data, country_order):

    # Use NaN for missing values
    df = pd.DataFrame(
        np.nan,
        index=country_order,
        columns=country_order
    )

    for (country_1, country_2), fst_value in fst_data.items():

        if country_1 not in country_order:
            print(f'Warning: {country_1} not in country order')
            continue

        if country_2 not in country_order:
            print(f'Warning: {country_2} not in country order')
            continue

        # Fill both sides of the matrix
        df.loc[country_1, country_2] = fst_value
        df.loc[country_2, country_1] = fst_value

    return df


# ============================================================
# Plot triangular heatmap
# ============================================================

def plot_fst_heatmap(df, name, output_file):

    # Mask upper triangle
    mask = np.triu(
        np.ones_like(df, dtype=bool)
    )

    plt.figure(figsize=(8, 6))

    sns.heatmap(
        df,
        annot=True,
        fmt='.3f',
        cmap='plasma',
        square=True,
        mask=mask,
        linewidths=0.5,
        linecolor='white',
        vmin=0,
        vmax=0.5,
        cbar_kws={'label': f'{name} FST'}
    )

    plt.title(f'Pairwise FST - {name}')

    plt.xlabel('Country')
    plt.ylabel('Country')

    plt.tight_layout()

    plt.savefig(
        output_file,
        format='svg',
        bbox_inches='tight'
    )

    plt.close()


# ============================================================
# Parse SLURM files
# ============================================================
snp_fst = parse_slurm_fst(SNPs_data)
reference_fst = parse_slurm_fst(reference)
non_reference_fst = parse_slurm_fst(non_reference)


# ============================================================
# Print values
# ============================================================

print('\nSNP SLURM FST values:')

for pair, value in snp_fst.items():
    print(f'{pair}: {value}')

print('\nReference SLURM FST values:')

for pair, value in reference_fst.items():
    print(f'{pair}: {value}')


print('\nNon-reference SLURM FST values:')

for pair, value in non_reference_fst.items():
    print(f'{pair}: {value}')


# ============================================================
# Create dataframes
# ============================================================

snp_df = make_fst_dataframe(
    snp_fst,
    country_order
)

reference_df = make_fst_dataframe(
    reference_fst,
    country_order
)

non_reference_df = make_fst_dataframe(
    non_reference_fst,
    country_order
)


# ============================================================
# Print matrices
# ============================================================

print('\nSNP FST matrix:')
print(snp_df)

print('\nReference FST matrix:')
print(reference_df)

print('\nNon-reference FST matrix:')
print(non_reference_df)


# ============================================================
# Plot heatmaps
# ============================================================

plot_fst_heatmap(
    snp_df,
    'SNPs',
    'heatmap_SNPs_slurm.svg'
)

plot_fst_heatmap(
    reference_df,
    'Reference',
    'heatmap_reference_slurm.svg'
)

plot_fst_heatmap(
    non_reference_df,
    'Non-reference',
    'heatmap_nonreference_slurm.svg'
)


print('\nHeatmaps saved:')
print('  heatmap_SNPs_slurm.svg')
print('  heatmap_reference_slurm.svg')
print('  heatmap_nonreference_slurm.svg')