import pandas as pd
import matplotlib.pyplot as plt
import polars as pl
import numpy as np
import math

NUM_BINS = 10  # number of frequency bins between 0 and 1


def count_TE_occurence(general_df):
    """Identify the broad TE class present in each TE_name and count occurrences."""
    general_df = general_df.with_columns(pl.col("TE_name").str.to_uppercase())

    possible_classes = ['LINE', 'SINE', 'LTR', 'TIR', 'MITE', 'UNKNOWN', 'HELITRON', 'PENELOPE']

    df_te = general_df.with_columns(
        pl.col("TE_name").map_elements(
            lambda x: next((cls for cls in possible_classes if cls in x.upper()), x),
            return_dtype=pl.Utf8
        ).alias("TE_name")
    )

    grouped = df_te.group_by('TE_name').len()
    sorted_grouped = grouped.sort('len', descending=True)

    return sorted_grouped


def split_by_te_type_ref(top_te, country_df):
    """Return only the rows whose TE_name matches the given TE class."""
    country_df = country_df.with_columns(pl.col('TE_name').str.split('|').list.get(0).str.to_uppercase())
    mask = country_df[:, 3].str.contains(top_te)
    selected_df = country_df.filter(mask).clone()
    return selected_df

def calculate_frequency_TE(df, base_cols):
    """
    For each TE insertion (row), calculate the derived allele frequency
    across all sample columns.
    Genotype encoding:
        1/1        -> 2 derived alleles
        0/1 / 1/0  -> 1 derived allele  (heterozygous, phasing-agnostic)
        0/0        -> 0 derived alleles
        anything else (missing, .'/.' etc.) -> excluded from denominator
    """
    meta_cols     = [c for c in base_cols if c in df.columns] + ["te_class", "te_size", "size_class"]
    genotype_cols = [c for c in df.columns if c not in meta_cols]
    df = df.copy()
    # calculate the frequencies
    df['count_1_1'] = (df[genotype_cols] == "1/1").sum(axis=1)
    df['count_0_1'] = (df[genotype_cols] == "0/1").sum(axis=1)
    df['count_0_0'] = (df[genotype_cols] == "0/0").sum(axis=1)
    df['called'] = (df['count_1_1'] + df['count_0_1'] + df['count_0_0'])*2
    df['frequency_1'] = (2 * df['count_1_1'] + df['count_0_1']) / (df['called'])
    df['frequency_0'] = (2 * df['count_0_0'] + df['count_0_1']) / (df['called'])
    # remove rows were all rows are missing (no called genotypes)
    df = df[df['called'] > 0]
    # remove rows where all genotypes are fixed (either all 1/1 or all 0/0)
    df = df[df['frequency_1'] < 1]
    df = df[df['frequency_0'] < 1]
    # drop the intermediate count columns
    df = df.drop(columns=['count_1_1', 'count_0_1', 'count_0_0', 'called', 'frequency_0'])
    return df

def bin_frequency_data(frequency_df, top_te, num_bins=NUM_BINS):
    """
    Bin the allele frequencies (0-1) for a single TE type into fixed-width bins
    and return the proportion of TEs falling into each bin.
    """
    bins = np.linspace(0, 1, num_bins + 1)
    bin_labels = [f"{bins[i]:.1f}-{bins[i + 1]:.1f}" for i in range(num_bins)]

    if frequency_df.empty:
        return pd.DataFrame({'bin': bin_labels, top_te: [0.0] * num_bins})

    binned = pd.cut(frequency_df['frequency_1'], bins=bins, labels=bin_labels, include_lowest=True)
    counts = binned.value_counts().reindex(bin_labels, fill_value=0)
    total = counts.sum()
    proportions = counts / total if total > 0 else counts.astype(float)

    result_df = pd.DataFrame({'bin': bin_labels, top_te: proportions.values})
    return result_df


def sfs_general_data(dataframes):
    """Merge the per-TE-type binned frequency dataframes into one master dataframe indexed by bin."""
    dataframes = [d for d in dataframes if not d.empty]
    if not dataframes:
        return pd.DataFrame()

    master_df = dataframes[0][['bin']].copy()
    for d in dataframes:
        master_df = master_df.merge(d, on='bin', how='left')

    master_df = master_df.set_index('bin').fillna(0)
    return master_df


def plot_sfs_master_df(ax, df, suf, country):
    """Plot the binned SFS master dataframe as a grouped bar plot."""
    if df.empty:
        ax.set_visible(False)
        return

    variable_columns = df.columns.tolist()
    n_vars = len(variable_columns)
    bar_width = 0.8 / max(n_vars, 1)

    x = np.arange(len(df.index))
    colormap = plt.get_cmap('Blues', n_vars + 2)

    for i, col in enumerate(variable_columns):
        ax.bar(x + (i - n_vars / 2) * bar_width, df[col], bar_width,
               label=col, color=colormap(i + 1), edgecolor='black')

    ax.set_xlabel('Allele frequency bin', fontsize=11)
    ax.set_ylabel('Proportion of TEs', fontsize=11)
    ax.set_xticks(x)
    ax.set_xticklabels(df.index, rotation=45, ha='right', fontsize=8)
    ax.set_title(f"{country} - {suf}", fontsize=13)

    legend = ax.legend(loc='upper right', fontsize=7)
    if legend:
        for item in legend.get_texts():
            item.set_fontsize(7)


# ---------------------------------------------------------------------------
# Load and prepare data
# ---------------------------------------------------------------------------

reference_path = '/work/users/g/a/gabivla/lab/SV_mosquitos/TE_paper_aedesaegytpi/reference_vcf/genotyping_reference_TEs/post_processing/merged_results_5_2_either_filtered_keep_small_TEs.csv'
df = pd.read_csv(reference_path)
# drop size column if it exists
df = df.drop(columns=['size'], errors='ignore')

location_path = '/nas/longleaf/home/gabivla/lab/SV_mosquitos/location_information/mcclintock_output_locations.csv'
location_df = pd.read_csv(location_path, usecols=['Genome_code', 'country'])
location_dict = dict(zip(location_df['Genome_code'], location_df['country']))

base_cols = ['chromosome', 'position_start', 'position_end', 'TE_name', 'strand', 'score']

genome_cols = [col for col in df.columns if col not in base_cols]

# standardize genome column names
for genome in genome_cols:
    if genome.startswith('SRR'):
        genome_updated = genome.replace('_1', '')
    elif genome.startswith('FEMALE') or genome.startswith('MALE'):
        genome_updated = genome.replace('_L001_L002_R1_001', '')
    else:
        genome_updated = genome
    df = df.rename(columns={genome: genome_updated})

genome_cols = [col for col in df.columns if col not in base_cols]

rename_dict = {col: f"{col}_{location_dict.get(col, '')}" for col in genome_cols}
df = df.loc[:, ~df.columns.str.contains("Unnamed")]
df = df.rename(columns=rename_dict)

countries = set()
for original_col in genome_cols:
    new_col = rename_dict[original_col]
    country = new_col.split('_')[-1]
    countries.add(country)
countries = sorted(countries)

df_pl = pl.DataFrame(df)
TE_types = count_TE_occurence(df_pl)
TE_types = TE_types['TE_name'].to_list()

# ---------------------------------------------------------------------------
# Main loop: one figure per TE length category, one subplot per country
# ---------------------------------------------------------------------------

te_length_categories = ['small', 'medium', 'large']

with open('150_df_sizes.txt', 'w') as outfile:

    for te_length_category in te_length_categories:
        print(f"\n--- Processing category: {te_length_category} ---")

        n_countries = len(countries)
        n_cols = 2
        n_rows = math.ceil(n_countries / n_cols)

        fig, axs = plt.subplots(nrows=n_rows, ncols=n_cols,
                                 figsize=(6 * n_cols, 4.5 * n_rows),
                                 constrained_layout=True)
        axs = np.array(axs).reshape(n_rows, n_cols)

        for i, country in enumerate(countries):
            country_cols = [col for col in df.columns
                             if col not in base_cols and col.endswith('_' + country)]

            row, col_idx = divmod(i, n_cols)
            ax = axs[row][col_idx]

            if not country_cols:
                ax.set_visible(False)
                continue

            country_df = df[base_cols + country_cols].copy()
            country_df = country_df.dropna(subset=country_cols, how='all')

            print(f"  Processing {country}")

            country_df['TE_length'] = (country_df['position_end'] - country_df['position_start']).astype(int)

            bins = [float('-inf'), 150, 600, float('inf')]
            labels = ['small', 'medium', 'large']
            country_df['TE_length_category'] = pd.cut(country_df['TE_length'], bins=bins, labels=labels, right=False)

            te_length_df = country_df[country_df['TE_length_category'] == te_length_category]

            outfile.write(f"Size of the dataframe for {country} {te_length_category}: {len(te_length_df)}\n")

            if len(te_length_df) == 0:
                print(f"    No {te_length_category} data for {country}")
                ax.set_visible(False)
                continue

            list_of_binned_dfs = []

            for top_te in TE_types:
                pl_te_df = pl.DataFrame(te_length_df.drop(columns=['TE_length', 'TE_length_category']))
                te_content_all_genomes = split_by_te_type_ref(top_te, pl_te_df)
                te_content_all_genomes = te_content_all_genomes.to_pandas()

                frequency_TEs = calculate_frequency_TE(te_content_all_genomes, base_cols)
                binned_df = bin_frequency_data(frequency_TEs, top_te, num_bins=NUM_BINS)
                list_of_binned_dfs.append(binned_df)

            master_df_ref = sfs_general_data(list_of_binned_dfs)

            if master_df_ref.empty:
                print(f"    No TE data for {country} {te_length_category}, skipping...")
                ax.set_visible(False)
                continue

            print(master_df_ref)

            plot_sfs_master_df(ax, master_df_ref, 'ref', country)

        # hide any leftover unused axes
        for j in range(n_countries, n_rows * n_cols):
            row, col_idx = divmod(j, n_cols)
            axs[row][col_idx].set_visible(False)

        fig.suptitle(f'Site Frequency Spectrum (binned) - {te_length_category} TEs', fontsize=16)

        plt.savefig(f'150_SiteFrequencySpectrum_reference_{te_length_category}.svg')
        plt.close(fig)