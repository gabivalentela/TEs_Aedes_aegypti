import pandas as pd
import matplotlib.pyplot as plt
import re
import polars as pl
import os
import numpy as np

# define main class name
class frequency_analysis:
    def count_TE_occurence(general_df):
        """
        Counts the occurrences of TEs and identifies the top 5 most frequent ones.

        Args:
            general_df (pl.DataFrame): Polars DataFrame containing TE data.

        Returns:
            pl.DataFrame: DataFrame with the top 5 TEs sorted by their occurrence count.
        """
        general_df = general_df.with_columns(pl.col("TE_name").str.to_uppercase())
        
        possible_classes = ['LINE', 'SINE', 'LTR', 'TIR', 'MITE', 'UNKNOWN']

        # Extract clean TE names - BASICALLY CHECK WHICH OF THE POSSIBLE CLASSES IS PRESENT IN THE TE NAME AND KEEP ONLY THAT
        df_te = general_df.with_columns(
            pl.col("TE_name").apply(lambda x: next((cls for cls in possible_classes if cls in x), x)).alias("TE_name")
        )
        grouped = df_te.group_by('TE_name').len()
        sorted_grouped = grouped.sort('len',  descending=True)
        
        return sorted_grouped

    def calculate_frequency_TE(df, country):
        """
        For each TE insertion (row), calculate the derived allele frequency
        across all sample columns.

        Genotype encoding:
            1/1        -> 2 derived alleles
            0/1 / 1/0  -> 1 derived allele  (heterozygous, phasing-agnostic)
            0/0        -> 0 derived alleles
            anything else (missing, .'/.' etc.) -> excluded from denominator
        """
        base_cols = ['chromosome', 'position_start', 'position_end', 'TE_name', 'strand']
        meta_cols     = [c for c in base_cols if c in df.columns] + ["te_class", "te_size", "size_class"]
        genotype_cols = [c for c in df.columns if c not in meta_cols]

        df = df.copy()

        # calculate the frequencies
        df['count_1_1'] = (df[genotype_cols] == "1/1").sum(axis=1)
        df['count_0_1'] = (df[genotype_cols] == "0/1").sum(axis=1)
        df['count_0_0'] = (df[genotype_cols] == "0/0").sum(axis=1)
        df['called'] = (df['count_1_1'] + df['count_0_1'] + df['count_0_0'])*2
        df['frequency'] = (2 * df['count_1_1'] + df['count_0_1']) / (df['called'])
        df['frequency_0'] = (2 * df['count_0_0'] + df['count_0_1']) / (df['called'])

        # remove rows were all rows are missing (no called genotypes)
        df = df[df['called'] > 0]

        # remove rows where all genotypes are fixed (either all 1/1 or all 0/0)
        df = df[df['frequency'] < 1]
        df = df[df['frequency_0'] < 1]

        # drop the intermediate count columns
        df = df.drop(columns=['count_1_1', 'count_0_1', 'count_0_0', 'called', 'frequency_0'])

        # getting the ones in the last bin equal of above 0.9
        high_freq_tes = df[df['frequency'] >= 0.9]

        return df, high_freq_tes

    def getting_sfs_data(df):
        # group by frequency and count how many TEs have that frequency
        sfs_data = df.groupby('frequency').size().reset_index(name='count')
        sfs_data = sfs_data.sort_values(by='frequency')
        return sfs_data
    def plot_binned_sfs(sfs_data_by_country, type_df):
        country_colors = {
            "Brazil":   "#60D08F",  # Teal Green
            "USA":      "#785EF0",  # Strong Red
            "Colombia": "#DC267F",  # Orange
            "Senegal":  "#FE6100",  # Blue
            "Gabon":    "#FFB000",  # Vermillion
            "Kenya":    "#648FFF"   # Reddish Purple
        }

        countries = [c for c in country_colors if c in sfs_data_by_country]

        n_bins = 10
        bin_edges = np.linspace(0, 1.0, n_bins + 1)

        bar_width = 0.1
        fig, ax = plt.subplots(figsize=(14, 6))

        x_bins = np.arange(n_bins)
        total_countries = len(countries)

        for i, country in enumerate(countries):
            sfs_data = sfs_data_by_country[country].copy()

            # Assign each TE's frequency to its correct range-based bin
            sfs_data['bin'] = pd.cut(
                sfs_data['frequency'], bins=bin_edges,
                labels=False, include_lowest=True
            )

            # Sum counts within each bin
            counts = (
                sfs_data.groupby('bin')['count']
                .sum()
                .reindex(range(n_bins), fill_value=0)
            )

            # Normalize to proportion of TEs in each bin
            counts = counts / counts.sum()

            offset = (i - total_countries / 2) * bar_width * 1.5
            x_te = x_bins + offset

            ax.bar(
                x_te, counts.values,
                width=bar_width, color=country_colors[country],
                label=f"{country}"
            )

        ax.set_xticks(x_bins)
        ax.set_xticklabels([f"{bin_edges[i]:.1f}-{bin_edges[i+1]:.1f}" for i in range(n_bins)])
        ax.tick_params(axis='x', labelsize=15)
        ax.tick_params(axis='y', labelsize=15)
        ax.set_xlabel('Frequency', fontsize=15)
        ax.set_ylabel('Proportion of TEs', fontsize=15)
        ax.set_title(f"Binned Site Frequency Spectrum - {type_df}")
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')

        plt.tight_layout()
        plt.savefig(f'binned_SiteFrequencySpectrum_{type_df}.svg')