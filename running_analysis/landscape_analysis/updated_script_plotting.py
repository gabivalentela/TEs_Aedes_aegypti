"""
Plot the TE landscape across chromosomes, in the same style as the MAPQ plots.

For each chromosome, bins TEs into fixed-size windows and plots TE density
(count per window) and mean allele frequency per window vs. genomic position,
as line plots stacked one-subplot-per-chromosome — the same layout as the
MAPQ script. Save these images next to your MAPQ plots and eyeball them
side-by-side (or overlay them) to check whether TE-dense regions line up
with low-MAPQ regions.

Requirements:
    pip install pandas numpy matplotlib

Notes:
    - WINDOW_SIZE defaults to 100 kb to match the MAPQ script's default, so
      x-axis bins line up 1:1 between the two sets of plots. Change both
      scripts together if you adjust resolution.
    - `chromosomes` is hardcoded to match the MAPQ script's chromosome list
      and order, again so the panels line up when compared side-by-side.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

WINDOW_SIZE = 100_000  # bp per bin — keep in sync with the MAPQ script
OUT_DIR = './'
CHROMOSOMES = ['AaegL5_1', 'AaegL5_2', 'AaegL5_3']

TE_ORDER = ['LINE', 'SINE', 'LTR', 'MITE', 'TIR', 'HELITRON', 'PENELOPE', 'UNKNOWN']


def classify_te(name):
    for cls in TE_ORDER:
        if cls in name:
            return cls
    return 'OTHER'


def calculate_frequency(df):
    base_columns = ['chromosome', 'position_start', 'position_end', 'TE_name',
                     'strand', 'TE_class', 'window']
    genotype_columns = [c for c in df.columns if c not in base_columns]
    df['count_11'] = df[genotype_columns].apply(lambda row: (row == '1/1').sum(), axis=1)
    df['count_01'] = df[genotype_columns].apply(lambda row: (row == '0/1').sum(), axis=1)
    df['frequency'] = (df['count_11'] * 2 + df['count_01']) / (2 * len(genotype_columns))
    return df


def compute_te_windows(df_chrom, window_size):
    """Bin TEs into fixed windows; return (overall counts/mean-freq per window,
    per-class counts per window)."""
    df_chrom = df_chrom.copy()
    df_chrom['TE_class'] = df_chrom['TE_name'].apply(classify_te)
    df_chrom['window'] = df_chrom['position_start'] // window_size
    df_chrom = calculate_frequency(df_chrom)

    max_window = int(df_chrom['window'].max())

    overall = (
        df_chrom.groupby('window')
        .agg(count=('frequency', 'count'), mean_frequency=('frequency', 'mean'))
        .reindex(range(max_window + 1), fill_value=0)
    )
    overall.loc[overall['count'] == 0, 'mean_frequency'] = np.nan

    by_class = (
        df_chrom.groupby(['window', 'TE_class'])
        .size()
        .unstack(fill_value=0)
        .reindex(range(max_window + 1), fill_value=0)
    )

    return overall, by_class


def plot_density_landscape(df, chromosomes, window_size, out_dir, type_df):
    """TE count per window across the genome — same panel layout as the MAPQ plot."""
    fig, axes = plt.subplots(len(chromosomes), 1, figsize=(12, 3 * len(chromosomes)))
    if len(chromosomes) == 1:
        axes = [axes]

    class_colors = cm.tab10(np.linspace(0, 1, len(TE_ORDER)))

    for ax, chrom in zip(axes, chromosomes):
        df_chrom = df[df['chromosome'] == chrom]
        if df_chrom.empty:
            ax.set_title(f'{chrom} (no data)', fontsize=10)
            continue

        overall, by_class = compute_te_windows(df_chrom, window_size)
        positions_mb = overall.index.values * window_size / 1e6

        ax.plot(positions_mb, overall['count'], lw=1.0, color='black',
                label='Total TEs', zorder=3)

        for cls, color in zip(TE_ORDER, class_colors):
            if cls in by_class.columns:
                ax.plot(positions_mb, by_class[cls].values, lw=0.6, alpha=0.7,
                         color=color, label=cls, zorder=2)

        ax.set_title(chrom, fontsize=10)
        ax.set_ylabel('TE count / window')
        ax.grid(True, linewidth=0.3, alpha=0.4)

    axes[-1].set_xlabel('Position (Mb)')

    handles, labels = axes[0].get_legend_handles_labels()
    axes[0].legend(handles, labels, loc='upper right', fontsize=7, ncol=4, framealpha=0.9)

    fig.suptitle(f'TE Landscape (density) — {type_df.capitalize()}', fontsize=12)
    plt.tight_layout()

    out_name = f'{out_dir}/TE_landscape_density_{type_df}.png'
    plt.savefig(out_name, dpi=150)
    plt.close(fig)
    print(f'Saved {out_name}')


def plot_frequency_landscape(df, chromosomes, window_size, out_dir, type_df):
    """Mean TE allele frequency per window — same panel layout as the MAPQ plot."""
    fig, axes = plt.subplots(len(chromosomes), 1, figsize=(12, 3 * len(chromosomes)),
                              sharey=True)
    if len(chromosomes) == 1:
        axes = [axes]

    for ax, chrom in zip(axes, chromosomes):
        df_chrom = df[df['chromosome'] == chrom]
        if df_chrom.empty:
            ax.set_title(f'{chrom} (no data)', fontsize=10)
            continue

        overall, _ = compute_te_windows(df_chrom, window_size)
        positions_mb = overall.index.values * window_size / 1e6

        ax.plot(positions_mb, overall['mean_frequency'], lw=0.8, color='darkorange')
        ax.set_title(chrom, fontsize=10)
        ax.set_ylabel('Mean TE frequency')
        ax.set_ylim(0, 1)
        ax.grid(True, linewidth=0.3, alpha=0.4)

    axes[-1].set_xlabel('Position (Mb)')

    fig.suptitle(f'TE Landscape (allele frequency) — {type_df.capitalize()}', fontsize=12)
    plt.tight_layout()

    out_name = f'{out_dir}/TE_landscape_frequency_{type_df}.png'
    plt.savefig(out_name, dpi=150)
    plt.close(fig)
    print(f'Saved {out_name}')


# ── Main ─────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    reference_path = '/work/users/g/a/gabivla/lab/SV_mosquitos/TE_paper_aedesaegytpi/reference_vcf/genotyping_reference_TEs/post_processing/merged_results_5_2_either_filtered.csv'
    non_reference_path = '/work/users/g/a/gabivla/lab/SV_mosquitos/TE_paper_aedesaegytpi/nonreference_vcf/post_processing_nonreference_TEs/post_processing_overlap_df/max_min_updated_strategy/final_df_overlap_retroseq_teforest_nonreference_1000.csv'

    df_ref = pd.read_csv(reference_path).drop(columns=['score'])
    df_nonref = pd.read_csv(non_reference_path)

    for df, type_df in [(df_ref, 'reference'), (df_nonref, 'non-reference')]:
        plot_density_landscape(df, CHROMOSOMES, WINDOW_SIZE, OUT_DIR, type_df)
        plot_frequency_landscape(df, CHROMOSOMES, WINDOW_SIZE, OUT_DIR, type_df)