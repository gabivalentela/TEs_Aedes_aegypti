"""
Plot the TE landscape across chromosomes, in the same style as the MAPQ plots.

For each chromosome, bins TEs into fixed-size windows and plots TE density
(count per window) and mean allele frequency per window vs. genomic position,
as line plots stacked one-subplot-per-chromosome — the same layout as the
MAPQ script.

Unlike the original version (which overlaid all TE classes on one set of
axes), this version saves a SEPARATE file for the total and for each
individual TE class, for both density and frequency. Save these images next
to your MAPQ plots and eyeball them side-by-side (or overlay them) to check
whether TE-dense regions line up with low-MAPQ regions.

Requirements:
    pip install pandas numpy matplotlib

Notes:
    - WINDOW_SIZE defaults to 100 kb to match the MAPQ script's default, so
      x-axis bins line up 1:1 between the two sets of plots. Change both
      scripts together if you adjust resolution.
    - `chromosomes` is hardcoded to match the MAPQ script's chromosome list
      and order, again so the panels line up when compared side-by-side.
    - Output files are named:
          TE_landscape_density_<CLASS>_<type_df>.png
          TE_landscape_frequency_<CLASS>_<type_df>.png
      where <CLASS> is "total" or one of the TE classes (LINE, SINE, LTR,
      MITE, TIR, HELITRON, PENELOPE, UNKNOWN, OTHER).
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

WINDOW_SIZE = 100_000  # bp per bin — keep in sync with the MAPQ script
OUT_DIR = './'
CHROMOSOMES = ['AaegL5_1', 'AaegL5_2', 'AaegL5_3']

TE_ORDER = ['LINE', 'SINE', 'LTR', 'MITE', 'TIR', 'HELITRON', 'PENELOPE', 'UNKNOWN']
ALL_CLASSES = TE_ORDER + ['OTHER']  # OTHER catches names classify_te() couldn't match


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
    """Bin TEs into fixed windows for one chromosome.

    Returns:
        overall       - total count & mean frequency per window
        by_class_count - count per window, one column per TE class
        by_class_freq  - mean frequency per window, one column per TE class
    """
    df_chrom = df_chrom.copy()
    df_chrom['TE_class'] = df_chrom['TE_name'].apply(classify_te)
    df_chrom['window'] = df_chrom['position_start'] // window_size
    df_chrom = calculate_frequency(df_chrom)

    max_window = int(df_chrom['window'].max())
    full_index = range(max_window + 1)

    overall = (
        df_chrom.groupby('window')
        .agg(count=('frequency', 'count'), mean_frequency=('frequency', 'mean'))
        .reindex(full_index, fill_value=0)
    )
    overall.loc[overall['count'] == 0, 'mean_frequency'] = np.nan

    by_class_count = (
        df_chrom.groupby(['window', 'TE_class'])
        .size()
        .unstack(fill_value=0)
        .reindex(full_index, fill_value=0)
    )

    by_class_freq = (
        df_chrom.groupby(['window', 'TE_class'])['frequency']
        .mean()
        .unstack()
        .reindex(full_index)
    )

    return overall, by_class_count, by_class_freq


def compute_landscape_data(df, chromosomes, window_size):
    """Precompute per-chromosome window data once, reused across all output files."""
    data = {}
    for chrom in chromosomes:
        df_chrom = df[df['chromosome'] == chrom]
        if df_chrom.empty:
            data[chrom] = None
            continue
        overall, by_class_count, by_class_freq = compute_te_windows(df_chrom, window_size)
        data[chrom] = {
            'overall': overall,
            'count': by_class_count,
            'freq': by_class_freq,
        }
    return data


def _extract_series(cd, value_key, te_class):
    """Pull the plotting series for one chromosome's worth of precomputed data."""
    if te_class is None:  # total
        return cd['overall']['count'] if value_key == 'count' else cd['overall']['mean_frequency']

    table = cd[value_key]
    if te_class in table.columns:
        return table[te_class]
    fill = 0 if value_key == 'count' else np.nan
    return pd.Series(fill, index=cd['overall'].index)


def plot_landscape_panel(chrom_data, chromosomes, window_size, out_path, title,
                          value_key, te_class, color, ylabel, ylim=None):
    """Render and save one stacked, one-subplot-per-chromosome figure."""
    share_y = value_key == 'freq'
    fig, axes = plt.subplots(len(chromosomes), 1, figsize=(12, 3 * len(chromosomes)),
                              sharey=share_y)
    if len(chromosomes) == 1:
        axes = [axes]

    for ax, chrom in zip(axes, chromosomes):
        cd = chrom_data.get(chrom)
        if cd is None:
            ax.set_title(f'{chrom} (no data)', fontsize=10)
            continue

        series = _extract_series(cd, value_key, te_class)
        positions_mb = cd['overall'].index.values * window_size / 1e6

        ax.plot(positions_mb, series.values, lw=0.9, color=color)
        ax.set_title(chrom, fontsize=10)
        ax.set_ylabel(ylabel)
        if ylim is not None:
            ax.set_ylim(*ylim)
        ax.grid(True, linewidth=0.3, alpha=0.4)

    axes[-1].set_xlabel('Position (Mb)')

    fig.suptitle(title, fontsize=12)
    plt.tight_layout()

    plt.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f'Saved {out_path}')


def plot_te_landscapes(df, chromosomes, window_size, out_dir, type_df):
    """Save one density file and one frequency file for the total, plus one
    density file and one frequency file for each individual TE class."""
    chrom_data = compute_landscape_data(df, chromosomes, window_size)
    class_colors = dict(zip(TE_ORDER, cm.tab10(np.linspace(0, 1, len(TE_ORDER)))))
    class_colors['OTHER'] = 'gray'

    # --- Total, across all TE classes ---
    plot_landscape_panel(
        chrom_data, chromosomes, window_size,
        f'{out_dir}/TE_landscape_density_total_{type_df}.png',
        f'TE Landscape (density) — Total — {type_df.capitalize()}',
        value_key='count', te_class=None, color='black',
        ylabel='TE count / window',
    )
    plot_landscape_panel(
        chrom_data, chromosomes, window_size,
        f'{out_dir}/TE_landscape_frequency_total_{type_df}.png',
        f'TE Landscape (allele frequency) — Total — {type_df.capitalize()}',
        value_key='freq', te_class=None, color='darkorange',
        ylabel='Mean TE frequency', ylim=(0, 1),
    )

    # --- One file per individual TE class ---
    for cls in ALL_CLASSES:
        color = class_colors[cls]
        plot_landscape_panel(
            chrom_data, chromosomes, window_size,
            f'{out_dir}/TE_landscape_density_{cls}_{type_df}.png',
            f'TE Landscape (density) — {cls} — {type_df.capitalize()}',
            value_key='count', te_class=cls, color=color,
            ylabel='TE count / window',
        )
        plot_landscape_panel(
            chrom_data, chromosomes, window_size,
            f'{out_dir}/TE_landscape_frequency_{cls}_{type_df}.png',
            f'TE Landscape (allele frequency) — {cls} — {type_df.capitalize()}',
            value_key='freq', te_class=cls, color=color,
            ylabel='Mean TE frequency', ylim=(0, 1),
        )


# ── Main ─────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    reference_path = '/work/users/g/a/gabivla/lab/SV_mosquitos/TE_paper_aedesaegytpi/reference_vcf/genotyping_reference_TEs/post_processing/merged_results_5_2_either_filtered.csv'
    non_reference_path = '/work/users/g/a/gabivla/lab/SV_mosquitos/TE_paper_aedesaegytpi/nonreference_vcf/post_processing_nonreference_TEs/post_processing_overlap_df/max_min_updated_strategy/final_df_overlap_retroseq_teforest_nonreference_1000.csv'

    df_ref = pd.read_csv(reference_path).drop(columns=['score'])
    df_nonref = pd.read_csv(non_reference_path)

    for df, type_df in [(df_ref, 'reference'), (df_nonref, 'non-reference')]:
        plot_te_landscapes(df, CHROMOSOMES, WINDOW_SIZE, OUT_DIR, type_df)