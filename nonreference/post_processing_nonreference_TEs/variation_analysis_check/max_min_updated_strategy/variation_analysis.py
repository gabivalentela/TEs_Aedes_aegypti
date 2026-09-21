import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

def calculate_frequency(df):
    base_columns = ['chromosome', 'position_start', 'position_end', 'TE_name', 'strand']
    genotype_columns = [col for col in df.columns if col not in base_columns]
    df['count_1/1'] = df[genotype_columns].apply(lambda x: (x == '1/1').sum(), axis=1)
    df['count_0/0'] = df[genotype_columns].apply(lambda x: (x == '0/0').sum(), axis=1)
    df['frequency'] = df['count_1/1'] / (df['count_1/1'] + df['count_0/0'])
    return df

def get_te_family(te_name, possible_names):
    """Assign a TE family label based on the TE name string."""
    te_upper = str(te_name).upper()
    for family in possible_names:
        if family in te_upper:
            return family
    return 'UNKNOWN'

# ── Load data ──────────────────────────────────────────────────────────────────
list_files = os.listdir('../../post_processing_overlap_df/max_min_updated_strategy')
list_dfs = {}
for file in list_files:
    if file.endswith('.csv'):
        variation_name = file.split('.')[0].split('_')[-1]
        df = pd.read_csv('../../post_processing_overlap_df/max_min_updated_strategy/' + file)
        df = calculate_frequency(df)
        list_dfs[variation_name] = df

possible_names = ['LINE', 'SINE', 'LTR', 'MITE', 'TIR', 'HELITRON', 'PENELOPE', 'UNKNOWN']
variations     = list(list_dfs.keys())
colors         = plt.cm.tab10(np.linspace(0, 1, len(variations)))

# ── 1. Number of TEs per variation ────────────────────────────────────────────
def plot_te_count(list_dfs):
    counts = {var: len(df) for var, df in list_dfs.items()}

    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar(counts.keys(), counts.values(), color=colors, edgecolor='black')
    ax.bar_label(bars, padding=3, fontsize=10)
    ax.set_title('Number of TEs per Variation', fontsize=14, fontweight='bold')
    ax.set_xlabel('Variation')
    ax.set_ylabel('TE Count')
    ax.set_ylim(0, max(counts.values()) * 1.15)
    plt.tight_layout()
    plt.savefig('te_count_per_variation.png', dpi=150)
    plt.close()
    print("\n[1] TE counts:")
    for var, n in counts.items():
        print(f"  {var}: {n}")

# ── 2. TE size distribution per variation ─────────────────────────────────────
def plot_te_size(list_dfs):
    """Box-plot of TE length (position_end – position_start) per variation."""
    sizes = {
        var: (df['position_end'] - df['position_start']).dropna()
        for var, df in list_dfs.items()
    }

    fig, ax = plt.subplots(figsize=(10, 6))
    bp = ax.boxplot(
        [sizes[v] for v in variations],
        labels=variations,
        patch_artist=True,
        notch=True,
        medianprops=dict(color='black', linewidth=2),
    )
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.75)

    ax.set_title('TE Size Distribution per Variation', fontsize=14, fontweight='bold')
    ax.set_xlabel('Variation')
    ax.set_ylabel('TE Length (bp)')
    plt.tight_layout()
    plt.savefig('te_size_distribution.png', dpi=150)
    plt.close()

    print("\n[2] TE size summary (bp):")
    for var, s in sizes.items():
        print(f"  {var}: median={s.median():.0f}, mean={s.mean():.0f}, "
              f"min={s.min():.0f}, max={s.max():.0f}")

# ── 3. Frequency distributions per variation ──────────────────────────────────
def plot_frequency_distribution(list_dfs):
    """Overlapping KDE + histogram of insertion frequencies."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # -- histogram (stacked, semi-transparent)
    ax = axes[0]
    for (var, df), color in zip(list_dfs.items(), colors):
        freq = df['frequency'].dropna()
        ax.hist(freq, bins=20, alpha=0.55, label=var, color=color, edgecolor='white')
    ax.set_title('Frequency Distribution (Histogram)', fontsize=13, fontweight='bold')
    ax.set_xlabel('Insertion Frequency')
    ax.set_ylabel('Count')
    ax.legend()

    # -- violin plot
    ax = axes[1]
    data   = [list_dfs[v]['frequency'].dropna().values for v in variations]
    vparts = ax.violinplot(data, positions=range(len(variations)), showmedians=True)
    for pc, color in zip(vparts['bodies'], colors):
        pc.set_facecolor(color)
        pc.set_alpha(0.75)
    ax.set_xticks(range(len(variations)))
    ax.set_xticklabels(variations)
    ax.set_title('Frequency Distribution (Violin)', fontsize=13, fontweight='bold')
    ax.set_xlabel('Variation')
    ax.set_ylabel('Insertion Frequency')

    plt.tight_layout()
    plt.savefig('te_frequency_distribution.png', dpi=150)
    plt.close()

    print("\n[3] Frequency summary:")
    for var, df in list_dfs.items():
        freq = df['frequency'].dropna()
        print(f"  {var}: mean={freq.mean():.3f}, median={freq.median():.3f}, "
              f"std={freq.std():.3f}, n={len(freq)}")

# ── 4. TE family composition per variation ────────────────────────────────────
def plot_te_families(list_dfs, possible_names):
    """Stacked-bar chart of TE family proportions per variation."""
    family_counts = {}
    for var, df in list_dfs.items():
        df['family'] = df['TE_name'].apply(lambda x: get_te_family(x, possible_names))
        counts = df['family'].value_counts().reindex(possible_names, fill_value=0)
        family_counts[var] = counts

    family_df = pd.DataFrame(family_counts).T   # rows = variations, cols = families

    # -- absolute counts stacked bar
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    cmap = plt.cm.Set2(np.linspace(0, 1, len(possible_names)))

    ax = axes[0]
    bottom = np.zeros(len(variations))
    for fam, color in zip(possible_names, cmap):
        vals = family_df[fam].values
        ax.bar(variations, vals, bottom=bottom, label=fam, color=color, edgecolor='white')
        bottom += vals
    ax.set_title('TE Family Counts per Variation', fontsize=13, fontweight='bold')
    ax.set_xlabel('Variation')
    ax.set_ylabel('Count')
    ax.legend(loc='upper right', fontsize=8)

    # -- proportions stacked bar
    ax = axes[1]
    prop_df = family_df.div(family_df.sum(axis=1), axis=0)
    bottom  = np.zeros(len(variations))
    for fam, color in zip(possible_names, cmap):
        vals = prop_df[fam].values
        ax.bar(variations, vals, bottom=bottom, label=fam, color=color, edgecolor='white')
        bottom += vals
    ax.set_title('TE Family Proportions per Variation', fontsize=13, fontweight='bold')
    ax.set_xlabel('Variation')
    ax.set_ylabel('Proportion')
    ax.legend(loc='upper right', fontsize=8)

    plt.tight_layout()
    plt.savefig('te_family_composition.png', dpi=150)
    plt.close()

    print("\n[4] TE family counts:")
    print(family_df.to_string())

# ── Run all comparisons ───────────────────────────────────────────────────────
if __name__ == '__main__':
    plot_te_count(list_dfs)
    plot_te_size(list_dfs)
    plot_frequency_distribution(list_dfs)
    plot_te_families(list_dfs, possible_names)