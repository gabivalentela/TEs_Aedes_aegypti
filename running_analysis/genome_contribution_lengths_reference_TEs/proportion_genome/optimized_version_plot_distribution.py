import gc
import gzip
from collections import defaultdict

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd


# ── Config ────────────────────────────────────────────────────────────────────

INPUT_FILE = '/work/users/g/a/gabivla/lab/SV_mosquitos/TE_paper_aedesaegytpi/reference_vcf/genotyping_reference_TEs/post_processing/merged_results_5_2_either_filtered.csv'
OUTPUT_DIR = '.'
NROWS      = None   # set to an int (e.g. 5000) for testing, None for full file

CHROM_SIZES = {
    'AaegL5_1': 310_827_022,
    'AaegL5_2': 474_425_716,
    'AaegL5_3': 409_777_670,
}

CLASS_COLORS = {
    'LINE':     '#e63946',
    'SINE':     '#457b9d',
    'LTR':      '#2a9d8f',
    'MITE':     '#e9c46a',
    'TIR':      '#f4a261',
    'HELITRON': '#8338ec',
    'PENELOPE': '#fb5607',
    'UNKNOWN':  '#adb5bd',
    'OTHER':    '#6c757d',
}
SNP_COLOR = '#3d405b'

TE_ORDER = ['LINE', 'SINE', 'LTR', 'MITE', 'TIR', 'HELITRON', 'PENELOPE', 'UNKNOWN', 'OTHER']

# Genotype string → int8 encoding (saves ~10x memory vs keeping strings)
# 0 = ref/ref, 1 = het (0/1 or 1/0), 2 = hom-alt (1/1), -1 = missing
GT_MAP = {'0/0': 0, '0/1': 1, '1/0': 1, '1/1': 2, './.': -1, '.': -1}

# Columns that are metadata, not sample genotypes
BASE_COLUMNS = frozenset({
    'chromosome', 'position_start', 'position_end', 'TE_name', 'strand',
    'TE_class', 'score', 'window', 'length_bp',
})


# ── SNP streaming (no full matrix ever held in memory) ────────────────────────

snp_path = '/users/g/a/gabivla/VCF_SNP_variant/whole_sample_filtered_snps_110722.vcf.gz'

_BASE_VCF_COLS = frozenset({
    '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'
})

# Genotypes counted as "carrying the ALT allele" per site/sample
_ALT_CARRYING_GTS = frozenset({'0/1', '1/0', '1/1'})

# snp_chrom_counts   : chrom -> total SNP count  (all samples, just site counts)
# snp_sample_carrier : sample -> # sites where that sample is 0/1, 1/0, or 1/1
snp_chrom_counts:   dict[str, int] = defaultdict(int)
snp_sample_carrier: dict[str, int] = defaultdict(int)
snp_sample_cols:    list[str]      = []

print('Streaming SNP VCF…')
_count = 0
with gzip.open(snp_path, 'rt') as _f:
    for _line in _f:
        if _line.startswith('##'):
            continue
        if _line.startswith('#CHROM'):
            _all_cols      = _line.strip().split('\t')
            snp_sample_cols = [c for c in _all_cols
                               if c not in _BASE_VCF_COLS and not c.startswith('JB')]
            _col_idx        = {c: i for i, c in enumerate(_all_cols)}
            _sample_indices = [_col_idx[s] for s in snp_sample_cols]
            continue

        _fields = _line.split('\t')
        _chrom  = _fields[0]
        snp_chrom_counts[_chrom] += 1

        for _s, _idx in zip(snp_sample_cols, _sample_indices):
            if _fields[_idx].split(':')[0] in _ALT_CARRYING_GTS:
                snp_sample_carrier[_s] += 1

        _count += 1
        if NROWS and _count >= NROWS:
            break

print(f'  {sum(snp_chrom_counts.values()):,} SNP sites | {len(snp_sample_cols)} samples')


# ── Helpers ───────────────────────────────────────────────────────────────────

def classify_te(name: str) -> str:
    u = str(name).upper()
    for cls in ['LINE', 'SINE', 'LTR', 'MITE', 'TIR', 'HELITRON', 'PENELOPE', 'UNKNOWN']:
        if cls in u:
            return cls
    return 'OTHER'


def bp_coverage(df: pd.DataFrame) -> pd.DataFrame:
    """Sum of (end - start) per (chrom, TE_class)."""
    tmp = df[['chromosome', 'TE_class', 'length_bp']].copy()
    cov = (
        tmp.groupby(['chromosome', 'TE_class'])['length_bp']
        .sum().reset_index()
        .rename(columns={'length_bp': 'coverage_bp'})
    )
    cov['coverage_mb'] = cov['coverage_bp'] / 1_000_000
    return cov


def snp_coverage(snp_chrom_counts: dict, chrom_sizes: dict) -> pd.DataFrame:
    """Build per-chromosome SNP coverage table from the pre-aggregated counts dict."""
    genome_bp = sum(chrom_sizes.values())
    rows = [
        {
            'chromosome':         chrom,
            'snp_count':     cnt,
            'pct_of_chrom':  cnt / chrom_sizes.get(chrom, np.nan) * 100,
            'pct_of_genome': cnt / genome_bp * 100,
        }
        for chrom, cnt in snp_chrom_counts.items()
    ]
    return pd.DataFrame(rows)


def _draw_violin(ax, xi, vals, color, rng):
    """Single violin + jitter strip."""
    if len(vals) == 0:
        return
    if len(vals) < 2:
        ax.scatter([xi], vals, color=color, s=40, zorder=3)
        return
    parts = ax.violinplot([vals], positions=[xi], widths=0.65,
                          showmedians=True, showextrema=True)
    for pc in parts['bodies']:
        pc.set_facecolor(color)
        pc.set_alpha(0.6)
    for part in ('cmedians', 'cmins', 'cmaxes', 'cbars'):
        parts[part].set_color('black')
        parts[part].set_linewidth(0.9)
    jitter = rng.uniform(-0.12, 0.12, size=len(vals))
    ax.scatter(xi + jitter, vals, color=color,
               s=12, alpha=0.5, zorder=3, linewidths=0, edgecolors='none')


# ── TE population-level plots ─────────────────────────────────────────────────

def plot_per_chromosome(coverage: pd.DataFrame, output_dir: str = '.') -> None:
    """Stacked horizontal bar: Mb covered by each TE class per chromosome."""
    chromosomes     = sorted(coverage['chromosome'].unique())
    present_classes = [c for c in TE_ORDER if c in coverage['TE_class'].values]

    pivot = (
        coverage.pivot_table(index='chromosome', columns='TE_class',
                             values='coverage_mb', aggfunc='sum', fill_value=0)
        .reindex(index=chromosomes, columns=present_classes, fill_value=0)
    )

    fig, ax = plt.subplots(figsize=(10, max(4, len(chromosomes) * 0.55 + 1)))
    left = np.zeros(len(chromosomes))
    for cls in present_classes:
        vals = pivot[cls].values
        ax.barh(chromosomes, vals, left=left,
                color=CLASS_COLORS.get(cls, '#999'), label=cls,
                edgecolor='white', linewidth=0.4)
        left += vals

    ax.set_xlabel('Genomic Coverage (Mb)', fontsize=11)
    ax.set_title('TE Coverage per Chromosome', fontsize=13, fontweight='bold')
    ax.legend(title='TE Class', bbox_to_anchor=(1.01, 1), loc='upper left', fontsize=9)
    ax.xaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax.grid(axis='x', linewidth=0.4, alpha=0.5, which='both')
    ax.set_axisbelow(True)
    plt.tight_layout()
    out = f'{output_dir}/te_coverage_per_chromosome.png'
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out}')


def plot_proportion_per_chromosome(coverage: pd.DataFrame, output_dir: str = '.') -> None:
    """Bars normalised to 100 % to show relative composition."""
    chromosomes     = sorted(coverage['chromosome'].unique())
    present_classes = [c for c in TE_ORDER if c in coverage['TE_class'].values]

    pivot = (
        coverage.pivot_table(index='chromosome', columns='TE_class',
                             values='coverage_mb', aggfunc='sum', fill_value=0)
        .reindex(index=chromosomes, columns=present_classes, fill_value=0)
    )
    pivot_pct = pivot.div(pivot.sum(axis=1).replace(0, np.nan), axis=0) * 100

    fig, ax = plt.subplots(figsize=(10, max(4, len(chromosomes) * 0.55 + 1)))
    left = np.zeros(len(chromosomes))
    for cls in present_classes:
        vals = pivot_pct[cls].fillna(0).values
        ax.barh(chromosomes, vals, left=left,
                color=CLASS_COLORS.get(cls, '#999'), label=cls,
                edgecolor='white', linewidth=0.4)
        left += vals

    ax.set_xlabel('Proportion of TE Coverage (%)', fontsize=11)
    ax.set_xlim(0, 100)
    ax.set_title('TE Class Composition per Chromosome (%)', fontsize=13, fontweight='bold')
    ax.legend(title='TE Class', bbox_to_anchor=(1.01, 1), loc='upper left', fontsize=9)
    ax.grid(axis='x', linewidth=0.4, alpha=0.5)
    ax.set_axisbelow(True)
    plt.tight_layout()
    out = f'{output_dir}/te_coverage_proportion_per_chromosome.png'
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out}')


def plot_pct_of_chromosome(coverage: pd.DataFrame,
                           chrom_sizes: dict,
                           output_dir: str = '.') -> None:
    """TE classes as % of each chromosome's length; unannotated fills to 100 %."""
    chromosomes     = sorted(coverage['chromosome'].unique())
    present_classes = [c for c in TE_ORDER if c in coverage['TE_class'].values]

    pivot = (
        coverage.pivot_table(index='chromosome', columns='TE_class',
                             values='coverage_bp', aggfunc='sum', fill_value=0)
        .reindex(index=chromosomes, columns=present_classes, fill_value=0)
    )

    chrom_lengths = pd.Series(chrom_sizes)
    missing = [c for c in chromosomes if c not in chrom_lengths.index]
    if missing:
        print(f'WARNING: no size for {missing} — skipping in this plot.')
        pivot       = pivot.loc[~pivot.index.isin(missing)]
        chromosomes = [c for c in chromosomes if c not in missing]

    pivot_pct   = pivot.div(chrom_lengths[chromosomes], axis=0) * 100
    unannotated = (100 - pivot_pct.sum(axis=1)).clip(lower=0)

    fig, ax = plt.subplots(figsize=(10, max(4, len(chromosomes) * 0.7 + 1.5)))
    left = np.zeros(len(chromosomes))
    for cls in present_classes:
        vals = pivot_pct[cls].values
        ax.barh(chromosomes, vals, left=left,
                color=CLASS_COLORS.get(cls, '#999'), label=cls,
                edgecolor='white', linewidth=0.4)
        left += vals
    ax.barh(chromosomes, unannotated.values, left=left,
            color='#dee2e6', label='Unannotated', edgecolor='white', linewidth=0.4)

    ax.set_xlabel('% of Chromosome Length', fontsize=11)
    ax.set_xlim(0, 100)
    ax.set_title('TE Coverage as % of Chromosome Length', fontsize=13, fontweight='bold')
    ax.legend(title='TE Class', bbox_to_anchor=(1.01, 1), loc='upper left', fontsize=9)
    ax.grid(axis='x', linewidth=0.4, alpha=0.5)
    ax.set_axisbelow(True)

    for y, (chrom, pct) in enumerate(zip(chromosomes, pivot_pct.sum(axis=1))):
        ax.text(pct + 0.5, y, f'{pct:.1f}%', va='center', fontsize=8, color='#333')

    plt.tight_layout()
    out = f'{output_dir}/te_coverage_pct_of_chromosome.png'
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out}')


def plot_genome_wide(coverage: pd.DataFrame,
                     chrom_sizes: dict,
                     output_dir: str = '.') -> None:
    """Two-panel: absolute Mb per TE class (bar) + pie of composition."""
    present_classes = [c for c in TE_ORDER if c in coverage['TE_class'].values]
    totals_bp = coverage.groupby('TE_class')['coverage_bp'].sum().reindex(present_classes, fill_value=0)
    totals_mb = totals_bp / 1_000_000
    genome_bp = sum(chrom_sizes.values())
    colors    = [CLASS_COLORS.get(c, '#999') for c in present_classes]

    fig, (ax_bar, ax_pie) = plt.subplots(
        1, 2, figsize=(14, max(5, len(present_classes) * 0.55 + 2))
    )
    fig.suptitle('Genome-wide TE Coverage', fontsize=14, fontweight='bold')

    bars = ax_bar.barh(present_classes, totals_mb.values,
                       color=colors, edgecolor='white', linewidth=0.4)
    ax_bar.bar_label(bars, fmt='%.2f Mb', padding=4, fontsize=8)
    ax_bar.set_xlabel('Coverage (Mb)', fontsize=11)
    ax_bar.set_title('Absolute Coverage per TE Class', fontsize=11)
    ax_bar.grid(axis='x', linewidth=0.4, alpha=0.5)
    ax_bar.set_axisbelow(True)

    # Summary stacked bar showing annotated vs total genome
    ax_bar.barh(['Total genome'], [genome_bp / 1_000_000],
                color='#dee2e6', edgecolor='white', linewidth=0.4)
    left = 0
    for cls, val, col in zip(present_classes, totals_mb.values, colors):
        ax_bar.barh(['Total genome'], [val], left=left, color=col,
                    edgecolor='white', linewidth=0.4)
        left += val

    wedge_labels = [f'{c}\n{v:.1f} Mb' for c, v in zip(present_classes, totals_mb.values)]
    _, _, autotexts = ax_pie.pie(
        totals_mb.values, labels=wedge_labels, colors=colors,
        autopct='%1.1f%%', startangle=140, pctdistance=0.75,
        wedgeprops=dict(edgecolor='white', linewidth=0.6),
    )
    for t in autotexts:
        t.set_fontsize(8)
    ax_pie.set_title('TE Class Composition\n(annotated bp only)', fontsize=11)

    plt.tight_layout()
    out = f'{output_dir}/te_coverage_genome_wide.png'
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out}')


def save_summary_table(coverage: pd.DataFrame,
                       chrom_sizes: dict,
                       output_dir: str = '.') -> None:
    """Wide-format CSV: rows = chromosomes, columns = TE classes (Mb) + totals."""
    present_classes = [c for c in TE_ORDER if c in coverage['TE_class'].values]
    pivot_mb = (
        coverage.pivot_table(index='chromosome', columns='TE_class',
                             values='coverage_mb', aggfunc='sum', fill_value=0)
        .reindex(columns=present_classes, fill_value=0)
        .round(4)
    )
    genome_mb                   = sum(chrom_sizes.values()) / 1_000_000
    chrom_mb                    = pd.Series({k: v / 1_000_000 for k, v in chrom_sizes.items()})
    pivot_mb['TOTAL_TE_MB']     = pivot_mb.sum(axis=1)
    pivot_mb['CHROM_LENGTH_MB'] = chrom_mb
    pivot_mb['PCT_OF_CHROM']    = (pivot_mb['TOTAL_TE_MB'] / pivot_mb['CHROM_LENGTH_MB'] * 100).round(2)
    pivot_mb['PCT_OF_GENOME']   = (pivot_mb['TOTAL_TE_MB'] / genome_mb * 100).round(2)

    out = f'{output_dir}/te_coverage_summary.csv'
    pivot_mb.to_csv(out)
    print(f'Saved: {out}')
    print('\nCoverage summary (Mb):')
    print(pivot_mb.to_string())


# ── SNP plots ─────────────────────────────────────────────────────────────────

def plot_snp_per_chromosome(snp_cov: pd.DataFrame, output_dir: str = '.') -> None:
    """Horizontal bar: % of each chromosome's length occupied by SNPs."""
    chromosomes = sorted(snp_cov['chromosome'].unique())
    pct = snp_cov.set_index('chromosome').reindex(chromosomes)['pct_of_chrom'].fillna(0)

    fig, ax = plt.subplots(figsize=(9, max(4, len(chromosomes) * 0.7 + 1.5)))
    bars = ax.barh(chromosomes, pct.values,
                   color=SNP_COLOR, edgecolor='white', linewidth=0.4, alpha=0.85)
    ax.bar_label(bars, fmt='%.4f%%', padding=4, fontsize=8)
    ax.set_xlabel('SNPs as % of Chromosome Length', fontsize=11)
    ax.set_title('SNP Density per Chromosome\n(each SNP = 1 bp)',
                 fontsize=13, fontweight='bold')
    ax.grid(axis='x', linewidth=0.4, alpha=0.5)
    ax.set_axisbelow(True)
    plt.tight_layout()
    out = f'{output_dir}/snp_pct_per_chromosome.png'
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out}')


# ── Per-sample genotype-aware TE coverage ────────────────────────────────────

def per_sample_coverage(df: pd.DataFrame, chrom_sizes: dict) -> pd.DataFrame:
    """
    For every sample, sum bp of TEs where the sample carries at least one
    ALT allele (genotype 0/1, 1/0, or 1/1), grouped by (sample, chrom, TE_class).

    Genotype columns are int8: 0=0/0, 1=0/1 or 1/0 (het), 2=1/1 (hom-alt),
    -1=missing. Mask is `>= 1`, so both het and hom-alt calls are counted;
    only ref/ref (0) and missing (-1) are excluded.
    """
    genome_bp   = sum(chrom_sizes.values())
    sample_cols = [c for c in df.columns if c not in BASE_COLUMNS]

    records = []
    for sample in sample_cols:
        mask    = df[sample] >= 1
        sub     = df.loc[mask, ['chromosome', 'TE_class', 'length_bp']]
        grouped = (
            sub.groupby(['chromosome', 'TE_class'])['length_bp']
            .sum().reset_index()
            .rename(columns={'length_bp': 'coverage_bp'})
        )
        grouped['sample']        = sample
        grouped['coverage_mb']   = grouped['coverage_bp'] / 1_000_000
        grouped['pct_of_chrom']  = grouped.apply(
            lambda r: r['coverage_bp'] / chrom_sizes.get(r['chromosome'], np.nan) * 100,
            axis=1,
        )
        grouped['pct_of_genome'] = grouped['coverage_bp'] / genome_bp * 100
        records.append(grouped)

    return pd.concat(records, ignore_index=True)


def plot_variation_boxplot(sample_cov: pd.DataFrame,
                           chrom_sizes: dict,
                           output_dir: str = '.') -> None:
    """Box + jitter per TE class per chromosome, distribution across samples."""
    chromosomes     = sorted(sample_cov['chromosome'].unique())
    present_classes = [c for c in TE_ORDER if c in sample_cov['TE_class'].values]
    n_classes       = len(present_classes)
    x_pos           = np.arange(n_classes)

    ncols = 2
    nrows = (len(chromosomes) + 1) // ncols
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(7 * ncols, 4.5 * nrows),
                             constrained_layout=True)
    axes = np.array(axes).flatten()

    for idx, chrom in enumerate(chromosomes):
        ax  = axes[idx]
        sub = sample_cov[sample_cov['chromosome'] == chrom]

        box_data = [
            sub[sub['TE_class'] == cls]['pct_of_chrom'].values
            if len(sub[sub['TE_class'] == cls]) else np.array([0.0])
            for cls in present_classes
        ]
        bp = ax.boxplot(
            box_data, positions=x_pos, widths=0.5, patch_artist=True,
            medianprops=dict(color='black', linewidth=1.5),
            whiskerprops=dict(linewidth=0.8),
            capprops=dict(linewidth=0.8),
            flierprops=dict(marker='', linestyle='none'),
            zorder=2,
        )
        for patch, cls in zip(bp['boxes'], present_classes):
            patch.set_facecolor(CLASS_COLORS.get(cls, '#999'))
            patch.set_alpha(0.7)

        rng = np.random.default_rng(42)
        for xi, cls in zip(x_pos, present_classes):
            vals = sub[sub['TE_class'] == cls]['pct_of_chrom'].values
            if len(vals):
                jitter = rng.uniform(-0.18, 0.18, size=len(vals))
                ax.scatter(xi + jitter, vals,
                           color=CLASS_COLORS.get(cls, '#999'),
                           s=14, alpha=0.55, zorder=3, linewidths=0)

        ax.set_xticks(x_pos)
        ax.set_xticklabels(present_classes, rotation=35, ha='right', fontsize=8)
        ax.set_ylabel('% of chromosome (per sample)', fontsize=9)
        ax.set_title(f'Chromosome {chrom}', fontsize=11, fontweight='bold')
        ax.grid(axis='y', linewidth=0.4, alpha=0.5)
        ax.set_axisbelow(True)

    for j in range(len(chromosomes), len(axes)):
        axes[j].set_visible(False)

    fig.suptitle('Per-sample TE Coverage Variation\n(0/1, 1/0, and 1/1 genotypes counted)',
                 fontsize=13, fontweight='bold')
    out = f'{output_dir}/te_variation_boxplot_per_chrom.png'
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out}')


def plot_variation_total_genome(sample_cov: pd.DataFrame,
                                chrom_sizes: dict,
                                snp_vals: np.ndarray,
                                output_dir: str = '.') -> None:
    """
    Left  — violin per TE class + SNP violin appended at the right.
    Right — histogram of total TE % of genome per sample.

    snp_vals: pre-computed array of per-sample SNP ALT-carrying count / genome_bp * 100,
              produced from snp_sample_carrier during the streaming pass.
    """
    genome_bp       = sum(chrom_sizes.values())
    present_classes = [c for c in TE_ORDER if c in sample_cov['TE_class'].values]

    total = (
        sample_cov.groupby(['sample', 'TE_class'])['coverage_bp']
        .sum().reset_index()
    )
    total['pct_of_genome'] = total['coverage_bp'] / genome_bp * 100

    total_per_sample = (
        total.groupby('sample')['pct_of_genome'].sum().reset_index()
        .rename(columns={'pct_of_genome': 'total_pct'})
    )

    all_labels = present_classes + ['SNP']

    fig, axes = plt.subplots(1, 2, figsize=(16, 6), constrained_layout=True)
    fig.suptitle(
        'Genome-wide Coverage Variation Across Samples\n(0/1, 1/0, and 1/1 genotypes counted)',
        fontsize=13, fontweight='bold',
    )

    ax  = axes[0]
    rng = np.random.default_rng(42)
    for xi, cls in enumerate(present_classes):
        vals = total[total['TE_class'] == cls]['pct_of_genome'].values
        _draw_violin(ax, xi, vals, CLASS_COLORS.get(cls, '#999'), rng)
    _draw_violin(ax, len(present_classes), snp_vals, SNP_COLOR, rng)

    ax.set_xticks(range(len(all_labels)))
    ax.set_xticklabels(all_labels, rotation=35, ha='right', fontsize=9)
    ax.set_ylabel('% of genome (per sample)', fontsize=10)
    ax.set_title('Distribution per TE Class + SNPs', fontsize=11)
    ax.grid(axis='y', linewidth=0.4, alpha=0.5)
    ax.set_axisbelow(True)

    from matplotlib.patches import Patch
    handles = [Patch(facecolor=CLASS_COLORS.get(c, '#999'), label=c, alpha=0.6)
               for c in present_classes]
    handles.append(Patch(facecolor=SNP_COLOR, label='SNP', alpha=0.6))
    ax.legend(handles=handles, title='Class', fontsize=8,
              bbox_to_anchor=(1.01, 1), loc='upper left')

    ax2  = axes[1]
    vals = total_per_sample['total_pct'].values
    ax2.hist(vals, bins=30, color='#457b9d', edgecolor='white', linewidth=0.4, alpha=0.85)
    ax2.axvline(vals.mean(),     color='#e63946', linewidth=1.5, linestyle='--',
                label=f'Mean   {vals.mean():.2f}%')
    ax2.axvline(np.median(vals), color='#2a9d8f', linewidth=1.5, linestyle=':',
                label=f'Median {np.median(vals):.2f}%')
    ax2.set_xlabel('Total TE % of genome (per sample)', fontsize=10)
    ax2.set_ylabel('Number of samples', fontsize=10)
    ax2.set_title('Total TE Genome Proportion\nDistribution Across Samples', fontsize=11)
    ax2.legend(fontsize=9)
    ax2.grid(axis='y', linewidth=0.4, alpha=0.5)
    ax2.set_axisbelow(True)
    ax2.text(0.97, 0.97,
             f'SD = {vals.std():.3f}%\nRange = {vals.max()-vals.min():.3f}%\nn = {len(vals):,}',
             transform=ax2.transAxes, va='top', ha='right', fontsize=9,
             bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.7))

    out = f'{output_dir}/te_variation_genome_wide.png'
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out}')


def save_variation_table(sample_cov: pd.DataFrame,
                         chrom_sizes: dict,
                         output_dir: str = '.') -> None:
    """
    per_sample_te_coverage.csv — long table
    per_sample_summary.csv     — wide: one row per sample, TE class % of genome
    per_sample_stats.csv       — summary stats across samples
    """
    genome_bp       = sum(chrom_sizes.values())
    present_classes = [c for c in TE_ORDER if c in sample_cov['TE_class'].values]

    out1 = f'{output_dir}/per_sample_te_coverage.csv'
    sample_cov.round(4).to_csv(out1, index=False)
    print(f'Saved: {out1}')

    total = sample_cov.groupby(['sample', 'TE_class'])['coverage_bp'].sum().reset_index()
    total['pct_of_genome'] = (total['coverage_bp'] / genome_bp * 100).round(4)

    wide = (
        total.pivot_table(index='sample', columns='TE_class',
                          values='pct_of_genome', aggfunc='sum', fill_value=0)
        .reindex(columns=present_classes, fill_value=0)
    )
    wide['TOTAL_TE_PCT'] = wide.sum(axis=1).round(4)

    stats = pd.DataFrame({
        cls: {
            'mean':   wide[cls].mean(),
            'median': wide[cls].median(),
            'std':    wide[cls].std(),
            'min':    wide[cls].min(),
            'max':    wide[cls].max(),
        }
        for cls in present_classes + ['TOTAL_TE_PCT']
    }).T.round(4)

    out2 = f'{output_dir}/per_sample_summary.csv'
    wide.to_csv(out2)
    print(f'Saved: {out2}')

    out3 = f'{output_dir}/per_sample_stats.csv'
    stats.to_csv(out3)
    print(f'Saved: {out3}')

    print('\nVariation summary (% of genome):')
    print(stats.to_string())


# ── Main ──────────────────────────────────────────────────────────────────────

print('\nLoading TE file…')
df = pd.read_csv(INPUT_FILE)
print(f'  {len(df):,} rows loaded from {INPUT_FILE}')

# Auto-detect chrom name style and remap CHROM_SIZES if needed
sample_chrom = str(df['chromosome'].iloc[0])
if not sample_chrom.startswith('chr') and any(k.startswith('chr') for k in CHROM_SIZES):
    CHROM_SIZES = {k.replace('chr', ''): v for k, v in CHROM_SIZES.items()}
    print('  Auto-remapped CHROM_SIZES keys')

# ── Prepare TE dataframe (memory optimisation) ────────────────────────────────

# Classify TEs up front so the column is available throughout
df['TE_class']  = df['TE_name'].apply(classify_te)

# Downcast coordinate columns
df['position_start'] = pd.to_numeric(df['position_start'], downcast='unsigned')
df['position_end']   = pd.to_numeric(df['position_end'],   downcast='unsigned')
df['length_bp'] = (df['position_end'] - df['position_start']).astype(np.uint32)

# Encode sample genotype columns as int8 (0/1/2/-1) — ~10x smaller than strings
te_sample_cols = [c for c in df.columns if c not in BASE_COLUMNS
                  and c not in ('TE_name', 'chromosome', 'position_start', 'position_end',
                                'strand', 'score', 'window')]
for col in te_sample_cols:
    df[col] = df[col].map(GT_MAP).fillna(-1).astype(np.int8)

# Drop any columns we won't use again (frees memory before heavy work)
keep = ['chromosome', 'TE_class', 'length_bp'] + te_sample_cols
df   = df[[c for c in keep if c in df.columns]]

# ── TE population-level plots ─────────────────────────────────────────────────

coverage = bp_coverage(df)

plot_per_chromosome(coverage,                          output_dir=OUTPUT_DIR)
plot_proportion_per_chromosome(coverage,               output_dir=OUTPUT_DIR)
plot_pct_of_chromosome(coverage, CHROM_SIZES,         output_dir=OUTPUT_DIR)
plot_genome_wide(coverage,       CHROM_SIZES,         output_dir=OUTPUT_DIR)
save_summary_table(coverage,     CHROM_SIZES,         output_dir=OUTPUT_DIR)

del coverage
gc.collect()

# ── SNP population-level plots ────────────────────────────────────────────────

print('\nBuilding SNP coverage table…')
snp_cov = snp_coverage(snp_chrom_counts, CHROM_SIZES)
plot_snp_per_chromosome(snp_cov, output_dir=OUTPUT_DIR)
del snp_cov
gc.collect()

# ── Per-sample genotype-aware variation ───────────────────────────────────────

print('\nCalculating per-sample TE coverage (0/1, 1/0, and 1/1)…')
sample_cov = per_sample_coverage(df, CHROM_SIZES)

# Pre-compute per-sample SNP % from the streamed counts dict
# (no snp_df matrix — just a small dict lookup)
_genome_bp = sum(CHROM_SIZES.values())
snp_vals   = np.array([
    snp_sample_carrier.get(s, 0) / _genome_bp * 100
    for s in snp_sample_cols
])

plot_variation_boxplot(sample_cov,                     CHROM_SIZES, output_dir=OUTPUT_DIR)
plot_variation_total_genome(sample_cov, CHROM_SIZES,   snp_vals,    output_dir=OUTPUT_DIR)
save_variation_table(sample_cov,        CHROM_SIZES,                output_dir=OUTPUT_DIR)