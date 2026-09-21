"""
plot_combined_gene_TE_sweeps.py
================================
One figure per gene containing three sections:
  1. Genomic view  – gene body, exons, and all overlapping non-ref TEs
  2. SweepFinder2  – LR curves (all countries except Colombia); countries whose
                     sweep score exceeds the 99th-percentile threshold inside
                     the TE region are highlighted with scatter markers
  3. TE frequencies – proportion of genomes per country that carry each TE

Usage
-----
Adjust the path constants in the "── paths ──" section below, then run:
    python plot_combined_gene_TE_sweeps.py
"""

import os
import re

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pandas as pd

# ═══════════════════════════════════════════════════════════════════════════════
# ── paths  (edit these) ────────────────────────────────────────────────────────
# ═══════════════════════════════════════════════════════════════════════════════

ANNOTATION_GFF   = '../../../../../../../../../../SV_mosquitos/TEs_Aedes_aegypti/Annotation_VectorBase/VectorBase-54_AaegyptiLVP_AGWG.gff'
TE_CSV           = '../getting_TEs_info/TE_information_all_genes_overlap_nonreference.csv'
SWEEPS_DIR       = '../../plot_sweep_findeer/updated_strategy_analysis/make_new_sweep_files/new_sweeps_files/'
OUTPUT_DIR       = './plots_combined_nonref'

# ═══════════════════════════════════════════════════════════════════════════════
# ── constants ──────────────────────────────────────────────────────────────────
# ═══════════════════════════════════════════════════════════════════════════════

TOTAL_GENOMES = {
    'Brazil':   18,
    'Colombia': 24,
    'USA':      28,
    'Senegal':  20,
    'Gabon':    13,
    'Kenya':    19,
}

# Consistent country colours across all panels
COUNTRY_COLORS = {
    'Brazil':   '#60D08F',
    'USA':      '#785EF0',
    'Colombia': '#DC267F',
    'Senegal':  '#FE6100',
    'Gabon':    '#FFB000',
    'Kenya':    '#648FFF',
}

MARGIN_BP    = 10_000   # flanking region for genomic + sweep panels
GENE_COLOR   = 'lightpink'


# ═══════════════════════════════════════════════════════════════════════════════
# ── helpers ────────────────────────────────────────────────────────────────────
# ═══════════════════════════════════════════════════════════════════════════════

def read_gff(path):
    df = pd.read_csv(path, sep='\t', comment='#', header=None)
    df.columns = ['seqid', 'source', 'type', 'start', 'end',
                  'score', 'strand', 'phase', 'attributes']
    return df


def get_exons_for_gene(gff_df, gene_code):
    """Return list of dicts {exon_name, start_position, end_position}."""
    exon_df = gff_df[gff_df['type'] == 'exon']
    gene_exons = exon_df[exon_df['attributes'].str.contains(gene_code, na=False)].copy()
    gene_exons['exon_id'] = gene_exons['attributes'].str.extract(r'ID=([^;]+)')[0]
    return [
        {'exon_name': row['exon_id'],
         'start_position': row['start'],
         'end_position':   row['end']}
        for _, row in gene_exons.iterrows()
    ]


def parse_te_coords(te_info_str):
    """AaegL5_<chr>_<start>_<end>_TE_... → (start, end)."""
    parts = te_info_str.split('_')
    try:
        return int(parts[2]), int(parts[3])
    except (IndexError, ValueError):
        return None, None


def get_te_chromosome(te_info_str):
    """AaegL5_3_... → 'AaegL5_3'."""
    parts = te_info_str.split('_')
    try:
        return f"{parts[0]}_{parts[1]}"
    except IndexError:
        return None


def load_sweep_data(sweep_dir, chrom, region_start, region_end):
    """
    Load sweep LR files for *chrom* (excluding Colombia).
    Returns dict[pop_code] = {'full_region': df, 'threshold_region': df}
    """
    sweep_files = [
        f for f in os.listdir(sweep_dir)
        if chrom in f and 'Colombia' not in f
    ]
    result = {}
    for fname in sweep_files:
        match = re.search(r'AaegL5_\d_no_relatives_(.+)_raw.+', fname)
        if not match:
            continue
        pop_code = match.group(1)
        df = pd.read_csv(os.path.join(sweep_dir, fname), sep='\t')
        threshold_99 = df['LR'].quantile(0.99)
        df_above     = df[df['LR'] >= threshold_99]
        df_region    = df[(df['location'] >= region_start) & (df['location'] <= region_end)]
        df_threshold_region = df_above[
            (df_above['location'] >= region_start) &
            (df_above['location'] <= region_end)
        ]
        result[pop_code] = {
            'full_region':      df_region,
            'threshold_region': df_threshold_region,
            'threshold_99':     threshold_99,
        }
    return result


# ═══════════════════════════════════════════════════════════════════════════════
# ── load & prepare data ────────────────────────────────────────────────────────
# ═══════════════════════════════════════════════════════════════════════════════

print("Loading annotation…")
gff_full = read_gff(ANNOTATION_GFF)

# Gene coordinates
gene_annot = gff_full[gff_full['type'] == 'protein_coding_gene'].copy()
gene_annot['gene_symbol'] = gene_annot['attributes'].str.extract(r'ID=([^;]+)')[0]

print("Loading TE data…")
df = pd.read_csv(TE_CSV)
df['gene_symbol']       = df['gene_symbol'].str.replace('ID=', '', regex=False)
df['gene_name_updated'] = df['gene_name'].str.extract(r'Name=([^;]+)')[0]

# Merge gene coordinates
subset_annot = gene_annot[gene_annot['gene_symbol'].isin(df['gene_symbol'])][
    ['gene_symbol', 'start', 'end']
].copy().rename(columns={'start': 'gene_start', 'end': 'gene_end'})

for col in ('start', 'end', 'start_x', 'end_x', 'start_y', 'end_y', 'gene_start', 'gene_end'):
    if col in df.columns:
        df = df.drop(columns=col)

df = df.merge(subset_annot, on='gene_symbol', how='left')
df['display_name'] = df['gene_name_updated'].where(df['gene_name_updated'].notna(), df['gene_symbol'])

# TE coordinates
df[['te_start', 'te_end']] = df['TE_information'].apply(
    lambda x: pd.Series(parse_te_coords(x))
)

# Country proportion columns
count_cols = [c for c in df.columns if c.startswith('frequency_')]
for col in count_cols:
    country = col.replace('frequency_', '')
    if country in TOTAL_GENOMES:
        df[f'Proportion_{country}'] = pd.to_numeric(df[col], errors='coerce').fillna(0)

print(f"Found {df['gene_symbol'].nunique()} unique genes.\n")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ═══════════════════════════════════════════════════════════════════════════════
# ── main loop: one figure per gene ────────────────────────────────────────────
# ═══════════════════════════════════════════════════════════════════════════════

for gene_symbol, gene_df in df.groupby('gene_symbol'):

    gene_start   = gene_df['gene_start'].iloc[0]
    gene_end     = gene_df['gene_end'].iloc[0]
    display_name = gene_df['display_name'].iloc[0]

    tes = gene_df[gene_df['te_start'].notna()].copy()
    if tes.empty:
        print(f"  Skipping {display_name}: no TE coordinates.")
        continue

    # ── decide chromosome from first TE ───────────────────────────────────────
    chrom = get_te_chromosome(tes['TE_information'].iloc[0])

    # ── plot region ───────────────────────────────────────────────────────────
    region_start = min(gene_start, tes['te_start'].min()) - MARGIN_BP
    region_end   = max(gene_end,   tes['te_end'].max())   + MARGIN_BP

    # ── sweep data ────────────────────────────────────────────────────────────
    sweep_data = load_sweep_data(SWEEPS_DIR, chrom, region_start, region_end)

    # Determine which countries have a sweep above threshold inside ANY TE
    countries_with_sweep = []
    for pop, data in sweep_data.items():
        for _, te_row in tes.iterrows():
            te_s, te_e = int(te_row['te_start']), int(te_row['te_end'])
            hits = data['threshold_region'][
                (data['threshold_region']['location'] >= te_s) &
                (data['threshold_region']['location'] <= te_e)
            ]
            if not hits.empty:
                countries_with_sweep.append(pop)
                break

    has_sweeps = bool(countries_with_sweep)

    # ── TE colour map (consistent across panels) ──────────────────────────────
    unique_te_names = tes['TE_information'].unique()
    cmap_te         = cm.get_cmap('tab10', len(unique_te_names))
    te_color_map    = {name: cmap_te(i) for i, name in enumerate(unique_te_names)}

    # ── figure layout ─────────────────────────────────────────────────────────
    # Rows: 1 genomic + (optional) 1 sweep + 1 frequency
    n_rows    = 3 if has_sweeps else 2
    row_ratios = [2, 3, 2] if has_sweeps else [2, 2]

    fig, axes = plt.subplots(
        n_rows, 1,
        figsize=(18, 5 * n_rows),
        gridspec_kw={'height_ratios': row_ratios},
    )
    if n_rows == 2:
        ax_genomic, ax_freq = axes
        ax_sweep = None
    else:
        ax_genomic, ax_sweep, ax_freq = axes

    fig.suptitle(
        f"Gene: {display_name}  ({gene_symbol})\n"
        f"Chrom: {chrom}   |   {len(tes)} non-ref TE insertion(s)",
        fontsize=16, fontweight='bold', y=1.01
    )

    # ══════════════════════════════════════════════════════════════════════════
    # PANEL 1 – Genomic view
    # ══════════════════════════════════════════════════════════════════════════
    ax = ax_genomic
    y  = 0.5

    # background band
    ax.add_patch(patches.Rectangle(
        (region_start, y - 0.35), region_end - region_start, 0.70,
        linewidth=0, facecolor='lightgray', alpha=0.25
    ))
    # gene body
    ax.add_patch(patches.Rectangle(
        (gene_start, y - 0.20), gene_end - gene_start, 0.40,
        linewidth=2, edgecolor='black', facecolor=GENE_COLOR, alpha=0.85
    ))
    # exons
    exon_list = get_exons_for_gene(gff_full, gene_symbol)
    for exon in exon_list:
        ax.add_patch(patches.Rectangle(
            (exon['start_position'], y - 0.10),
            exon['end_position'] - exon['start_position'], 0.20,
            linewidth=1, edgecolor='black', facecolor='black'
        ))
    # TEs
    for _, te_row in tes.iterrows():
        te_name  = te_row['TE_information']
        te_color = te_color_map[te_name]
        te_s     = int(te_row['te_start'])
        te_e     = int(te_row['te_end'])
        ax.add_patch(patches.Rectangle(
            (te_s, y - 0.15), te_e - te_s, 0.30,
            linewidth=2, edgecolor=te_color, facecolor=te_color, alpha=0.85
        ))

    # legend
    legend_handles = [
        patches.Patch(color=GENE_COLOR, label=f'Gene ({display_name})'),
        patches.Patch(color='black',    label='Exons'),
    ]
    for te_name, tc in te_color_map.items():
        short = '_'.join(te_name.split('_')[4:]) if '_' in te_name else te_name
        legend_handles.append(patches.Patch(color=tc, label=f'TE: {short}'))

    ax.legend(handles=legend_handles, loc='lower right', fontsize=8,
              ncol=min(3, len(legend_handles)))
    ax.set_xlim(region_start, region_end)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_xlabel('Genomic position (bp)', fontsize=11)
    ax.set_title('Genomic view', fontsize=12, loc='left', pad=4)
    ax.ticklabel_format(style='plain', axis='x')
    plt.setp(ax.get_xticklabels(), rotation=30, ha='right', fontsize=9)
    ax.grid(True, alpha=0.3, axis='x')

    # ══════════════════════════════════════════════════════════════════════════
    # PANEL 2 – SweepFinder LR  (only drawn when there are sweeps)
    # ══════════════════════════════════════════════════════════════════════════
    if has_sweeps and ax_sweep is not None:
        ax = ax_sweep

        for pop_code, data in sweep_data.items():
            color = COUNTRY_COLORS.get(pop_code, 'gray')
            if data['full_region'].empty:
                continue
            ax.plot(
                data['full_region']['location'],
                data['full_region']['LR'],
                label=pop_code, color=color, alpha=0.7, linewidth=1.5
            )
            # scatter for above-threshold within the whole region
            if not data['threshold_region'].empty:
                ax.scatter(
                    data['threshold_region']['location'],
                    data['threshold_region']['LR'],
                    color=color, s=60, zorder=5,
                    label=f'{pop_code} ≥99th pct'
                )

        # shade TE regions
        for i, (_, te_row) in enumerate(tes.iterrows()):
            ax.axvspan(
                int(te_row['te_start']), int(te_row['te_end']),
                color='yellow', alpha=0.25,
                label='TE region' if i == 0 else ''
            )

        # shade exon regions
        for i, exon in enumerate(exon_list):
            ax.axvspan(
                exon['start_position'], exon['end_position'],
                color='pink', alpha=0.25,
                label='Exon region' if i == 0 else ''
            )

        ax.set_xlim(region_start, region_end)
        ax.set_xlabel('Genomic position (bp)', fontsize=11)
        ax.set_ylabel('LR', fontsize=11)
        ax.set_title(
            f'SweepFinder2 LR  |  Countries with sweep inside TE: '
            f'{", ".join(sorted(set(countries_with_sweep)))}',
            fontsize=12, loc='left', pad=4
        )
        ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left', fontsize=9)
        ax.ticklabel_format(style='plain', axis='x')
        plt.setp(ax.get_xticklabels(), rotation=30, ha='right', fontsize=9)
        ax.grid(True, alpha=0.25)

    # ══════════════════════════════════════════════════════════════════════════
    # PANEL 3 – TE frequencies per country
    # ══════════════════════════════════════════════════════════════════════════
    ax = ax_freq

    countries_ordered = list(TOTAL_GENOMES.keys())
    n_tes     = len(unique_te_names)
    bar_width = 0.8 / n_tes
    x_base    = np.arange(len(countries_ordered))

    for i, te_name in enumerate(unique_te_names):
        te_row   = tes[tes['TE_information'] == te_name].iloc[0]
        te_color = te_color_map[te_name]
        short    = '_'.join(te_name.split('_')[4:]) if '_' in te_name else te_name

        proportions = [
            te_row.get(f'Proportion_{c}', 0) or 0
            for c in countries_ordered
        ]
        offset = (i - (n_tes - 1) / 2) * bar_width
        ax.bar(
            x_base + offset, proportions,
            width=bar_width * 0.9,
            color=te_color, edgecolor='black', linewidth=0.8,
            label=f'TE: {short}'
        )

    ax.set_xticks(x_base)
    ax.set_xticklabels(countries_ordered, rotation=30, ha='right', fontsize=11)
    ax.set_ylabel('Proportion of genomes with TE', fontsize=11)
    ax.set_ylim(0, 1.05)
    ax.set_title('TE frequency per country', fontsize=12, loc='left', pad=4)
    ax.legend(loc='upper right', fontsize=8,
              ncol=min(3, n_tes))
    ax.grid(True, alpha=0.25, axis='y')

    # ── colour the x-tick labels by country ──────────────────────────────────
    for label, country in zip(ax.get_xticklabels(), countries_ordered):
        label.set_color(COUNTRY_COLORS.get(country, 'black'))
        label.set_fontweight('bold')

    # ── save ──────────────────────────────────────────────────────────────────
    plt.tight_layout()
    out_path = os.path.join(OUTPUT_DIR, f'{display_name}_{gene_symbol}_combined.svg')
    plt.savefig(out_path, format='svg', bbox_inches='tight')
    plt.close()
    print(f"Saved: {out_path}  "
          f"({'sweep panel included' if has_sweeps else 'no sweeps above threshold – sweep panel omitted'})")

print("\nAll done.")