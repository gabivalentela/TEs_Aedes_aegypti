"""
plot_combined_gene_TE_sweeps_ref_nonref.py
===========================================
One figure per gene, combining BOTH reference (ref) and non-reference (nonref)
TE insertions in the same plot. Reference and non-reference TEs are always
visually distinguishable:

  - In the genomic panel, reference TEs are drawn in a lane ABOVE the gene
    body and non-reference TEs in a lane BELOW it, each lane labeled.
  - Reference TEs are drawn with a solid fill; non-reference TEs are drawn
    with diagonal hatching. This hatch/solid convention is repeated in every
    panel and legend.
  - In the SweepFinder panel, reference-TE regions are shaded yellow and
    non-reference-TE regions are shaded orange.
  - In the frequency panel, reference bars are solid and non-reference bars
    are hatched, grouped by TE and by country.

Three panels per figure:
  1. Genomic view  – gene body, exons, and all overlapping TEs (ref + nonref)
  2. SweepFinder2  – LR curves (all countries except Colombia); countries
                     whose sweep score exceeds the 99th-percentile threshold
                     inside a TE region are highlighted with scatter markers
  3. TE frequencies – proportion of genomes per country that carry each TE

Usage
-----
Adjust the path constants in the "── paths ──" section below, then run:
    python plot_combined_gene_TE_sweeps_ref_nonref.py
"""

import os
import re

import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib import colormaps
import numpy as np
import pandas as pd

# ═══════════════════════════════════════════════════════════════════════════════
# ── paths  (edit these) ────────────────────────────────────────────────────────
# ═══════════════════════════════════════════════════════════════════════════════

ANNOTATION_GFF  = '/work/users/g/a/gabivla/lab/SV_mosquitos/TEs_Aedes_aegypti/Annotation_VectorBase/VectorBase-54_AaegyptiLVP_AGWG.gff'
TE_CSV_REF      = '../getting_TEs_info/TE_information_all_genes_overlap_reference.csv'
TE_CSV_NONREF   = '../getting_TEs_info/TE_information_all_genes_overlap_nonreference.csv'
SWEEPS_DIR      = '../../plot_sweep_findeer/updated_strategy_analysis/make_new_sweep_files/new_sweeps_files/'
OUTPUT_DIR      = './plots_combined_ref_nonref'

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

MARGIN_BP  = 10_000   # flanking region for genomic + sweep panels
GENE_COLOR = 'lightpink'

# Reference vs non-reference styling — reused in every panel/legend
ORIGIN_LABEL = {'ref': 'Reference',     'nonref': 'Non-reference'}
ORIGIN_HATCH = {'ref': '',              'nonref': '///'}
ORIGIN_ALPHA = {'ref': 0.85,            'nonref': 0.55}
TE_REGION_COLOR = {'ref': 'gold',       'nonref': 'darkorange'}


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
    parts = str(te_info_str).split('_')
    try:
        return int(parts[2]), int(parts[3])
    except (IndexError, ValueError):
        return None, None


def get_te_chromosome(te_info_str):
    """AaegL5_3_... → 'AaegL5_3'."""
    parts = str(te_info_str).split('_')
    try:
        return f"{parts[0]}_{parts[1]}"
    except IndexError:
        return None


def short_te_label(te_info_str):
    parts = str(te_info_str).split('_')
    return '_'.join(parts[4:]) if len(parts) > 4 else te_info_str


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
        d = pd.read_csv(os.path.join(sweep_dir, fname), sep='\t')
        threshold_99 = d['LR'].quantile(0.99)
        df_above     = d[d['LR'] >= threshold_99]
        df_region    = d[(d['location'] >= region_start) & (d['location'] <= region_end)]
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


def load_te_csv(path, origin_key):
    """Load a ref or nonref TE csv and tag every row with its origin."""
    d = pd.read_csv(path)
    d['origin'] = origin_key
    d['gene_symbol'] = d['gene_symbol'].str.replace('ID=', '', regex=False)
    d['gene_name_updated'] = d['gene_name'].str.extract(r'Name=([^;]+)')[0]
    return d


# ═══════════════════════════════════════════════════════════════════════════════
# ── load & prepare data ────────────────────────────────────────────────────────
# ═══════════════════════════════════════════════════════════════════════════════

print("Loading annotation…")
gff_full = read_gff(ANNOTATION_GFF)

# Gene coordinates
gene_annot = gff_full[gff_full['type'] == 'protein_coding_gene'].copy()
gene_annot['gene_symbol'] = gene_annot['attributes'].str.extract(r'ID=([^;]+)')[0]

print("Loading TE data (reference + non-reference)…")
df_ref    = load_te_csv(TE_CSV_REF, 'ref')
df_nonref = load_te_csv(TE_CSV_NONREF, 'nonref')
df = pd.concat([df_ref, df_nonref], ignore_index=True)

# Merge gene coordinates
subset_annot = gene_annot[gene_annot['gene_symbol'].isin(df['gene_symbol'])][
    ['gene_symbol', 'start', 'end']
].copy().rename(columns={'start': 'gene_start', 'end': 'gene_end'})

for col in ('start', 'end', 'start_x', 'end_x', 'start_y', 'end_y', 'gene_start', 'gene_end'):
    if col in df.columns:
        df = df.drop(columns=col)

df = df.merge(subset_annot, on='gene_symbol', how='left')
df['display_name'] = df['gene_name_updated'].where(df['gene_name_updated'].notna(), df['gene_symbol'])

# TE coordinates (re-derived from TE_information so ref and nonref are consistent)
df[['te_start', 'te_end']] = df['TE_information'].apply(
    lambda x: pd.Series(parse_te_coords(x))
)

# Country proportion columns — force numeric so bar plots never see stray strings
count_cols = [c for c in df.columns if c.startswith('frequency_')]
for col in count_cols:
    country = col.replace('frequency_', '')
    if country in TOTAL_GENOMES:
        df[f'Proportion_{country}'] = pd.to_numeric(df[col], errors='coerce').fillna(0)

print(f"Found {df['gene_symbol'].nunique()} unique genes "
      f"({(df['origin'] == 'ref').sum()} ref rows, "
      f"{(df['origin'] == 'nonref').sum()} nonref rows).\n")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ═══════════════════════════════════════════════════════════════════════════════
# ── main loop: one figure per gene (ref + nonref combined) ───────────────────
# ═══════════════════════════════════════════════════════════════════════════════

for gene_symbol, gene_df in df.groupby('gene_symbol'):

    gene_start   = gene_df['gene_start'].iloc[0]
    gene_end     = gene_df['gene_end'].iloc[0]
    display_name = gene_df['display_name'].iloc[0]

    tes = gene_df[gene_df['te_start'].notna()].copy()
    if tes.empty:
        print(f"  Skipping {display_name}: no TE coordinates.")
        continue

    n_ref    = int((tes['origin'] == 'ref').sum())
    n_nonref = int((tes['origin'] == 'nonref').sum())

    # ── decide chromosome from first TE ───────────────────────────────────────
    chrom = get_te_chromosome(tes['TE_information'].iloc[0])

    # ── plot region ───────────────────────────────────────────────────────────
    region_start = min(gene_start, tes['te_start'].min()) - MARGIN_BP
    region_end   = max(gene_end,   tes['te_end'].max())   + MARGIN_BP

    # ── sweep data ────────────────────────────────────────────────────────────
    sweep_data = load_sweep_data(SWEEPS_DIR, chrom, region_start, region_end)

    # Determine which countries have a sweep above threshold inside ANY TE
    # (kept per-origin so the sweep title can report ref vs nonref separately)
    countries_with_sweep = {'ref': set(), 'nonref': set()}
    for pop, data in sweep_data.items():
        for _, te_row in tes.iterrows():
            te_s, te_e = int(te_row['te_start']), int(te_row['te_end'])
            hits = data['threshold_region'][
                (data['threshold_region']['location'] >= te_s) &
                (data['threshold_region']['location'] <= te_e)
            ]
            if not hits.empty:
                countries_with_sweep[te_row['origin']].add(pop)

    has_sweeps = bool(countries_with_sweep['ref'] or countries_with_sweep['nonref'])

    # ── TE colour map (consistent across panels; ref & nonref TEs get their
    #    own colour since TE_information encodes distinct coordinates) ────────
    unique_te_names = tes['TE_information'].unique()
    cmap_te = colormaps['tab20'].resampled(max(len(unique_te_names), 1))
    te_color_map = {name: cmap_te(i) for i, name in enumerate(unique_te_names)}

    # ── figure layout ─────────────────────────────────────────────────────────
    n_rows     = 3 if has_sweeps else 2
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
        f"Chrom: {chrom}   |   {n_ref} reference TE(s)   |   {n_nonref} non-reference TE(s)",
        fontsize=16, fontweight='bold', y=1.01
    )

    # ══════════════════════════════════════════════════════════════════════════
    # PANEL 1 – Genomic view
    #   Layout (bottom → top): non-ref TE lane / gene+exons / ref TE lane
    # ══════════════════════════════════════════════════════════════════════════
    ax = ax_genomic

    Y_NONREF_TE = 0.15
    Y_GENE      = 0.50
    Y_REF_TE    = 0.85
    LANE_H      = 0.16

    # background band
    ax.add_patch(patches.Rectangle(
        (region_start, 0), region_end - region_start, 1.0,
        linewidth=0, facecolor='lightgray', alpha=0.20
    ))

    # gene body
    ax.add_patch(patches.Rectangle(
        (gene_start, Y_GENE - 0.09), gene_end - gene_start, 0.18,
        linewidth=2, edgecolor='black', facecolor=GENE_COLOR, alpha=0.85, zorder=3
    ))
    # exons
    exon_list = get_exons_for_gene(gff_full, gene_symbol)
    for exon in exon_list:
        ax.add_patch(patches.Rectangle(
            (exon['start_position'], Y_GENE - 0.045),
            exon['end_position'] - exon['start_position'], 0.09,
            linewidth=1, edgecolor='black', facecolor='black', zorder=4
        ))

    # TEs — ref lane above the gene, non-ref lane below it
    for _, te_row in tes.iterrows():
        te_name  = te_row['TE_information']
        origin   = te_row['origin']
        te_color = te_color_map[te_name]
        te_s     = int(te_row['te_start'])
        te_e     = int(te_row['te_end'])
        y_lane   = Y_REF_TE if origin == 'ref' else Y_NONREF_TE

        ax.add_patch(patches.Rectangle(
            (te_s, y_lane - LANE_H / 2), te_e - te_s, LANE_H,
            linewidth=2, edgecolor='black', facecolor=te_color,
            hatch=ORIGIN_HATCH[origin], alpha=ORIGIN_ALPHA[origin], zorder=2
        ))

    # lane labels
    ax.text(region_start, Y_REF_TE, ' Reference TEs', va='center', ha='left',
            fontsize=9, fontweight='bold', color='dimgray',
            transform=ax.transData, zorder=5)
    ax.text(region_start, Y_NONREF_TE, ' Non-reference TEs', va='center', ha='left',
            fontsize=9, fontweight='bold', color='dimgray',
            transform=ax.transData, zorder=5)

    # legend: one patch per TE colour, plus origin key (hatch vs solid)
    legend_handles = [
        patches.Patch(color=GENE_COLOR, label=f'Gene ({display_name})'),
        patches.Patch(color='black',    label='Exons'),
        patches.Patch(facecolor='white', edgecolor='black', hatch='',
                      label='Reference TE (solid)'),
        patches.Patch(facecolor='white', edgecolor='black', hatch='///',
                      label='Non-reference TE (hatched)'),
    ]
    for te_name, tc in te_color_map.items():
        legend_handles.append(patches.Patch(color=tc, label=f'TE: {short_te_label(te_name)}'))

    ax.legend(handles=legend_handles, loc='lower right', fontsize=7,
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
            if not data['threshold_region'].empty:
                ax.scatter(
                    data['threshold_region']['location'],
                    data['threshold_region']['LR'],
                    color=color, s=60, zorder=5,
                    label=f'{pop_code} ≥99th pct'
                )

        # shade TE regions — colour encodes ref vs nonref
        ref_labeled = False
        nonref_labeled = False
        for _, te_row in tes.iterrows():
            origin = te_row['origin']
            lbl = None
            if origin == 'ref' and not ref_labeled:
                lbl, ref_labeled = 'Reference TE region', True
            elif origin == 'nonref' and not nonref_labeled:
                lbl, nonref_labeled = 'Non-reference TE region', True
            ax.axvspan(
                int(te_row['te_start']), int(te_row['te_end']),
                color=TE_REGION_COLOR[origin], alpha=0.25,
                label=lbl or ''
            )

        # shade exon regions
        for i, exon in enumerate(exon_list):
            ax.axvspan(
                exon['start_position'], exon['end_position'],
                color='pink', alpha=0.20,
                label='Exon region' if i == 0 else ''
            )

        ax.set_xlim(region_start, region_end)
        ax.set_xlabel('Genomic position (bp)', fontsize=11)
        ax.set_ylabel('LR', fontsize=11)

        sweep_summary_parts = []
        if countries_with_sweep['ref']:
            sweep_summary_parts.append(f"ref: {', '.join(sorted(countries_with_sweep['ref']))}")
        if countries_with_sweep['nonref']:
            sweep_summary_parts.append(f"nonref: {', '.join(sorted(countries_with_sweep['nonref']))}")
        ax.set_title(
            f"SweepFinder2 LR  |  Countries with sweep inside TE ({'; '.join(sweep_summary_parts)})",
            fontsize=12, loc='left', pad=4
        )
        ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left', fontsize=9)
        ax.ticklabel_format(style='plain', axis='x')
        plt.setp(ax.get_xticklabels(), rotation=30, ha='right', fontsize=9)
        ax.grid(True, alpha=0.25)

    # ══════════════════════════════════════════════════════════════════════════
    # PANEL 3 – TE frequencies per country (ref = solid, nonref = hatched)
    # ══════════════════════════════════════════════════════════════════════════
    ax = ax_freq

    countries_ordered = list(TOTAL_GENOMES.keys())
    unique_te_full = list(tes[['TE_information', 'origin']].drop_duplicates().itertuples(index=False, name=None))
    n_bars    = len(unique_te_full)
    bar_width = 0.8 / max(n_bars, 1)
    x_base    = np.arange(len(countries_ordered))

    seen_origin_in_legend = set()
    for i, (te_name, origin) in enumerate(unique_te_full):
        te_row   = tes[(tes['TE_information'] == te_name) & (tes['origin'] == origin)].iloc[0]
        te_color = te_color_map[te_name]

        proportions = [
            float(te_row.get(f'Proportion_{c}', 0) or 0)
            for c in countries_ordered
        ]
        offset = (i - (n_bars - 1) / 2) * bar_width

        bar_label = f'TE: {short_te_label(te_name)} ({ORIGIN_LABEL[origin]})'
        ax.bar(
            x_base + offset, proportions,
            width=bar_width * 0.9,
            color=te_color, edgecolor='black', linewidth=0.8,
            hatch=ORIGIN_HATCH[origin], alpha=ORIGIN_ALPHA[origin],
            label=bar_label
        )

    ax.set_xticks(x_base)
    ax.set_xticklabels(countries_ordered, rotation=30, ha='right', fontsize=11)
    ax.set_ylabel('Proportion of genomes with TE', fontsize=11)
    ax.set_ylim(0, 1.05)
    ax.set_title('TE frequency per country (solid = reference, hatched = non-reference)',
                 fontsize=12, loc='left', pad=4)
    ax.legend(loc='upper right', fontsize=7, ncol=min(2, max(n_bars, 1)))
    ax.grid(True, alpha=0.25, axis='y')

    # ── colour the x-tick labels by country ──────────────────────────────────
    for label, country in zip(ax.get_xticklabels(), countries_ordered):
        label.set_color(COUNTRY_COLORS.get(country, 'black'))
        label.set_fontweight('bold')

    # ── save ──────────────────────────────────────────────────────────────────
    plt.tight_layout()
    out_path = os.path.join(OUTPUT_DIR, f'{display_name}_{gene_symbol}_combined_ref_nonref.svg')
    plt.savefig(out_path, format='svg', bbox_inches='tight')
    plt.close()
    print(f"Saved: {out_path}  "
          f"({'sweep panel included' if has_sweeps else 'no sweeps above threshold – sweep panel omitted'})")

print("\nAll done.")