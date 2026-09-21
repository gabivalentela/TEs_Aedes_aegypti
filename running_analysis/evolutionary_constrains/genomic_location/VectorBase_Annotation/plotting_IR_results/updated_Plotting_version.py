import pandas as pd
import numpy as np
import re
import os

def load_te_csv(path, origin_key):
    """Load a ref or nonref TE csv and tag every row with its origin."""
    d = pd.read_csv(path)
    d['origin'] = origin_key
    d['gene_symbol'] = d['gene_symbol'].str.replace('ID=', '', regex=False)
    # VGSC uses its VectorBase gene code; all other genes use the Name field.
    d['gene_name_updated'] = d['gene_name'].str.extract(r'Name=([^;]+)')[0]
    d.loc[d['gene_symbol'].eq('AAEL023266'), 'gene_name_updated'] = 'VGSC'
    return d

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

def get_te_chromosome(te_info_str):
    """AaegL5_3_... → 'AaegL5_3'."""
    parts = str(te_info_str).split('_')
    try:
        return f"{parts[0]}_{parts[1]}"
    except IndexError:
        return None

def parse_te_coords(te_info_str):
    """AaegL5_<chr>_<start>_<end>_TE_... → (start, end)."""
    parts = str(te_info_str).split('_')
    try:
        return int(parts[2]), int(parts[3])
    except (IndexError, ValueError):
        return None, None

import matplotlib.pyplot as plt
def get_gene_strand(gff_df, gene_code):
    """Return the strand ('+' or '-') for a protein_coding_gene."""
    gene_rows = gff_df[
        (gff_df['type'] == 'protein_coding_gene') &
        (gff_df['attributes'].str.contains(gene_code, na=False))
    ]
    if gene_rows.empty:
        return None
    return gene_rows['strand'].iloc[0]


def get_gene_chromosome(gff_df, gene_code):
    """Return the seqid (chromosome) for a protein_coding_gene."""
    gene_rows = gff_df[
        (gff_df['type'] == 'protein_coding_gene') &
        (gff_df['attributes'].str.contains(gene_code, na=False))
    ]
    if gene_rows.empty:
        return None
    return gene_rows['seqid'].iloc[0]


def parse_te_family(te_info_str):
    """
    Extract the TE family/type from a TE_information string of the form
    AaegL5_<chr>_<start>_<end>_TE_<family...>. Returns everything after
    the '_TE_' marker. Adjust this if your TE_information format differs
    from that pattern.
    """
    s = str(te_info_str)
    marker = '_TE_'
    idx = s.find(marker)
    if idx == -1:
        return 'TE'
    return s[idx + len(marker):]


def build_te_color_map(significant_sweep_df):
    """Assign each unique TE family a fixed color, shared across all genes."""
    families = sorted(significant_sweep_df['TE_information'].apply(parse_te_family).unique())
    cmap = plt.get_cmap('tab20')
    return {fam: cmap(i % 20) for i, fam in enumerate(families)}


from matplotlib.ticker import FuncFormatter

TE_CLASS_COLORS = {
    'SINE':     '#56B4E9',  # sky blue
    'MITE':     '#D55E00',  # vermillion - bold, draws the eye
    'TIR':      '#009E73',  # bluish green
    'LINE':     '#0072B2',  # blue
    'Unknown':  '#FF0000',  # red
    'Helitron': '#E69F00',  # orange
    'LTR':      '#CC79A7',  # reddish purple - bold, draws the eye
    'other':    '#666666',  # dark gray - muted, low priority
}

def get_te_name(te_row):
    """Prefer an explicit TE_name column if present; otherwise fall back
    to parsing it out of TE_information."""
    if 'TE_name' in te_row.index and pd.notna(te_row.get('TE_name')):
        return te_row['TE_name']
    return parse_te_family(te_row['TE_information'])

def get_te_main_class(te_name):
    """
    Collapse a full TE name down to one of the broad classes in
    TE_CLASS_COLORS, so the legend doesn't need the full name.
    """
    main_classes = ['SINE', 'MITE', 'TIR', 'LINE', 'Unknown', 'Helitron', 'LTR']
    if pd.isna(te_name):
        return 'other'
    name_parts = str(te_name).split('_')
    for main_class in main_classes:
        if main_class in name_parts:
            return main_class
    return 'other'


def compute_te_lanes(te_rows, pad_frac=0.01):
    """
    Greedy interval packing so TEs that overlap (or sit very close
    together) get stacked into separate lanes instead of drawn on top of
    each other - this is what was making the VGSC panel unreadable.

    te_rows: list of (te_start, te_end) tuples, in the same order as the
    rows they came from.
    pad_frac: extra horizontal buffer (as a fraction of the full TE span)
    so TEs that are merely close together, not just overlapping, still
    get separated.

    Returns (lanes, n_lanes) where lanes is a list of lane indices
    parallel to te_rows.
    """
    if not te_rows:
        return [], 0
    span = max(e for _, e in te_rows) - min(s for s, _ in te_rows)
    pad = span * pad_frac if span > 0 else 1

    order = sorted(range(len(te_rows)), key=lambda i: te_rows[i][0])
    lanes_end = []
    lanes = [0] * len(te_rows)
    for i in order:
        s, e = te_rows[i]
        placed = False
        for lane_i, lane_end in enumerate(lanes_end):
            if s > lane_end + pad:
                lanes[i] = lane_i
                lanes_end[lane_i] = e
                placed = True
                break
        if not placed:
            lanes[i] = len(lanes_end)
            lanes_end.append(e)
    return lanes, len(lanes_end)

def draw_te_frequency_bars(ax, gene_df):
    """
    Left-side bar plot: for every TE belonging to this gene, show its
    frequency in each country. If the gene has more than one TE, each
    TE gets its own color (by TE class, matching the gene panel) and
    edge style (solid = ref, dashed = nonref), offset side-by-side
    within each country group, and numbered (1, 2, ...) above its
    bars - matching the numbers drawn under each TE box in the gene
    panel next to it.
    """
    freq_cols = [c for c in gene_df.columns if c.startswith('frequency_') and c != 'frequency_diff']
    countries = [c.replace('frequency_', '') for c in freq_cols]
    te_rows = list(gene_df.iterrows())
    n_te = len(te_rows)

    x = np.arange(len(countries))
    total_width = 0.8
    bar_width = total_width / max(n_te, 1)

    for te_i, (_, te_row) in enumerate(te_rows):
        te_class = get_te_main_class(get_te_name(te_row))
        color = TE_CLASS_COLORS.get(te_class, TE_CLASS_COLORS['other'])
        linestyle = 'dashed' if te_row['origin'] == 'nonref' else 'solid'
        values = [te_row[c] if pd.notna(te_row[c]) else 0 for c in freq_cols]
        offset = (te_i - (n_te - 1) / 2) * bar_width

        ax.bar(x + offset, values, width=bar_width * 0.9,
               facecolor=color, edgecolor='black', linestyle=linestyle, linewidth=1.2)

        if n_te > 1:
            max_val = max(values) if values else 0
            ax.text(x[0] + offset, max_val + 0.03, str(te_i + 1),
                     ha='center', va='bottom', fontsize=12, fontweight='bold')

    ax.set_xticks(x)
    ax.set_xticklabels(countries, rotation=45, ha='right', fontsize=18)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel('TE frequency', fontsize=16)
    ax.axhline(0.5, color='gray', linestyle=':', linewidth=1)
    ax.tick_params(axis='y', labelsize=16)


def plot_gene_region_on_ax(ax, gene_symbol, gene_df, exons, strand, chrom, te_lanes, n_te_lanes):
    """
    Draw one gene's region onto the given axis:
    - backbone line at y=0 with exon boxes on top of it
    - a small direction arrow centered on the gene's midpoint
    - a TE track above showing every TE inserted in the gene (colored by
      main class, solid edge for 'ref', dashed edge for 'nonref'),
      stacked into lanes (te_lanes/n_te_lanes) so overlapping TEs don't
      sit on top of each other. If there's more than one TE, each gets a
      number label under its box matching the frequency bar plot.
    - the gene name above everything
    x-axis is shown in Mb rather than raw bp.
    """
    gene_start = gene_df['gene_start'].iloc[0]
    gene_end   = gene_df['gene_end'].iloc[0]
    gene_length = gene_end - gene_start
    margin = gene_length * 0.1

    region_start = min(gene_start, gene_df['te_start'].min()) - margin
    region_end   = max(gene_end,   gene_df['te_end'].max())   + margin

    gene_name = gene_df['gene_name_updated'].iloc[0]
    if pd.isna(gene_name):
        gene_name = gene_symbol

    gene_y = 0

    # backbone line for the gene body (exons sit on top of this)
    ax.plot([gene_start, gene_end], [gene_y, gene_y], color='black', linewidth=1.5, zorder=1)

    # each exon as a black box, drawn on top of the backbone
    exon_height = 0.3
    for exon in exons:
        ax.add_patch(plt.Rectangle(
            (exon['start_position'], gene_y - exon_height / 2),
            exon['end_position'] - exon['start_position'],
            exon_height,
            facecolor='black',
            zorder=2
        ))

    # small direction arrow, centered on the gene's midpoint
    arrow_y = gene_y + 0.45
    mid = (gene_start + gene_end) / 2
    arrow_half = gene_length * 0.04
    if strand == '-':
        arrow_start, arrow_end = mid + arrow_half, mid - arrow_half
    else:
        arrow_start, arrow_end = mid - arrow_half, mid + arrow_half

    ax.annotate(
        '', xy=(arrow_end, arrow_y), xytext=(arrow_start, arrow_y),
        arrowprops=dict(arrowstyle='-|>', color='black', lw=1.5, mutation_scale=15),
        zorder=3
    )

    # TE track: every TE inserted in this gene, colored by main class,
    # solid edge = ref, dashed edge = nonref, stacked into lanes so
    # overlapping TEs are still individually visible
    n_te = len(gene_df)
    te_base_y = 0.85
    lane_height = 0.32
    te_height = 0.25
    for te_i, (lane, (_, te_row)) in enumerate(zip(te_lanes, gene_df.iterrows())):
        te_s, te_e = te_row['te_start'], te_row['te_end']
        te_class = get_te_main_class(get_te_name(te_row))
        color = TE_CLASS_COLORS.get(te_class, TE_CLASS_COLORS['other'])
        linestyle = 'dashed' if te_row['origin'] == 'nonref' else 'solid'
        te_y = te_base_y + lane * lane_height
        te_width = max(te_e - te_s, gene_length * 0.001)  # guarantee a visible sliver for tiny TEs
        ax.add_patch(plt.Rectangle(
            (te_s, te_y - te_height / 2),
            te_width,
            te_height,
            facecolor=color,
            edgecolor='black',
            linestyle=linestyle,
            linewidth=1.2,
            zorder=2
        ))
        if n_te > 1:
            ax.text(te_s + te_width / 2, te_y - te_height / 2 - 0.06, str(te_i + 1),
                     ha='center', va='top', fontsize=12, fontweight='bold')

    te_top = te_base_y + max(n_te_lanes - 1, 0) * lane_height

    # gene name above the TE stack
    ax.text(mid, te_top + te_height / 2 + 0.15, gene_name,
             ha='center', va='bottom', fontsize=17, fontweight='bold')

    ax.set_xlim(region_start, region_end)
    ax.set_ylim(-1, te_top + 1.0)
    ax.set_yticks([])

    # x-axis in Mb instead of raw bp
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x / 1e6:.2f}'))
    chrom_label = chrom if chrom else ''
    ax.set_xlabel(f'{chrom_label} position (Mb)', fontsize=17)
    ax.tick_params(axis='both', labelsize=16)


import string


def plot_all_genes(significant_sweep_df, gff_full):
    """
    Plot every gene's region stacked one on top of the other in a single
    figure, ordered by chromosome (chr 1 first, then chr 2, chr 3, ...).
    Each row has two panels: the gene region (with TEs, stacked into
    lanes when they overlap) on the left, and a frequency bar plot
    (per-TE, per-country) on the right. Every panel gets a bold letter
    label (A, B, C, ...) in reading order. One shared legend sits on the
    side of the whole figure. Panel heights scale with how many TE lanes
    each gene needs, so dense panels like VGSC get proportionally more
    room.
    """
    gene_groups = []
    for gene_symbol, gene_df in significant_sweep_df.groupby('gene_symbol'):
        chrom = get_gene_chromosome(gff_full, gene_symbol)
        gene_groups.append((chrom, gene_symbol, gene_df))
    gene_groups.sort(key=lambda item: (str(item[0]), item[2]['gene_start'].iloc[0]))

    # precompute TE lane packing per gene so panel heights can be set
    # before the figure is created
    panel_info = []
    for chrom, gene_symbol, gene_df in gene_groups:
        te_rows = list(zip(gene_df['te_start'], gene_df['te_end']))
        lanes, n_lanes = compute_te_lanes(te_rows)
        panel_info.append((chrom, gene_symbol, gene_df, lanes, n_lanes))

    n_genes = len(panel_info)
    height_ratios = [2 + 0.6 * n_lanes for *_, n_lanes in panel_info]
    total_height = sum(height_ratios) * 1.4

    fig, axes = plt.subplots(nrows=n_genes, ncols=2, figsize=(30, total_height),
                              gridspec_kw={'height_ratios': height_ratios, 'width_ratios': [2.6, 2.4]})
    if n_genes == 1:
        axes = axes.reshape(1, 2)

    panel_letters = iter(string.ascii_uppercase)

    for row_i, (chrom, gene_symbol, gene_df, lanes, n_lanes) in enumerate(panel_info):
        gene_ax = axes[row_i, 0]
        freq_ax = axes[row_i, 1]
        exons  = get_exons_for_gene(gff_full, gene_symbol)
        strand = get_gene_strand(gff_full, gene_symbol)
        plot_gene_region_on_ax(gene_ax, gene_symbol, gene_df, exons, strand, chrom, lanes, n_lanes)
        draw_te_frequency_bars(freq_ax, gene_df)

        for ax in (gene_ax, freq_ax):
            letter = next(panel_letters)
            ax.text(-0.05, 1.08, letter, transform=ax.transAxes,
                     fontsize=22, fontweight='bold', va='bottom', ha='right')

    present_classes = {
        get_te_main_class(get_te_name(te_row))
        for _, te_row in significant_sweep_df.iterrows()
    }
    class_handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=color, edgecolor='black', label=cls)
        for cls, color in TE_CLASS_COLORS.items()
        if cls in present_classes
    ]
    origin_handles = []
    if 'ref' in set(significant_sweep_df['origin']):
        origin_handles.append(
            plt.Rectangle((0, 0), 1, 1, facecolor='white', edgecolor='black',
                           linestyle='solid', label='ref')
        )
    if 'nonref' in set(significant_sweep_df['origin']):
        origin_handles.append(
            plt.Rectangle((0, 0), 1, 1, facecolor='white', edgecolor='black',
                           linestyle='dashed', label='nonref')
        )
    legend_handles = class_handles + origin_handles
    fig.legend(handles=legend_handles,
               loc='upper center', bbox_to_anchor=(0.5, 0.012),
               ncol=max(1, len(legend_handles)),
               frameon=False, fontsize=20)

    plt.tight_layout(rect=(0, 0.02, 1, 1))
    output_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        'all_genes_regions.png'
    )
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f'Saved figure to: {output_path}')
    plt.close(fig)

# define margin in base pairs
margin_bp  = 10000
# reading the annotation
path_annotation = '/work/users/g/a/gabivla/lab/SV_mosquitos/TEs_Aedes_aegypti/Annotation_VectorBase/VectorBase-54_AaegyptiLVP_AGWG.gff'
gff_full = read_gff(path_annotation)
# Gene coordinates
gene_annot = gff_full[gff_full['type'] == 'protein_coding_gene'].copy()
gene_annot['gene_symbol'] = gene_annot['attributes'].str.extract(r'ID=([^;]+)')[0]
'''# loading the te information
path_ref = '../IR_genes_Freq/plot_IR_genes/getting_TEs_info/TE_information_all_genes_overlap_reference.csv'
path_nonref = '../IR_genes_Freq/plot_IR_genes/getting_TEs_info/TE_information_all_genes_overlap_nonreference.csv'
df_ref = pd.read_csv(path_ref)
df_nonref = pd.read_csv(path_nonref)
df_ref    = load_te_csv(path_ref, 'ref')
df_nonref = load_te_csv(path_nonref, 'nonref')
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
# TE coordinates
df[['te_start', 'te_end']] = df['TE_information'].apply(
    lambda x: pd.Series(parse_te_coords(x))
)
# need to find the ones that are significant according to sweepfinder data
sweep_info_path = '../IR_genes_Freq/plot_sweep_findeer/updated_strategy_analysis/make_new_sweep_files/new_sweeps_files/'

# Rows (one per TE) that have at least one significant CLR hit inside their
# own coordinates get collected here, then concatenated into the final df.
significant_dfs = []

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
    chrom = get_te_chromosome(tes['TE_information'].iloc[0])
    region_start = min(gene_start, tes['te_start'].min()) - margin_bp
    region_end   = max(gene_end,   tes['te_end'].max())   + margin_bp
    sweep_data = load_sweep_data(sweep_info_path, chrom, region_start, region_end)

    # For every TE belonging to this gene, check each population for
    # above-threshold (significant) sweep locations that fall inside that
    # TE's own coordinates, and record the actual CLR/LR value(s) found.
    te_sig_values = []
    te_has_sweep  = []
    for _, te_row in tes.iterrows():
        te_s, te_e = int(te_row['te_start']), int(te_row['te_end'])
        sig_entries = []
        for pop, data in sweep_data.items():
            hits = data['threshold_region'][
                (data['threshold_region']['location'] >= te_s) &
                (data['threshold_region']['location'] <= te_e)
            ]
            for lr_val in hits['LR']:
                sig_entries.append(f"{pop}:{lr_val}")
        te_sig_values.append(';'.join(sig_entries))
        te_has_sweep.append(len(sig_entries) > 0)

    tes['highest_frequency_country_CLR'] = te_sig_values
    tes['has_sweep'] = te_has_sweep

    # Keep only the TE/gene rows that actually had a significant CLR hit
    sig_tes = tes[tes['has_sweep']].copy()
    if not sig_tes.empty:
        significant_dfs.append(sig_tes)

    # Determine which countries have a sweep above threshold inside ANY TE
    # for this gene (kept from the original script for the summary print)
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

# Final dataframe: only TE/gene rows with at least one significant CLR hit,
# with the new highest_frequency_country_CLR column (';'-joined "pop:CLR"
# values when more than one population/hit is significant for that TE).
if significant_dfs:
    significant_sweep_df = pd.concat(significant_dfs, ignore_index=True)
else:
    significant_sweep_df = pd.DataFrame(columns=list(df.columns) +
                                         ['highest_frequency_country_CLR', 'has_sweep'])

# make all frequency_ columns numeric - non-numeric values like "missing_info" become NaN
for col in significant_sweep_df.columns:
    if col.startswith('frequency_'):
        significant_sweep_df[col] = pd.to_numeric(significant_sweep_df[col], errors='coerce')

# get only the ones that have a sweep - are in the significant_sweep_df and the sweep is correspondent to a country with frequency > 0.5
significant_sweep_df['sweep_country'] = significant_sweep_df['highest_frequency_country_CLR'].str.split(';').apply(
    lambda x: [entry.split(':')[0] for entry in x if entry]
)
# Rio_Claro is a Colombia population - treat any Rio_Claro sweep hit as a Colombia hit
significant_sweep_df['sweep_country'] = significant_sweep_df['sweep_country'].apply(
    lambda countries: ['Colombia' if c == 'Rio_Claro' else c for c in countries]
)
# among all the columns that start with frequency_* - get the country or countries
# tied for the highest frequency
freq_cols = [c for c in significant_sweep_df.columns if c.startswith('frequency_') and c != 'frequency_diff']
significant_sweep_df['freq_country'] = significant_sweep_df[freq_cols].apply(
    lambda row: [
        c.replace('frequency_', '')
        for c in freq_cols
        if pd.notna(row[c]) and row[c] == row[freq_cols].max(skipna=True)
    ],
    axis=1
)
print(significant_sweep_df)
# now keep the ones where freq_country and sweep_country are the same
significant_sweep_df = significant_sweep_df[significant_sweep_df.apply(
    lambda row: any(c in row['sweep_country'] for c in row['freq_country']), axis=1
)]
# drop the temporary columns
significant_sweep_df = significant_sweep_df.drop(columns=['sweep_country', 'freq_country'])
significant_sweep_df.to_csv('TEs_with_significant_sweeps.csv', index=False)'''

significant_sweep_df = pd.read_csv('TEs_with_significant_sweeps.csv')
# getting info for each exon
for gene_symbol, gene_df in significant_sweep_df.groupby('gene_symbol'):
    exons = get_exons_for_gene(gff_full, gene_symbol)
    significant_sweep_df.loc[gene_df.index, 'exons'] = str(exons)
    print(f"Gene {gene_symbol} has {len(exons)} exons.")

plot_all_genes(significant_sweep_df, gff_full)