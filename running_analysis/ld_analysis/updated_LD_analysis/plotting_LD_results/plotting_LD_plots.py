import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch
import pandas as pd
import numpy as np

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 10,
    'pdf.fonttype': 42,
    'svg.fonttype': 'none',
})

# ── Colour palette ──────────────────────────────────────────────────────────
C_SNP_VGSC = '#2F72B8'
C_SNP_RDL  = '#E42626'
C_EXON     = '#111111'

# TE colours, split by (type, LD_direction) — 4 per gene
C_VGSC_REF_POS    = '#6C3483'   # dark purple  – Reference, positive LD
C_VGSC_REF_NEG    = '#D7BDE2'   # light purple – Reference, negative LD
C_VGSC_NONREF_POS = '#1B4F72'   # dark blue    – Non-Reference, positive LD
C_VGSC_NONREF_NEG = '#AED6F1'   # light blue   – Non-Reference, negative LD

C_RDL_REF_POS     = '#922B21'   # dark red     – Reference, positive LD
C_RDL_REF_NEG     = '#F1948A'   # light red    – Reference, negative LD
C_RDL_NONREF_POS  = '#AD1457'   # dark pink    – Non-Reference, positive LD
C_RDL_NONREF_NEG  = '#F8BBD0'   # light pink   – Non-Reference, negative LD

C_UNKNOWN_LD = '#B0B0B0'   # gray – no overlapping non-missing genotype calls

VGSC_COLOR_MAP = {
    ('Reference', 'positive'):     C_VGSC_REF_POS,
    ('Reference', 'negative'):     C_VGSC_REF_NEG,
    ('Non-Reference', 'positive'): C_VGSC_NONREF_POS,
    ('Non-Reference', 'negative'): C_VGSC_NONREF_NEG,
    ('Reference', 'unknown'):      C_UNKNOWN_LD,
    ('Non-Reference', 'unknown'):  C_UNKNOWN_LD,
}
RDL_COLOR_MAP = {
    ('Reference', 'positive'):     C_RDL_REF_POS,
    ('Reference', 'negative'):     C_RDL_REF_NEG,
    ('Non-Reference', 'positive'): C_RDL_NONREF_POS,
    ('Non-Reference', 'negative'): C_RDL_NONREF_NEG,
    ('Reference', 'unknown'):      C_UNKNOWN_LD,
    ('Non-Reference', 'unknown'):  C_UNKNOWN_LD,
}

C_VGSC_OUTER = '#EAF4FD'
C_VGSC_INNER = '#C8DFF4'
C_RDL_OUTER  = '#FDECEF'
C_RDL_INNER  = '#F4CAD3'

GENE_Y, GENE_TOP = 0.35, 0.82
GENE_H = GENE_TOP - GENE_Y
PANEL_Y0, PANEL_Y1 = 0.12, 0.96


# ══════════════════════════════════════════════════════════════════════════
# EXON LOADER
# ══════════════════════════════════════════════════════════════════════════
def get_exons(gff_path, gene_id):
    """Extract exon start/end positions for a given gene ID from a GFF file."""
    gff_df = pd.read_csv(gff_path, sep='\t', comment='#', header=None)
    exon_df = gff_df[gff_df[2] == 'exon']
    gene_exon_df = exon_df[exon_df[8].str.contains(gene_id)].copy()
    gene_exon_df['exon_info'] = gene_exon_df[8].str.extract(r'ID=([^;]+)')[0].tolist()
    return [{'exon_name': row['exon_info'],
             'start_position': int(row[3]),
             'end_position': int(row[4])}
            for _, row in gene_exon_df.iterrows()]


# ══════════════════════════════════════════════════════════════════════════
# GENE + TE + EXON PLOTTER
# ══════════════════════════════════════════════════════════════════════════
def plot_gene_TEs(ax, panel_x0, panel_x1, gene_x0, gene_x1,
                   exons, te_df, snps,
                   title, outer_bg_color, gene_fill_color, border_color,
                   snp_color, color_map,
                   legend_entries, legend_dark_pairs=None):
    """
    Draw one gene panel with exons, SNP lines, and TE bars styled by
    (type, LD_direction) — 4 colours per gene — with r squared shown in each label.

    Parameters
    ----------
    te_df : DataFrame with columns TE_name, position_start, position_end,
            type ('Reference' / 'Non-Reference'), LD_direction
            ('positive' / 'negative'), r_squared (float).
            This is the TE-only output of checking_genotypes(), filtered
            to one gene region.
    snps  : list of dicts with TE_name, start_position, end_position
    color_map : dict keyed by (type, LD_direction) -> hex color, e.g.
            {('Reference', 'positive'): '#6C3483', ...}
    """
    span = panel_x1 - panel_x0

    # Coloured panel background
    ax.add_patch(FancyBboxPatch(
        (panel_x0, PANEL_Y0), span, PANEL_Y1 - PANEL_Y0,
        boxstyle='round,pad=0.006,rounding_size=0.22',
        lw=2.0, edgecolor=border_color, facecolor=outer_bg_color,
        zorder=1, clip_on=False))

    # White gene body (shows intergenic flanks)
    inner_pad = span * 0.005
    ax.add_patch(FancyBboxPatch(
        (panel_x0 + inner_pad, GENE_Y),
        span - 2 * inner_pad, GENE_H,
        boxstyle='round,pad=0.004,rounding_size=0.08',
        lw=1.1, edgecolor='#444', facecolor='white',
        zorder=2))

    # Gene body fill (only over the actual gene span)
    ax.add_patch(FancyBboxPatch(
        (gene_x0, GENE_Y), gene_x1 - gene_x0, GENE_H,
        boxstyle='round,pad=0,rounding_size=0.06',
        lw=0, edgecolor='none', facecolor=gene_fill_color,
        zorder=3))

    # Exons
    for ex in exons:
        ax.add_patch(patches.Rectangle(
            (ex['start_position'], GENE_Y),
            ex['end_position'] - ex['start_position'], GENE_H,
            lw=0, fc=C_EXON, zorder=4))

    # ── TE bars, styled by (type, LD_direction) ────────────────────────────
    te_annot = []
    for _, te in te_df.iterrows():
        start = int(te['position_start'])
        end = int(te['position_end'])
        w = max(end - start, span * 0.003)
        ld = te['LD_direction']
        te_type = te['type']

        face_color = color_map.get((te_type, ld), C_UNKNOWN_LD)

        if ld == 'negative':
            edge_color = color_map.get((te_type, 'positive'), C_UNKNOWN_LD)
            linestyle, lw = '--', 1.4
        elif ld == 'positive':
            edge_color, linestyle, lw = face_color, '-', 0.0
        else:  # 'unknown' or missing
            edge_color, linestyle, lw = '#666666', ':', 1.0

        ax.add_patch(patches.Rectangle(
            (start, GENE_Y), w, GENE_H,
            lw=lw, ec=edge_color, fc=face_color, ls=linestyle,
            zorder=5, alpha=0.92))

        te_annot.append({
            'TE_name': te['TE_name'],
            'start_position': start,
            'end_position': end,
            'LD_direction': ld,
            'type': te_type,
            'r_squared': te['r_squared'],
            '_kind': 'te'
        })

    # SNP / mutation lines
    for snp in snps:
        cx = (snp['start_position'] + snp['end_position']) / 2
        ax.plot([cx, cx], [GENE_Y, GENE_TOP], color=snp_color, lw=2.0, zorder=6)

    # Panel title
    ax.text((panel_x0 + panel_x1) / 2, PANEL_Y1 - 0.01,
            title, ha='center', va='top',
            fontsize=13, fontstyle='italic', fontweight='bold',
            color='#1a1a1a', zorder=8)

    # ── Annotation arrows + labels below panel ─────────────────────────────
    all_annot = sorted(
        [dict(e, _kind='snp') for e in snps] + te_annot,
        key=lambda e: e['start_position'])

    ROW_H = 0.095
    arrow_tip_y = PANEL_Y0 - 0.004

    for ri, e in enumerate(all_annot):
        if e['_kind'] == 'snp':
            col = snp_color
        else:
            col = color_map.get((e['type'], e['LD_direction']), C_UNKNOWN_LD)

        cx = (e['start_position'] + e['end_position']) / 2
        label_y = arrow_tip_y - 0.040 - ri * ROW_H

        ax.plot([cx, cx], [GENE_Y - 0.005, PANEL_Y0],
                color=col, lw=1.2, ls='--', zorder=6, clip_on=False)
        ax.annotate('', xy=(cx, arrow_tip_y), xytext=(cx, label_y + 0.035),
                    arrowprops=dict(arrowstyle='->', color=col,
                                     lw=1.6, mutation_scale=11),
                    zorder=7, annotation_clip=False)

        name = e['TE_name']
        label = name if e['_kind'] == 'snp' else f"{name} (r\u00b2={e['r_squared']:.2f})"
        ax.text(cx, label_y, label,
                ha='right', va='top', fontsize=9, color=col,
                fontweight='bold' if e['_kind'] == 'snp' else 'normal',
                clip_on=False, zorder=7)

    # ── Legend ──────────────────────────────────────────────────────────────
    legend_dark_pairs = legend_dark_pairs or {}
    lx, ly = panel_x1 + span * 0.055, PANEL_Y1 - 0.04
    sw, gap = span * 0.028, span * 0.014
    for (lbl, lcol, lstyle) in legend_entries:
        if lstyle == 'line':
            ax.plot([lx, lx + sw], [ly, ly], color=lcol, lw=3.5,
                    solid_capstyle='butt', clip_on=False, zorder=8)
        elif lstyle == 'dashed_patch':
            # lcol is the light fill; pair it with its own type's dark color as the outline
            edge = legend_dark_pairs.get(lbl, '#444')
            ax.add_patch(patches.Rectangle(
                (lx, ly - 0.030), sw, 0.060,
                lw=1.4, ec=edge, fc=lcol, ls='--', zorder=8, clip_on=False))
        elif lstyle == 'dotted_patch':
            ax.add_patch(patches.Rectangle(
                (lx, ly - 0.030), sw, 0.060,
                lw=1.0, ec='#666666', fc=lcol, ls=':', zorder=8, clip_on=False))
        else:
            ec = '#888' if lcol == 'white' else '#444'
            ax.add_patch(patches.Rectangle(
                (lx, ly - 0.030), sw, 0.060,
                lw=0.8, ec=ec, fc=lcol, zorder=8, clip_on=False))
        ax.text(lx + sw + gap, ly, lbl, ha='left', va='center',
                fontsize=10, color='#222', fontweight='bold',
                clip_on=False, zorder=8)
        ly -= 0.120


# ══════════════════════════════════════════════════════════════════════════
# GENOTYPE / LD DIRECTION / R² EXTRACTION (reads your actual CSVs)
# ══════════════════════════════════════════════════════════════════════════
import ast

def ld_direction(ref_genos, te_genos):
    """Compare two numeric genotype arrays and call the LD sign."""
    mask = ~np.isnan(ref_genos) & ~np.isnan(te_genos)
    if mask.sum() == 0:
        return np.nan, 'unknown'
    agreement = np.mean(ref_genos[mask] == te_genos[mask])
    direction = 'positive' if agreement > 0.5 else 'negative'
    return agreement, direction


def parse_genotypes(geno_str):
    """Turn the stringified list into numeric codes for comparison."""
    if pd.isna(geno_str):
        return None
    genos = ast.literal_eval(geno_str)
    mapping = {'0/0': 0, '1/1': 1, '0/1': 0.5, '1/0': 0.5, './.': np.nan}
    return np.array([mapping.get(g, np.nan) for g in genos], dtype=float)


def checking_genotypes(list_dfs, path_df_ref, path_df_nonref,
                        snp_col_1, snp_col_2):
    """
    Merge significant-r² dataframes with genotype calls, filter out the SNP
    marker rows themselves, and compute LD direction + r_squared for every
    remaining TE row.

    For each TE, the reference-SNP encoding used for the LD comparison is
    chosen from whichever of encoding_1 / encoding_2 is actually *available*
    (has real genotype data in the matching genotype_info file) — preferring
    the one with the higher r² (from the significant_r_squared_* df) when
    both encodings are available. An encoding with a significant r² but no
    genotype data for this TE set is skipped rather than crashing.
    """
    genotypes_ref = pd.read_csv(path_df_ref)
    genotypes_nonref = pd.read_csv(path_df_nonref)
    genotypes_ref['TEs'] = genotypes_ref['TE']
    genotypes_nonref['TEs'] = genotypes_nonref['TE']
    genotypes_ref = genotypes_ref.drop(columns=['TE', 'Unnamed: 0'])
    genotypes_nonref = genotypes_nonref.drop(columns=['TE', 'Unnamed: 0'])

    results = []
    for df in list_dfs:
        is_ref = df['type'].iloc[0] == 'Reference'
        genotype_table = genotypes_ref if is_ref else genotypes_nonref
        df = df.merge(genotype_table, on='TEs', how='left')

        snp_encoding_1_name = snp_col_1
        snp_encoding_2_name = snp_col_2

        v_row_1 = df[df['TEs'] == snp_encoding_1_name]
        v_row_2 = df[df['TEs'] == snp_encoding_2_name]

        # None if the encoding is missing from this genotype table, or
        # present but with a blank/NaN Genotypes cell.
        ref_genos_1 = (parse_genotypes(v_row_1['Genotypes'].iloc[0])
                        if not v_row_1.empty else None)
        ref_genos_2 = (parse_genotypes(v_row_2['Genotypes'].iloc[0])
                        if not v_row_2.empty else None)

        te_mask = ~df['TEs'].str.endswith(('encoding_1', 'encoding_2'))
        te_only_df = df[te_mask].copy()
        agreements, directions, used_encoding, r_squared_values = [], [], [], []

        for _, row in te_only_df.iterrows():
            te_genos = parse_genotypes(row['Genotypes'])
            r2_1 = row.get(snp_encoding_1_name, np.nan)
            r2_2 = row.get(snp_encoding_2_name, np.nan)

            # only treat an encoding as usable if it has BOTH a real r²
            # value from the significance df AND actual genotype data
            candidates = []
            if ref_genos_1 is not None and pd.notna(r2_1):
                candidates.append(('encoding_1', ref_genos_1, r2_1))
            if ref_genos_2 is not None and pd.notna(r2_2):
                candidates.append(('encoding_2', ref_genos_2, r2_2))

            if candidates:
                enc, ref_genos, r2_value = max(candidates, key=lambda c: c[2])
            else:
                ref_genos, enc, r2_value = None, 'unknown', np.nan

            if ref_genos is None or te_genos is None:
                agreements.append(np.nan)
                directions.append('unknown')
                used_encoding.append(enc)
                r_squared_values.append(r2_value)
                continue

            agreement, direction = ld_direction(ref_genos, te_genos)
            agreements.append(agreement)
            directions.append(direction)
            used_encoding.append(enc)
            r_squared_values.append(r2_value)

        te_only_df['agreement_fraction'] = agreements
        te_only_df['LD_direction'] = directions
        te_only_df['compared_to_encoding'] = used_encoding
        te_only_df['r_squared'] = r_squared_values

        results.append(te_only_df)

    merged = pd.concat(results, ignore_index=True)

    pattern = r'(AaegL5_\d)_(\d+)_(\d+)_(.+)'
    merged[['chromosome', 'position_start', 'position_end', 'TE_name']] = \
        merged['TEs'].str.split(pattern, expand=True).iloc[:, 1:5]
    merged['position_start'] = merged['position_start'].astype(int)
    merged['position_end'] = merged['position_end'].astype(int)

    return merged


# ══════════════════════════════════════════════════════════════════════════
# LD / TE / GENE SCHEME PLOTTER — core drawing routine (one gene "column")
# ══════════════════════════════════════════════════════════════════════════
def _plot_gene_column(ax_top, ax_snp, ax_bottom,
                       results, exons, snps, gene_name,
                       panel_x0, panel_x1, chromosome,
                       gene_start, gene_end, letter_offset=0):
    """
    Draw the LD/TE panel, optional variant (SNP) panel, and gene scheme
    panel into a set of already-created axes (one "column" of a figure).

    This holds all the actual drawing logic; both the single-gene
    ``plot_LD_scheme`` and the multi-gene ``plot_LD_scheme_combined``
    call into this.

    Parameters
    ----------
    ax_top, ax_bottom : matplotlib Axes
        Axes for the LD/TE panel and the gene-scheme panel respectively.
    ax_snp : matplotlib Axes or None
        Axes for the no-r² variant panel. Pass None if there are no
        SNPs without r² for this gene (in which case no variant panel
        is drawn, and ax_bottom is used directly below ax_top).
    letter_offset : int
        0 -> panels are lettered starting at 'A'; 3 -> starting at 'D'; etc.
        Lets multiple genes in one figure get consecutive letters
        (gene 1: A, B, C -> gene 2: D, E, F).

    Returns
    -------
    n_letters : int
        How many panel letters this gene used (2 or 3), so the caller
        can advance its letter cursor for the next gene.
    """

    def _letter(i):
        return chr(ord('A') + letter_offset + i)

    # =========================================================
    # Get TE information
    # =========================================================

    main_classes = ['SINE', 'TIR', 'LINE', 'Unknown']

    # Make a copy so the original dataframe is not modified
    results = results.copy()

    # ---------------------------------------------------------
    # Identify the main TE class from TE_name
    # ---------------------------------------------------------

    def get_main_class(te_name):
        if pd.isna(te_name):
            return 'other'

        # Split TE name on "_"
        name_parts = str(te_name).split("_")

        # Check whether any part matches a main class
        for main_class in main_classes:
            if main_class in name_parts:
                return main_class

        return 'other'

    # Use the actual TE_name to determine the class
    results['main_class'] = results['TE_name'].apply(get_main_class)

    # ---------------------------------------------------------
    # Extract TE information
    # ---------------------------------------------------------

    te_positions = results[['position_start', 'position_end']].values
    te_names = results['main_class'].values
    te_r_squared = results['r_squared'].values
    te_directions = results['LD_direction'].values
    if 'type' in results.columns:
        te_types = results['type'].values
    else:
        te_types = np.array(['Reference'] * len(results))

    # ---------------------------------------------------------
    # Split the SNPs/variants into those that DO carry a usable
    # r² value (plotted with the TEs) and those that DON'T
    # (plotted in their own panel).
    # ---------------------------------------------------------

    def _has_r2(snp):
        r2 = snp.get('r_squared', None)
        return r2 is not None and not pd.isna(r2)

    snps_with_r2 = [s for s in snps if _has_r2(s)]
    snps_without_r2 = [s for s in snps if not _has_r2(s)]

    show_variant_panel = ax_snp is not None and len(snps_without_r2) > 0

    # =========================================================
    # Label placement helper — keeps text horizontal and spreads
    # labels across vertical "tiers" so close-together points
    # don't collide, without ever moving a label off its point's
    # true x position (so genomic position stays unambiguous).
    # =========================================================

    def _stagger_labels(ax, items, base_gap, tier_step, x_span,
                         fontsize=7, min_sep_frac=0.040, max_tiers=5):
        """
        items: list of dicts with 'x' (Mb), 'y' (baseline, e.g. the r²
            value or marker y), 'label', 'color', optional 'weight'.
        Places each label directly above its point, in the lowest
        vertical tier that doesn't collide (in x) with an already-
        placed label in that tier. Returns the max tier index used,
        so the caller can extend the y-axis to make room.
        """
        min_sep = x_span * min_sep_frac
        last_x_per_tier = {}
        max_tier_used = 0

        for item in sorted(items, key=lambda d: d['x']):
            tier = 0
            while tier < max_tiers and (item['x'] - last_x_per_tier.get(tier, -1e18)) < min_sep:
                tier += 1
            last_x_per_tier[tier] = item['x']
            max_tier_used = max(max_tier_used, tier)

            label_y = item['y'] + base_gap + tier * tier_step

            ax.text(
                item['x'], label_y, item['label'],
                ha='center', va='bottom',
                fontsize=fontsize,
                fontweight=item.get('weight', 'normal'),
                color=item.get('color', 'black'),
                rotation=0,
                zorder=6,
            )

            # thin leader line so a label pushed up to a higher tier
            # still reads unambiguously as belonging to its point
            if tier > 0:
                ax.plot([item['x'], item['x']], [item['y'] + base_gap * 0.3, label_y - 0.01],
                        color=item.get('color', 'black'), lw=0.5, alpha=0.5, zorder=5)

        return max_tier_used

    # =========================================================
    # Convert genomic coordinates: bp -> Mb
    # =========================================================

    panel_x0_mb = panel_x0 / 1e6
    panel_x1_mb = panel_x1 / 1e6
    x_span_mb = panel_x1_mb - panel_x0_mb

    gene_start_mb = gene_start / 1e6
    gene_end_mb = gene_end / 1e6
    bottom_panel_y0, bottom_panel_y1 = 0.05, 0.95

    # =========================================================
    # TOP PANEL: LD / TEs (+ any SNPs that have their own r²)
    # =========================================================

    ax_top.set_xlim(panel_x0_mb, panel_x1_mb)

    # r² is limited to the top half of the scale (0.5–1), where the
    # signal of interest lives
    TOP_YMIN = 0.5
    TOP_YMAX = 1.0

    # ---------------------------------------------------------
    # Draw outer rectangle
    # ---------------------------------------------------------

    top_rect = patches.Rectangle(
        (panel_x0_mb, TOP_YMIN),
        panel_x1_mb - panel_x0_mb,
        TOP_YMAX - TOP_YMIN,
        linewidth=1,
        edgecolor="black",
        facecolor="none",
        zorder=10
    )

    ax_top.add_patch(top_rect)

    # Light reference gridlines behind the data
    for gy in (0.5, 0.6, 0.7, 0.8, 0.9, 1.0):
        ax_top.axhline(gy, color='#dddddd', lw=0.6, zorder=0)

    # ---------------------------------------------------------
    # Colors for LD direction
    # Direction ONLY controls color
    # ---------------------------------------------------------

    # Okabe-Ito colorblind-safe pair
    positive_color = "#0072B2"   # blue
    negative_color = "#D55E00"   # vermillion

    def _direction_color(direction):
        if direction == 'positive':
            return positive_color
        elif direction == 'negative':
            return negative_color
        return "gray"

    # ---------------------------------------------------------
    # Height of each TE rectangle
    # ---------------------------------------------------------

    te_height = 0.05

    # ---------------------------------------------------------
    # Plot each TE, collecting label placement info as we go
    # ---------------------------------------------------------

    top_label_items = []

    # Minimum visual width so very short TEs are still visible at
    # this genomic scale
    min_te_width_mb = x_span_mb * 0.004

    for position, te_name, r2, direction, te_type in zip(
        te_positions,
        te_names,
        te_r_squared,
        te_directions,
        te_types
    ):

        # -----------------------------------------------------
        # Convert TE coordinates bp -> Mb
        # -----------------------------------------------------

        te_start_mb = position[0] / 1e6
        te_end_mb = position[1] / 1e6

        # Make sure start < end
        if te_end_mb < te_start_mb:
            te_start_mb, te_end_mb = te_end_mb, te_start_mb

        te_width_mb = max(te_end_mb - te_start_mb, min_te_width_mb)
        te_end_mb = te_start_mb + te_width_mb

        # -----------------------------------------------------
        # r² value
        # -----------------------------------------------------

        r2 = float(r2)

        # Keep r² within 0-1
        r2 = max(0, min(r2, 1))

        # This panel only shows r² in [0.5, 1] — skip anything below
        # that so nothing gets drawn (or labeled) outside the visible box
        if r2 < TOP_YMIN:
            continue

        # -----------------------------------------------------
        # Direction ONLY determines color
        # -----------------------------------------------------

        color = _direction_color(direction)

        # -----------------------------------------------------
        # Shape ONLY determines TE type:
        #   Reference     -> rectangle
        #   Non-Reference -> triangle
        # Both are drawn with the same footprint (same x-span and
        # the same height), positioned AT the r² value.
        #
        # For example:
        # r² = 0.7 -> shape centered at y = 0.7
        # r² = 0.3 -> shape centered at y = 0.3
        # r² = 0.9 -> shape centered at y = 0.9
        # -----------------------------------------------------

        rect_bottom = r2 - te_height / 2
        rect_top = r2 + te_height / 2

        if te_type == 'Non-Reference':
            te_shape = patches.Polygon(
                [(te_start_mb, rect_bottom),
                 (te_end_mb, rect_bottom),
                 ((te_start_mb + te_end_mb) / 2, rect_top)],
                closed=True,
                linewidth=0.5,
                edgecolor="black",
                facecolor=color,
                alpha=0.85,
                zorder=3
            )
        else:
            te_shape = patches.Rectangle(
                (te_start_mb, rect_bottom),
                te_width_mb,
                te_height,
                linewidth=0.5,
                edgecolor="black",
                facecolor=color,
                alpha=0.85,
                zorder=3
            )

        ax_top.add_patch(te_shape)

        # Trace each TE position into the gene scheme using its LD color.
        te_x_mb = (te_start_mb + te_end_mb) / 2
        ax_bottom.plot(
            [te_x_mb, te_x_mb], [bottom_panel_y0, bottom_panel_y1],
            color=color, ls='--', lw=0.7, alpha=0.75, zorder=1,
        )

        top_label_items.append({
            'x': (te_start_mb + te_end_mb) / 2,
            'y': r2 + te_height / 2,
            'label': te_name,
            'color': '#333333',
        })

    # ---------------------------------------------------------
    # Plot any SNPs that DO have their own r² as circles,
    # positioned at their r² value just like the TE rectangles.
    # ---------------------------------------------------------

    snp_marker_size = 65

    for snp in snps_with_r2:
        snp_start_mb = snp['start_position'] / 1e6
        snp_end_mb = snp['end_position'] / 1e6
        snp_x_mb = (snp_start_mb + snp_end_mb) / 2

        r2 = float(snp['r_squared'])
        r2 = max(0, min(r2, 1))

        # This panel only shows r² in [0.5, 1] — skip anything below
        if r2 < TOP_YMIN:
            continue

        color = _direction_color(snp.get('LD_direction', 'unknown'))

        ax_top.scatter(
            snp_x_mb, r2,
            s=snp_marker_size,
            marker='o',
            facecolor=color,
            edgecolor='black',
            linewidth=0.9,
            alpha=0.95,
            zorder=4
        )

        top_label_items.append({
            'x': snp_x_mb,
            'y': r2 + te_height / 2,
            'label': snp['TE_name'],
            'color': 'black',
            'weight': 'bold',
        })

    # ---------------------------------------------------------
    # Place all top-panel labels (TEs + SNPs together) with
    # collision avoidance, horizontal orientation
    # ---------------------------------------------------------

    max_tier_top = _stagger_labels(
        ax_top, top_label_items,
        base_gap=0.025, tier_step=0.11, x_span=x_span_mb, fontsize=7.5,
    )

    # Keep the r² interval identical across gene panels and make it fill
    # the complete top-panel height.
    ax_top.set_ylim(TOP_YMIN, TOP_YMAX)

    # =========================================================
    # Top-panel Y-axis
    # =========================================================
    variant_vgsc = 'V410L'
    variant_rdl = 'A301S'
    variant_name = {
        'VGSC': variant_vgsc,
        'RDL': variant_rdl,
    }.get(gene_name.upper(), gene_name)
    ax_top.set_yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    ax_top.tick_params(axis='both', labelsize=9)
    ax_top.set_ylabel(rf"$r^2$ to {variant_name}", fontsize=11)

    # =========================================================
    # Legend — placed INSIDE the top panel, next to/within the
    # r²-rectangle, rather than floating outside the axes.
    # =========================================================

    positive_patch = patches.Patch(
        facecolor=positive_color, edgecolor='black', linewidth=0.5,
        alpha=0.85, label="Positive LD"
    )

    negative_patch = patches.Patch(
        facecolor=negative_color, edgecolor='black', linewidth=0.5,
        alpha=0.85, label="Negative LD"
    )

    # Shape legend: TE type is shown by shape, independent of color
    ref_shape_handle = plt.Line2D(
        [0], [0],
        marker='s', linestyle='None',
        markerfacecolor='#999999', markeredgecolor='black',
        markersize=8, label='Reference TE'
    )
    nonref_shape_handle = plt.Line2D(
        [0], [0],
        marker='^', linestyle='None',
        markerfacecolor='#999999', markeredgecolor='black',
        markersize=9, label='Non-Reference TE'
    )

    legend_handles = [positive_patch, negative_patch,
                       ref_shape_handle, nonref_shape_handle]

    if snps_with_r2:
        snp_r2_handle = plt.Line2D(
            [0], [0],
            marker='o', linestyle='None',
            markerfacecolor='white', markeredgecolor='black',
            markersize=7, label='Variant (with $r^2$)'
        )
        legend_handles.append(snp_r2_handle)

    ax_top.legend(
        handles=legend_handles,
        loc="upper right",
        frameon=True,
        framealpha=0.9,
        edgecolor='#cccccc',
        fontsize=8,
        handlelength=1.4,
        borderaxespad=0.6,
        labelspacing=0.4,
    )

    # Keep top panel box
    for spine in ax_top.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.9)

    # Panel letter
    ax_top.text(-0.045, 1.02, _letter(0), transform=ax_top.transAxes,
                fontsize=13, fontweight='bold', va='bottom', ha='left')

    # Gene name over the whole column
    ax_top.set_title(gene_name, fontsize=13, fontstyle='italic',
                      fontweight='bold', pad=10)

    # =========================================================
    # MIDDLE PANEL (optional): variants with no r² value,
    # drawn as circles on their own since they can't share the
    # r²-based y-axis of the top panel.
    # =========================================================

    if show_variant_panel:

        ax_snp.set_xlim(panel_x0_mb, panel_x1_mb)

        # Draw the variant-panel border without a top edge.
        ax_snp.plot(
            [panel_x0_mb, panel_x0_mb], [0, 1],
            color="black", lw=1, zorder=10,
        )
        ax_snp.plot(
            [panel_x1_mb, panel_x1_mb], [0, 1],
            color="black", lw=1, zorder=10,
        )
        ax_snp.plot(
            [panel_x0_mb, panel_x1_mb], [0, 0],
            color="black", lw=1, zorder=10,
        )

        variant_color = "#222222"
        variant_y = 0.42

        variant_items = []

        for snp in snps_without_r2:
            snp_start_mb = snp['start_position'] / 1e6
            snp_end_mb = snp['end_position'] / 1e6
            snp_x_mb = (snp_start_mb + snp_end_mb) / 2

            ax_snp.scatter(
                snp_x_mb, variant_y,
                s=snp_marker_size,
                marker='o',
                facecolor=variant_color,
                edgecolor='black',
                linewidth=0.9,
                zorder=4
            )

            variant_items.append({
                'x': snp_x_mb, 'y': variant_y,
                'label': snp['TE_name'], 'color': 'black',
            })

            # Dashed guide line down through the gene panel so the
            # variant's genomic position stays easy to trace.
            ax_bottom.plot(
                [snp_x_mb, snp_x_mb],
                [bottom_panel_y0, bottom_panel_y1],
                color='#aaaaaa', ls='--', lw=0.7, zorder=1,
            )

        max_tier_snp = _stagger_labels(
            ax_snp, variant_items,
            base_gap=0.10, tier_step=0.30, x_span=x_span_mb, fontsize=7.5,
        )

        ax_snp.set_ylim(0, max(1.0, variant_y + 0.20 + max_tier_snp * 0.30 + 0.15))
        ax_snp.set_yticks([])
        ax_snp.set_ylabel("IR related SNPs", fontsize=9.5)

        for spine in ax_snp.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(0.9)

        ax_snp.text(-0.045, 1.05, _letter(1), transform=ax_snp.transAxes,
                    fontsize=13, fontweight='bold', va='bottom', ha='left')

    # =========================================================
    # BOTTOM PANEL: Gene scheme
    # =========================================================

    gene_y = 0.50

    # ---------------------------------------------------------
    # Outer rectangle
    # ---------------------------------------------------------

    panel_rect = patches.Rectangle(
        (panel_x0_mb, bottom_panel_y0),
        panel_x1_mb - panel_x0_mb,
        bottom_panel_y1 - bottom_panel_y0,
        linewidth=1,
        edgecolor="black",
        facecolor="none",
        zorder=10
    )

    ax_bottom.add_patch(panel_rect)

    # ---------------------------------------------------------
    # Gene line
    # ---------------------------------------------------------

    ax_bottom.plot(
        [gene_start_mb, gene_end_mb],
        [gene_y, gene_y],
        color="#1a1a1a",
        linewidth=1.6,
        zorder=1
    )

    # ---------------------------------------------------------
    # Exons
    # ---------------------------------------------------------

    exon_height = 0.16

    for exon in exons:

        exon_start_mb = exon["start_position"] / 1e6
        exon_end_mb = exon["end_position"] / 1e6

        if exon_end_mb < exon_start_mb:
            exon_start_mb, exon_end_mb = exon_end_mb, exon_start_mb

        exon_rect = patches.Rectangle(
            (exon_start_mb, gene_y - exon_height / 2),
            exon_end_mb - exon_start_mb,
            exon_height,
            linewidth=0.5,
            edgecolor="#1a1a1a",
            facecolor="black",
            zorder=2
        )

        ax_bottom.add_patch(exon_rect)

    # ---------------------------------------------------------
    # Bottom X-axis
    # ---------------------------------------------------------

    ax_bottom.set_xlim(panel_x0_mb, panel_x1_mb)
    ax_bottom.set_ylim(0, 1)

    ax_bottom.set_xlabel(
        f"Position on chr {chromosome} (Mb)", fontsize=10.5
    )
    ax_bottom.xaxis.set_major_locator(plt.MaxNLocator(6))
    ax_bottom.xaxis.set_major_formatter(plt.FormatStrFormatter('%.2f'))
    ax_bottom.tick_params(axis='x', labelsize=9)

    # Hide bottom Y-axis
    ax_bottom.set_yticks([])

    ax_bottom.spines["left"].set_visible(False)
    ax_bottom.spines["right"].set_visible(False)
    ax_bottom.spines["top"].set_visible(False)
    ax_bottom.spines["bottom"].set_linewidth(0.9)

    n_letters = 3 if show_variant_panel else 2
    panel_letter = _letter(n_letters - 1)
    ax_bottom.text(-0.045, 1.05, panel_letter, transform=ax_bottom.transAxes,
                    fontsize=13, fontweight='bold', va='bottom', ha='left')

    return n_letters


# ══════════════════════════════════════════════════════════════════════════
# SINGLE-GENE ENTRY POINT (backward compatible)
# ══════════════════════════════════════════════════════════════════════════
def plot_LD_scheme(results, exons, snps, gene_name,
                    panel_x0, panel_x1, chromosome,
                    gene_start, gene_end):
    """
    Draw the LD/TE panel, an optional variant (SNP) panel, and the gene
    scheme panel, stacked vertically, for a single gene. Saves
    ``LD_plots_{gene_name}.png``.

    Parameters are as in ``_plot_gene_column`` (see that function's
    docstring for what ``results``, ``exons`` and ``snps`` need to
    contain).
    """

    def _has_r2(snp):
        r2 = snp.get('r_squared', None)
        return r2 is not None and not pd.isna(r2)

    show_variant_panel = any(not _has_r2(s) for s in snps)

    if show_variant_panel:
        fig, (ax_top, ax_snp, ax_bottom) = plt.subplots(
            3, 1,
            figsize=(11, 7.2),
            sharex=True,
            gridspec_kw={"height_ratios": [1.6, 0.55, 1]}
        )
    else:
        fig, (ax_top, ax_bottom) = plt.subplots(
            2, 1,
            figsize=(11, 6.2),
            sharex=True,
            gridspec_kw={"height_ratios": [1.6, 1]}
        )
        ax_snp = None

    fig.patch.set_facecolor('white')

    _plot_gene_column(
        ax_top, ax_snp, ax_bottom,
        results=results, exons=exons, snps=snps, gene_name=gene_name,
        panel_x0=panel_x0, panel_x1=panel_x1, chromosome=chromosome,
        gene_start=gene_start, gene_end=gene_end, letter_offset=0,
    )

    fig.align_ylabels()
    plt.tight_layout()
    plt.savefig(f'LD_plots_{gene_name}.png', dpi=300, bbox_inches='tight')

    if show_variant_panel:
        return fig, (ax_top, ax_snp, ax_bottom)
    return fig, (ax_top, ax_bottom)


# ══════════════════════════════════════════════════════════════════════════
# MULTI-GENE ENTRY POINT — side-by-side columns in ONE figure
# ══════════════════════════════════════════════════════════════════════════
def plot_LD_scheme_combined(genes, output_name='LD_plots_combined.png',
                             col_width=11, fig_height=7.2):
    """
    Plot several genes' LD/TE + variant + gene-scheme panels side by side
    as columns of ONE figure, with panel letters running consecutively
    across genes (gene 1 -> A, B, C; gene 2 -> D, E, F; ...).

    Parameters
    ----------
    genes : list of dict
        One dict per gene, each with the same keys ``plot_LD_scheme``
        takes: 'results', 'exons', 'snps', 'gene_name', 'panel_x0',
        'panel_x1', 'chromosome', 'gene_start', 'gene_end'.
    output_name : str
        Filename the combined figure is saved to.
    col_width, fig_height : float
        Figure size is (col_width * n_genes, fig_height) inches.

    Returns
    -------
    fig, axes_per_gene
        ``axes_per_gene`` is a list of (ax_top, ax_snp_or_None, ax_bottom)
        tuples, one per gene, in the same order as ``genes``.
    """
    n_genes = len(genes)

    def _has_r2(snp):
        r2 = snp.get('r_squared', None)
        return r2 is not None and not pd.isna(r2)

    # Work out ahead of time whether each gene needs a variant panel,
    # so each column's gridspec rows can be set up correctly.
    needs_variant_panel = [
        any(not _has_r2(s) for s in gene['snps']) for gene in genes
    ]

    fig = plt.figure(figsize=(col_width * n_genes, fig_height))
    fig.patch.set_facecolor('white')

    gs = fig.add_gridspec(
        3, n_genes,
        height_ratios=[1.6, 0.55, 1],
        wspace=0.12, hspace=0.35,
    )

    axes_per_gene = []
    letter_cursor = 0

    for col, (gene, show_variant) in enumerate(zip(genes, needs_variant_panel)):
        ax_top = fig.add_subplot(gs[0, col])
        if show_variant:
            ax_snp = fig.add_subplot(gs[1, col], sharex=ax_top)
            ax_bottom = fig.add_subplot(gs[2, col], sharex=ax_top)
        else:
            ax_snp = None
            # No variant panel needed for this gene -> let its gene
            # scheme panel span the two lower rows instead of leaving
            # an empty gap.
            ax_bottom = fig.add_subplot(gs[1:3, col], sharex=ax_top)

        n_letters = _plot_gene_column(
            ax_top, ax_snp, ax_bottom,
            results=gene['results'], exons=gene['exons'], snps=gene['snps'],
            gene_name=gene['gene_name'],
            panel_x0=gene['panel_x0'], panel_x1=gene['panel_x1'],
            chromosome=gene['chromosome'],
            gene_start=gene['gene_start'], gene_end=gene['gene_end'],
            letter_offset=letter_cursor,
        )

        axes_per_gene.append((ax_top, ax_snp, ax_bottom))
        letter_cursor += n_letters

    fig.align_ylabels()
    plt.tight_layout()
    plt.savefig(output_name, dpi=300, bbox_inches='tight')

    return fig, axes_per_gene


# ══════════════════════════════════════════════════════════════════════════
# MULTI-GENE ENTRY POINT — vertically stacked genes in ONE figure
# ══════════════════════════════════════════════════════════════════════════
def plot_LD_scheme_stacked(genes, output_name='LD_plots_stacked.png',
                            fig_width=11, gene_height=7.2):
    """Plot several genes vertically, with each gene using three rows."""
    n_genes = len(genes)

    def _has_r2(snp):
        r2 = snp.get('r_squared', None)
        return r2 is not None and not pd.isna(r2)

    needs_variant_panel = [
        any(not _has_r2(s) for s in gene['snps']) for gene in genes
    ]

    fig = plt.figure(figsize=(fig_width, gene_height * n_genes))
    fig.patch.set_facecolor('white')

    rows_per_gene = 4
    gs = fig.add_gridspec(
        rows_per_gene * n_genes, 1,
        height_ratios=(
            [1.6, 0.55, 1, 0.4] * n_genes
        ),
        hspace=0.35,
    )

    axes_per_gene = []
    letter_cursor = 0

    for row, (gene, show_variant) in enumerate(zip(genes, needs_variant_panel)):
        row_start = row * rows_per_gene
        ax_top = fig.add_subplot(gs[row_start, 0])
        if show_variant:
            ax_snp = fig.add_subplot(gs[row_start + 1, 0], sharex=ax_top)
            ax_bottom = fig.add_subplot(gs[row_start + 2, 0], sharex=ax_top)
        else:
            ax_snp = None
            ax_bottom = fig.add_subplot(gs[row_start + 1:row_start + 3, 0],
                                        sharex=ax_top)

        n_letters = _plot_gene_column(
            ax_top, ax_snp, ax_bottom,
            results=gene['results'], exons=gene['exons'], snps=gene['snps'],
            gene_name=gene['gene_name'],
            panel_x0=gene['panel_x0'], panel_x1=gene['panel_x1'],
            chromosome=gene['chromosome'],
            gene_start=gene['gene_start'], gene_end=gene['gene_end'],
            letter_offset=letter_cursor,
        )

        axes_per_gene.append((ax_top, ax_snp, ax_bottom))
        letter_cursor += n_letters

    fig.align_ylabels()
    plt.tight_layout()
    plt.savefig(output_name, dpi=300, bbox_inches='tight')

    return fig, axes_per_gene


# ══════════════════════════════════════════════════════════════════════════
# BUILD REAL RESULTS DATAFRAMES + PLOT
# ══════════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    ANNOTATION = '../../../../../../../TEs_Aedes_aegypti/Annotation_VectorBase/VectorBase-54_AaegyptiLVP_AGWG.gff'

    # ── VGSC ────────────────────────────────────────────────────────────────
    chromosome_vgsc = 'AaegL5 3'
    gVs, gVe = 315926360, 316405639
    pVs, pVe = gVs - 1000, gVe + 1000
    try:
        vgsc_exons = get_exons(ANNOTATION, 'AAEL023266')
    except FileNotFoundError:
        print('GFF not found – using representative exon positions for VGSC.')
        vgsc_exons = [{'start_position': gVs + i * 13000,
                       'end_position':   gVs + i * 13000 + 2800} for i in range(37)]

    vgsc_snps = [
        {'TE_name': 'V1016I', 'start_position': 315983763, 'end_position': 315983764},
        {'TE_name': 'I915K',  'start_position': 315999297, 'end_position': 315999298},
        {'TE_name': 'S723T',  'start_position': 316014588, 'end_position': 316014589},
        {'TE_name': 'V410L',  'start_position': 316080722, 'end_position': 316080723},
    ]

    significant_TEs_ref_vgsc = pd.read_csv('../../../vgsc_gene/significant_r_squared_ref.csv')
    significant_TEs_ref_vgsc['type'] = 'Reference'
    significant_TEs_nonref_vgsc = pd.read_csv('../../../vgsc_gene/significant_r_squared_nonref.csv')
    significant_TEs_nonref_vgsc['type'] = 'Non-Reference'
    dfs_vgsc = [significant_TEs_ref_vgsc, significant_TEs_nonref_vgsc]

    results_vgsc = checking_genotypes(
        dfs_vgsc,
        path_df_ref='../../../vgsc_gene/genotype_info_Reference.csv',
        path_df_nonref='../../../vgsc_gene/genotype_info_Non-Reference.csv',
        snp_col_1='AaegL5_3_316080722_V410L_encoding_1',
        snp_col_2='AaegL5_3_316080722_V410L_encoding_2',
    )
    results_vgsc.to_csv('vgsc_results_with_LD_direction.csv', index=False)

    # ── RDL ─────────────────────────────────────────────────────────────────
    chromosome_rdl = 'AaegL5 2'
    gRs, gRe = 41628484, 41861946
    pRs, pRe = gRs - 1000, gRe + 1000
    try:
        rdl_exons = get_exons(ANNOTATION, 'AAEL008354')
    except FileNotFoundError:
        print('GFF not found – using representative exon positions for RDL.')
        rdl_exons = [{'start_position': gRs + i * 6800,
                      'end_position':   gRs + i * 6800 + 1900} for i in range(33)]

    rdl_snps = [
        {'TE_name': 'A301S', 'start_position': 41847790, 'end_position': 41847791},
    ]

    significant_TEs_ref_rdl = pd.read_csv('../../../rdl_gene/significant_r_squared_ref.csv')
    significant_TEs_ref_rdl['type'] = 'Reference'
    significant_TEs_nonref_rdl = pd.read_csv('../../../rdl_gene/significant_r_squared_nonref.csv')
    significant_TEs_nonref_rdl['type'] = 'Non-Reference'
    dfs_rdl = [significant_TEs_ref_rdl, significant_TEs_nonref_rdl]

    results_rdl = checking_genotypes(
        dfs_rdl,
        path_df_ref='../../../rdl_gene/genotype_info_Reference.csv',
        path_df_nonref='../../../rdl_gene/genotype_info_Non-Reference.csv',
        snp_col_1='AaegL5_2_41847790_A301S_encoding_1',
        snp_col_2='AaegL5_2_41847790_A301S_encoding_2',
    )
    results_rdl.to_csv('rdl_results_with_LD_direction.csv', index=False)

    # ── Plot both genes side by side in ONE figure ───────────────────────────
    # VGSC -> panels A, B, C ; RDL -> panels D, E, F
    plot_LD_scheme_combined(
        [
            dict(results=results_vgsc, exons=vgsc_exons, snps=vgsc_snps,
                 gene_name='VGSC', panel_x0=pVs, panel_x1=pVe,
                 chromosome=chromosome_vgsc, gene_start=gVs, gene_end=gVe),
            dict(results=results_rdl, exons=rdl_exons, snps=rdl_snps,
                 gene_name='RDL', panel_x0=pRs, panel_x1=pRe,
                 chromosome=chromosome_rdl, gene_start=gRs, gene_end=gRe),
        ],
        output_name='LD_plots_VGSC_RDL_combined.png',
    )

    # ── Plot both genes stacked, with RDL below VGSC ───────────────────────
    plot_LD_scheme_stacked(
        [
            dict(results=results_vgsc, exons=vgsc_exons, snps=vgsc_snps,
                 gene_name='VGSC', panel_x0=pVs, panel_x1=pVe,
                 chromosome=chromosome_vgsc, gene_start=gVs, gene_end=gVe),
            dict(results=results_rdl, exons=rdl_exons, snps=rdl_snps,
                 gene_name='RDL', panel_x0=pRs, panel_x1=pRe,
                 chromosome=chromosome_rdl, gene_start=gRs, gene_end=gRe),
        ],
        output_name='LD_plots_VGSC_RDL_stacked.png',
    )