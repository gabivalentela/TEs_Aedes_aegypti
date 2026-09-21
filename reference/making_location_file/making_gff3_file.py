#!/usr/bin/env python3
"""
elem_sorted_to_gff3.py  –  Convert OneCodeToFindThemAll .elem_sorted.csv files
                            to a GFF3 annotation file and a matching taxonomy TSV.

Usage:
    python elem_sorted_to_gff3.py <tolerance_bp>

The taxonomy TSV is derived exclusively from records that survive into the
final GFF3, so the two files are guaranteed to be consistent.
"""

import os
import sys
from collections import defaultdict

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if len(sys.argv) < 2:
    print("Usage: python elem_sorted_to_gff3.py <tolerance_bp>")
    sys.exit(1)

tolerance   = int(sys.argv[1])
input_dir = ("/users/g/a/gabivla/New_library_test/onecodetofindthemall/Onecodetofindthemall/final_run_reference/results_onecodetofindthemall/chromosome")
output_dir  = "./results"
output_gff  = os.path.join(output_dir, f"all_chromosomes_TEs_{tolerance}.gff3")
output_tax  = os.path.join(output_dir, f"all_chromosomes_taxonomy_{tolerance}.tsv")

os.makedirs(output_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def parse_elem_sorted(filepath, preview=False):
    """
    Parse a single .elem_sorted.csv file produced by OneCodeToFindThemAll.

    Lines of interest start with '###'.  Returns a list of record dicts.
    If *preview* is True the first three records are printed to stdout.
    """
    records   = []
    previewed = 0

    with open(filepath) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line.startswith("###"):
                continue

            parts = line[3:].split("\t")

            if len(parts) < 16:
                print(f"  [WARN] Skipping short line ({len(parts)} fields): {line}")
                continue

            try:
                score_raw     = parts[0].strip()
                pct_div       = parts[1].strip()
                query         = parts[4].strip()
                beg           = int(parts[5].strip())
                end           = int(parts[6].strip())
                sense         = parts[8].strip()
                element       = parts[9].strip()
                family        = parts[10].strip()
                te_id         = parts[14].strip()
                num_assembled = parts[15].strip()
                pct_ref       = parts[16].strip() if len(parts) > 16 else "."
            except (ValueError, IndexError) as exc:
                print(f"  [WARN] Parse error ({exc}): {line}")
                continue

            strand    = "-" if sense == "C" else "+"
            gff_score = pct_div if pct_div not in ("", ".") else "."

            try:
                pct_div_num = float(pct_div)
            except (ValueError, TypeError):
                pct_div_num = float("inf")

            record = {
                "seqname":       query,
                "source":        "OneCodeToFindThemAll",
                "feature":       family,
                "start":         beg,
                "end":           end,
                "score":         gff_score,
                "strand":        strand,
                "frame":         ".",
                "element":       element,
                "family":        family,
                "te_id":         te_id,
                "num_assembled": num_assembled,
                "pct_div":       pct_div,
                "pct_div_num":   pct_div_num,
                "pct_ref":       pct_ref,
                "score_raw":     score_raw,
            }
            records.append(record)

            if preview and previewed < 3:
                print(f"\n  --- Preview record {previewed + 1} ---")
                print(f"  Raw line  : {line}")
                print(f"  seqname   : {query}")
                print(f"  feature   : {family}")
                print(f"  start/end : {beg} / {end}")
                print(f"  strand    : {strand}   score: {gff_score}")
                print(f"  element   : {element}")
                print(f"  te_id     : {te_id}   num_assembled: {num_assembled}   pct_ref: {pct_ref}")
                previewed += 1

    return records


# ---------------------------------------------------------------------------
# Deduplication
# ---------------------------------------------------------------------------

def deduplicate(records, threshold):
    """
    Remove duplicate records in two passes:

    Pass 1 – Exact duplicates
        Same (seqname, start, end, element) → keep the one with lowest pct_div.

    Pass 2 – Near-identical duplicates
        Same (seqname, element), start AND end both within *threshold* bp of a
        cluster representative → keep the one with lowest pct_div.
        Uses a sorted-sweep so no records are skipped.
    """
    n_input = len(records)

    # --- Pass 1: exact duplicates ---
    records.sort(key=lambda r: r["pct_div_num"])
    seen_exact: dict = {}
    for r in records:
        key = (r["seqname"], r["start"], r["end"], r["element"])
        if key not in seen_exact:
            seen_exact[key] = r
    deduped = list(seen_exact.values())
    n_exact_removed = n_input - len(deduped)
    print(f"  [DEDUP] Exact duplicates removed : {n_exact_removed}")

    # --- Pass 2: near-identical duplicates ---
    groups: dict = defaultdict(list)
    for r in deduped:
        groups[(r["seqname"], r["element"])].append(r)

    final            = []
    n_near_removed   = 0

    for group in groups.values():
        group.sort(key=lambda r: (r["start"], r["end"]))

        # Sweep: cluster records whose start/end are both within threshold of
        # the cluster's first (representative) record.
        clusters      = []
        current       = [group[0]]

        for r in group[1:]:
            rep = current[0]
            if (abs(r["start"] - rep["start"]) <= threshold and
                    abs(r["end"]   - rep["end"])   <= threshold):
                current.append(r)
            else:
                clusters.append(current)
                current = [r]
        clusters.append(current)

        for cluster in clusters:
            best = min(cluster, key=lambda r: r["pct_div_num"])
            n_near_removed += len(cluster) - 1
            final.append(best)

    print(f"  [DEDUP] Near-identical removed (±{threshold} bp) : {n_near_removed}")
    print(f"  [DEDUP] Records remaining                        : {len(final)}")
    return final


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def write_gff3(records, filepath):
    """Write records to a GFF3 file. Returns the set of Name= values written."""
    gff_names: set = set()

    with open(filepath, "w") as out:
        out.write("##gff-version 3\n")
        for r in records:
            gff_id = f"{r['element']}_{r['seqname']}_{r['start']}_{r['end']}"
            gff_names.add(r["element"])

            attrs = (
                f"ID={gff_id};"
                f"Name={r['element']};"
                f"Family={r['family']};"
                f"Fragment_IDs={r['te_id']};"
                f"Num_Assembled={r['num_assembled']};"
                f"Pct_Div={r['pct_div']};"
                f"Pct_Ref={r['pct_ref']};"
                f"Score={r['score_raw']}"
            )

            out.write(
                f"{r['seqname']}\t"
                f"{r['source']}\t"
                f"{r['feature']}\t"
                f"{r['start']}\t"
                f"{r['end']}\t"
                f"{r['score']}\t"
                f"{r['strand']}\t"
                f"{r['frame']}\t"
                f"{attrs}\n"
            )

    return gff_names


def write_taxonomy(records, filepath):
    with open(filepath, "w") as out:
        for r in records:
            gff_id = f"{r['element']}_{r['seqname']}_{r['start']}_{r['end']}"
            out.write(f"{gff_id}\t{r['element']}\n")  # <-- r['element'] not r['family']

    return len(records)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    all_records: list = []
    first_file        = True

    for fname in sorted(os.listdir(input_dir)):
        if not fname.endswith(".elem_sorted.csv"):
            continue

        fpath = os.path.join(input_dir, fname)
        print(f"\nProcessing: {fname}")
        recs = parse_elem_sorted(fpath, preview=first_file)
        first_file = False
        print(f"  -> {len(recs)} assembled TE records found")
        all_records.extend(recs)

    print(f"\nTotal assembled TE records (all chromosomes): {len(all_records)}")

    # --- Deduplication ---
    print(f"\nRunning deduplication (tolerance = {tolerance} bp) ...")
    all_records = deduplicate(all_records, threshold=tolerance)

    # --- GFF3 ---
    gff_names = write_gff3(all_records, output_gff)
    print(f"\nGFF3 written to    : {output_gff}")
    print(f"Records written    : {len(all_records)}")
    print(f"Unique Name= values: {len(gff_names)}")

    # --- Taxonomy TSV ---
    # One row per GFF3 record using the full ID= value (element_seqname_start_end)
    # so every taxonomy entry matches exactly one GFF3 ID= field.
    n_written = write_taxonomy(all_records, output_tax)
    print(f"\nTaxonomy TSV written to : {output_tax}")
    print(f"Rows written            : {n_written}")
    print("  [OK] Taxonomy IDs match GFF3 ID= fields by construction")


if __name__ == "__main__":
    main()