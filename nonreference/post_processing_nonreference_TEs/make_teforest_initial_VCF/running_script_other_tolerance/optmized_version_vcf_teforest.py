import pandas as pd
import os
import re
import numpy as np
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

# ──────────────────────────────────────────────
# 1. LOAD TE NAME DICTIONARY
# ──────────────────────────────────────────────
def updating_TE_names():
    dict_names = (
        "/work/users/g/a/gabivla/lab/SV_mosquitos/TEs_Aedes_aegypti/"
        "running_pipeline_new_library/consensus_TEs_dict/"
        "consensus_file_dictionary_TEforest.txt"
    )
    name_dict = {}
    with open(dict_names, "r") as f:
        for line in f:
            if line.startswith("Original:"):
                parts = line.strip().split("|")
                original_name = parts[0].split(":")[1].strip()
                updated_name  = parts[1].split(":")[1].strip()
                name_dict[original_name] = updated_name
    return name_dict


# ──────────────────────────────────────────────
# 2. FAST WITHIN-SAMPLE FILE READER
# ──────────────────────────────────────────────
def load_one_file(path, genome_code, chromosome_variant, dict_TE_names):
    """Read one BED file, filter, rename, return a tidy DataFrame."""
    cols = ["chromosome", "position_start", "position_end", "info", "score", "strand"]
    #print(f"Reading: {path}")
    df = pd.read_csv(path, sep="\t", header=None, skiprows=1, names=cols,
                     usecols=[0, 1, 2, 3, 4, 5],
                     dtype={"chromosome": "category",
                            "score": "float32",
                            "position_start": "int32",
                            "position_end": "int32"})

    #print(f"  [{genome_code}] shape after read: {df.shape}")
    #print(f"  [{genome_code}] chromosomes found: {df['chromosome'].unique().tolist()}")

    df = df[df["chromosome"] == chromosome_variant]
    #print(f"  [{genome_code}] shape after chrom filter: {df.shape}")
    if df.empty:
        #print(f"  [{genome_code}] DROPPED: empty after chromosome filter")
        return None

    df = df[df["chromosome"].str.startswith("AaegL5_")]
    if df.empty:
        #print(f"  [{genome_code}] DROPPED: empty after AaegL5_ filter")
        return None

    info_split = df["info"].str.split("|", n=2, expand=True)
    df["TE_name"] = info_split[0].map(dict_TE_names).fillna(info_split[0])
    df["TE_type"] = info_split[1]
    #print(f"  [{genome_code}] TE_types found: {df['TE_type'].unique().tolist()}")

    df = df[df["TE_type"] == "non-reference"].drop(columns=["info", "TE_type"])
    #print(f"  [{genome_code}] shape after non-reference filter: {df.shape}")
    if df.empty:
        #print(f"  [{genome_code}] DROPPED: empty after non-reference filter")
        return None

    df = df.sort_values(["chromosome", "position_start", "position_end"]).reset_index(drop=True)
    df[genome_code] = "1/1"
    print(f"  [{genome_code}] SUCCESS: {len(df)} rows")
    return df


# ──────────────────────────────────────────────
# 3. FAST CROSS-INDIVIDUAL MERGE
#    Replace n-1 outer merges with concat + groupby
# ──────────────────────────────────────────────
KEY_COLS = ["chromosome", "position_start", "position_end", "TE_name", "score", "strand"]

def build_merged_df(list_of_dfs):
    """
    Stack all per-sample DataFrames, pivot genotype columns.
    Much faster than chained outer-merges.
    """
    # collect all sample column names
    sample_cols = []
    for df in list_of_dfs:
        extra = [c for c in df.columns if c not in KEY_COLS]
        sample_cols.extend(extra)

    stacked = pd.concat(list_of_dfs, ignore_index=True)

    # melt to long form: one row per (locus, sample)
    stacked_long = stacked.melt(id_vars=KEY_COLS, value_vars=sample_cols,
                                var_name="sample", value_name="genotype")
    stacked_long = stacked_long.dropna(subset=["genotype"])

    # pivot back to wide form; fill missing with 0/0
    df_merged = stacked_long.pivot_table(
        index=KEY_COLS, columns="sample", values="genotype",
        aggfunc="first"
    ).reset_index()
    df_merged.columns.name = None
    df_merged[sample_cols] = df_merged[sample_cols].fillna("0/0")

    return df_merged


# ──────────────────────────────────────────────
# 4. FAST WITHIN-TOLERANCE MERGING
#    Replace O(n²) nested loop with sort + iterative scan (O(n log n))
# ──────────────────────────────────────────────
def merging_individuals_insertions(df_merged, TOLERANCE):
    """
    For each (chromosome, TE_name, strand, score) group, merge rows whose
    position_start and position_end are within TOLERANCE of each other.

    Runs in O(n log n) via sort + single-pass scan instead of O(n²).
    Genotype rule: if either row is 1/1 → 1/1, else 0/0.
    """
    key_cols   = ["chromosome", "TE_name", "strand", "score"]
    geno_cols  = [c for c in df_merged.columns if c not in
                  ["chromosome", "position_start", "position_end",
                   "TE_name", "score", "strand"]]

    df_merged  = df_merged.sort_values(key_cols + ["position_start", "position_end"]).reset_index(drop=True)

    keep_mask  = np.ones(len(df_merged), dtype=bool)
    # work on numpy arrays for speed
    starts  = df_merged["position_start"].to_numpy(dtype=np.int64)
    ends    = df_merged["position_end"].to_numpy(dtype=np.int64)
    geno_np = df_merged[geno_cols].to_numpy()          # object array of "0/0" / "1/1"

    # build a group-id array so we only compare within the same group
    group_ids = df_merged.groupby(key_cols, sort=False).ngroup().to_numpy()

    i = 0
    n = len(df_merged)
    while i < n:
        if not keep_mask[i]:
            i += 1
            continue
        j = i + 1
        while j < n and group_ids[j] == group_ids[i]:
            if keep_mask[j]:
                if (abs(starts[i] - starts[j]) <= TOLERANCE and
                        abs(ends[i]   - ends[j])   <= TOLERANCE):
                    # merge j into i
                    starts[i] = min(starts[i], starts[j])
                    ends[i]   = max(ends[i],   ends[j])
                    # genotype: 1/1 wins
                    geno_np[i] = np.where(
                        (geno_np[i] == "1/1") | (geno_np[j] == "1/1"),
                        "1/1", "0/0"
                    )
                    keep_mask[j] = False
            j += 1
        i += 1

    df_merged = df_merged[keep_mask].copy()
    df_merged["position_start"] = starts[keep_mask]
    df_merged["position_end"]   = ends[keep_mask]
    df_merged[geno_cols] = geno_np[keep_mask]
    df_merged.reset_index(drop=True, inplace=True)
    return df_merged


# ──────────────────────────────────────────────
# 5. MAIN
# ──────────────────────────────────────────────
def main():
    chromosome_variant = sys.argv[1]
    TOLERANCE = int(sys.argv[2])

    dict_TE_names = updating_TE_names()

    bed_dir   = "../../../results_TEforest/"
    bed_files = [
        f for f in os.listdir(bed_dir)
        if f.endswith(".bed") and not f.startswith("JB")
    ]
    print(f'number of genomes to process: {len(bed_files)}')
    # ── parallel file loading ──────────────────
    tasks = []
    for fname in bed_files:
        m = re.search(r"(.+)_TEforest_bps_nonredundant\.bed", fname)
        if not m:
            continue
        genome_code = m.group(1)
        full_path   = os.path.join(bed_dir, fname)
        tasks.append((full_path, genome_code))
    list_of_dfs = []
    # ThreadPoolExecutor is sufficient because pandas I/O releases the GIL
    with ThreadPoolExecutor(max_workers=min(16, len(tasks))) as pool:
        futures = {
            pool.submit(load_one_file, path, code, chromosome_variant, dict_TE_names): code
            for path, code in tasks
        }
        for fut in as_completed(futures):
            result = fut.result()
            if result is not None:
                list_of_dfs.append(result)

    if not list_of_dfs:
        print("No data found for chromosome:", chromosome_variant)
        return

    # ── merge individuals ──────────────────────
    df_merged = build_merged_df(list_of_dfs)
    df_merged = df_merged.sort_values(["chromosome", "position_start", "position_end"]).reset_index(drop=True)

    # ── collapse nearby insertions ─────────────
    df_merged = merging_individuals_insertions(df_merged, TOLERANCE)

    # ── save ───────────────────────────────────
    out = f"VCF_TEforest_nonreference_{chromosome_variant}_op_{TOLERANCE}.csv"
    df_merged.to_csv(out, index=False)
    print(f"Saved {len(df_merged)} rows → {out}")

if __name__ == "__main__":
    main()