import os
import gzip
from itertools import combinations

import numpy as np
import pandas as pd
import allel
import matplotlib.pyplot as plt

GT_MAP = {"0/0": (0, 0), "0/1": (0, 1), "1/0": (0, 1), "1/1": (1, 1), "0|0": (0, 0), "0|1": (0, 1), "1|0": (0, 1), "1|1": (1, 1)}

def genotype_df_to_array(df, genotype_cols):
    n_variants = len(df)
    n_samples = len(genotype_cols)
    arr = np.full((n_variants, n_samples, 2), -1, dtype=np.int8)

    for j, col in enumerate(genotype_cols):
        col_vals = df[col].values
        for code, (a1, a2) in GT_MAP.items():
            mask = col_vals == code
            arr[mask, j, 0] = a1
            arr[mask, j, 1] = a2

    return allel.GenotypeArray(arr)


def calculate_FST_allel(df, cols_p1, cols_p2, meta_cols, method="hudson"):
    genotype_cols = cols_p1 + cols_p2
    g = genotype_df_to_array(df, genotype_cols)

    subpops = [
        list(range(0, len(cols_p1))),
        list(range(len(cols_p1), len(cols_p1) + len(cols_p2))),
    ]

    if method == "hudson":
        ac1, ac2 = g.count_alleles_subpops({"p1": subpops[0], "p2": subpops[1]}).values()
        num, den = allel.hudson_fst(ac1, ac2)
        fst = np.full_like(num, np.nan, dtype=np.float64)
        np.divide(num, den, out=fst, where=den > 0)
    else:
        raise ValueError("method must be 'hudson'")

    out = df[meta_cols].copy()
    out["FST"] = fst
    return out


def calculate_FST_genomewide(df, cols_p1, cols_p2, method="hudson"):
    genotype_cols = cols_p1 + cols_p2
    g = genotype_df_to_array(df, genotype_cols)

    subpops = [
        list(range(0, len(cols_p1))),
        list(range(len(cols_p1), len(cols_p1) + len(cols_p2))),
    ]

    if method == "hudson":
        ac1, ac2 = g.count_alleles_subpops({"p1": subpops[0], "p2": subpops[1]}).values()
        num, den = allel.hudson_fst(ac1, ac2)
        return np.nansum(num) / np.nansum(den)
    else:
        raise ValueError("method must be 'hudson'")


def normalize_genome_name(name):
    """Strip dataset-specific suffixes so sample IDs match the location file."""
    if name.startswith("SRR"):
        return name.replace("_1", "")
    elif name.startswith("FEMALE") or name.startswith("MALE"):
        return name.replace("_L001_L002_R1_001", "")
    return name


def load_location_dict(location_path):
    location_df = pd.read_csv(location_path, usecols=["Genome_code", "country"])
    return dict(zip(location_df["Genome_code"], location_df["country"]))


def add_country_suffix(df, genome_cols, location_dict):
    """Rename each genome column to '{genome}_{country}' and report countries seen."""
    df = df.loc[:, ~df.columns.str.contains("Unnamed")].copy()
    rename_dict = {col: f"{col}_{location_dict.get(col, 'Unknown')}" for col in genome_cols}
    df = df.rename(columns=rename_dict)
    countries = {rename_dict[c].split("_")[-1] for c in genome_cols}
    return df, countries


def read_te_reference(csv_path):
    """Read the wide TE-genotype CSV (one row per TE call, one col per sample)."""
    df = pd.read_csv(csv_path)
    df = df.drop(columns=["score"])

    meta_cols = ["chromosome", "position_start", "position_end", "strand", "TE_name"]
    genome_cols = [c for c in df.columns if c not in meta_cols]

    for genome in genome_cols:
        df = df.rename(columns={genome: normalize_genome_name(genome)})
    genome_cols = [c for c in df.columns if c not in meta_cols]

    return df, meta_cols, genome_cols


def read_snp_vcf(vcf_path):
    """Read a gzipped VCF into the same wide (meta + genotype) shape as the TE CSV."""
    with gzip.open(vcf_path, "rt") as f:
        header = None
        for line in f:
            if line.startswith("#CHROM"):
                header = line.strip("#").strip().split("\t")
                break
        df = pd.read_csv(f, sep="\t", comment="#", header=None, names=header)

    drop_cols = ["ID", "REF", "ALT", "QUAL", "FILTER", "FORMAT", "INFO"]
    df = df.drop(columns=[c for c in drop_cols if c in df.columns])

    meta_cols = ["CHROM", "POS"]
    sample_cols = [c for c in df.columns if c not in meta_cols]

    # keep only the GT portion of each sample cell (e.g. "0/1:30:..." -> "0/1")
    df[sample_cols] = df[sample_cols].map(
        lambda x: x.split(":")[0] if isinstance(x, str) else x
    )

    for col in sample_cols:
        df = df.rename(columns={col: normalize_genome_name(col)})
    genome_cols = [c for c in df.columns if c not in meta_cols]

    return df, meta_cols, genome_cols

def run_fst_analysis(
    df,
    meta_cols,
    genome_cols,
    location_dict,
    output_dir,
    label,
    exclude_countries=None,
    make_histogram=False,
):
    os.makedirs(output_dir, exist_ok=True)

    df, countries = add_country_suffix(df, genome_cols, location_dict)
    if exclude_countries:
        countries = countries - set(exclude_countries)

    country_cols = {}
    for country in countries:
        cols = [c for c in df.columns if c not in meta_cols and c.endswith("_" + country)]
        if cols:
            country_cols[country] = cols

    country_pairs = list(combinations(country_cols.keys(), 2))

    fig, axes = None, None
    if make_histogram and country_pairs:
        ncols = 4
        nrows = int(np.ceil(len(country_pairs) / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(10 * ncols, 7 * nrows))
        axes = np.atleast_2d(axes).reshape(-1)

    summary_rows = []
    plot_idx = 0

    for country1, country2 in country_pairs:
        cols1, cols2 = country_cols[country1], country_cols[country2]
        sub_df = df.dropna(subset=cols1 + cols2, how="all").copy()

        if sub_df.empty:
            continue

        print(f"\n--- {label} FST: {country1} vs {country2} ({len(sub_df)} candidate sites) ---")

        fst_result = calculate_FST_allel(
            sub_df, cols_p1=cols1, cols_p2=cols2, meta_cols=meta_cols, method="hudson"
        )
        fst_result = fst_result.dropna(subset=["FST"])

        genomewide_fst = calculate_FST_genomewide(sub_df, cols1, cols2, method="hudson")
        print(f"Genome-wide weighted FST ({country1} vs {country2}): {genomewide_fst:.4f}")

        summary_rows.append(
            {
                "country1": country1,
                "country2": country2,
                "genome_wide_FST": genomewide_fst,
                "n_sites": len(fst_result),
            }
        )

        out_path = os.path.join(output_dir, f"FST_values_{country1}_{country2}_{label}.csv")
        fst_result.to_csv(out_path, index=False)

        if make_histogram and axes is not None:
            ax = axes[plot_idx]
            fst_values = fst_result["FST"].values
            counts, bins = np.histogram(fst_values, bins=20, range=(0, 1))
            total = counts.sum()
            percentages = (counts / total) * 100 if total > 0 else counts
            ax.bar(bins[:-1], percentages, width=np.diff(bins), edgecolor="black", alpha=0.7)
            ax.set_title(f"FST: {country1} vs {country2}")
            ax.set_xlabel("FST value")
            ax.set_ylabel("Percentage")
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 100)
            plot_idx += 1

    if make_histogram and axes is not None:
        for j in range(plot_idx, len(axes)):
            axes[j].axis("off")
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"fst_histograms_{label}.png"))
        plt.close(fig)

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(os.path.join(output_dir, f"FST_genomewide_summary_{label}.csv"), index=False)
    return summary_df

def main():
    location_path = (
        "/nas/longleaf/home/gabivla/lab/SV_mosquitos/location_information/"
        "mcclintock_output_locations.csv"
    )
    location_dict = load_location_dict(location_path)

    snp_vcf_path = (
        "/work/users/g/a/gabivla/lab/SV_mosquitos/neutral_data_vcf/"
        "updated_AaegL5_full_filtered_snps_110122_rep_map_100kbgenes_"
        "centromere_filtered_sites_mask.vcf.gz"
    )
    snp_df, snp_meta_cols, snp_genome_cols = read_snp_vcf(snp_vcf_path)
    run_fst_analysis(
        snp_df,
        meta_cols=snp_meta_cols,
        genome_cols=snp_genome_cols,
        location_dict=location_dict,
        output_dir="./FST_analysis/SNPs",
        label="SNPs",
        exclude_countries=["Cali"],
        make_histogram=True,
    )

if __name__ == "__main__":
    main()