import pandas as pd
import numpy as np
import allel
from itertools import combinations

GT_MAP = {"0/0": (0, 0), "0/1": (0, 1), "1/0": (0, 1), "1/1": (1, 1)}

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
        # count alleles in each subpopulation
        ac1, ac2 = g.count_alleles_subpops({"p1": subpops[0], "p2": subpops[1]}).values()
        # numerador and denominador and then can calculate per window or variant individually from this
        num, den = allel.hudson_fst(ac1, ac2)
        # if the denominador is zero, FST is zero (heterozygosity in the whole population is zero, meaning no variation)
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
        raise ValueError("method must be 'hudson' or 'weir_cockerham'")

def adding_location(df, genome_cols):
    location_path = (
        "/nas/longleaf/home/gabivla/lab/SV_mosquitos/location_information/"
        "mcclintock_output_locations.csv"
    )
    location_df = pd.read_csv(location_path, usecols=["Genome_code", "country"])
    location_dict = dict(zip(location_df["Genome_code"], location_df["country"]))

    rename_dict = {col: f"{col}_{location_dict.get(col, '')}" for col in genome_cols}
    df = df.loc[:, ~df.columns.str.contains("Unnamed")]
    df = df.rename(columns=rename_dict)

    countries = set()
    for original_col in genome_cols:
        new_col = rename_dict[original_col]
        country = new_col.split("_")[-1]
        countries.add(country)

    return df, countries

def main():
    reference_path = (
        "/work/users/g/a/gabivla/lab/SV_mosquitos/TE_paper_aedesaegytpi/reference_vcf/genotyping_reference_TEs/post_processing/merged_results_5_2_either_filtered.csv"
    )

    df = pd.read_csv(reference_path)
    df = df.drop(columns=["score"])

    base_columns = ["chromosome", "position_start", "position_end", "strand", "TE_name"]
    genome_cols = [c for c in df.columns if c not in base_columns]

    for genome in genome_cols:
        if genome.startswith("SRR"):
            genome_updated = genome.replace("_1", "")
        elif genome.startswith("FEMALE") or genome.startswith("MALE"):
            genome_updated = genome.replace("_L001_L002_R1_001", "")
        else:
            genome_updated = genome
        df = df.rename(columns={genome: genome_updated})

    genome_cols = [c for c in df.columns if c not in base_columns]

    df, countries = adding_location(df, genome_cols)

    country_cols = {}
    for country in countries:
        cols = [c for c in df.columns if c not in base_columns and c.endswith("_" + country)]
        if cols:
            country_cols[country] = cols

    country_pairs = list(combinations(country_cols.keys(), 2))

    for country1, country2 in country_pairs:
        cols1 = country_cols[country1]
        cols2 = country_cols[country2]

        sub_df = df.dropna(subset=cols1 + cols2, how="all").copy()

        if sub_df.empty:
            continue

        print(f"\n--- FST: {country1} vs {country2} ({len(sub_df)} candidate sites) ---")

        fst_result = calculate_FST_allel(
            sub_df,
            cols_p1=cols1,
            cols_p2=cols2,
            meta_cols=base_columns,
            method="hudson"
        )
        # remove the ones that are NaN
        fst_result = fst_result.dropna(subset=["FST"])

        genomewide_fst = calculate_FST_genomewide(sub_df, cols1, cols2, method="hudson")
        print(f"Genome-wide weighted FST ({country1} vs {country2}): {genomewide_fst:.4f}")

        out_path = f"./FST_analysis/FST_values_{country1}_{country2}_ref.csv"
        fst_result.to_csv(out_path, index=False)

if __name__ == "__main__":
    main()