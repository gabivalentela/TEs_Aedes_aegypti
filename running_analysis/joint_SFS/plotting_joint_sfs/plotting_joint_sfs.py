import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import itertools
import allel


def get_count_derived(df, country_name, base_cols):
    genotype_cols = [col for col in df.columns if col not in base_cols]
    df[f'count_1/1_{country_name}'] = df[genotype_cols].apply(lambda row: (row == '1/1').sum(), axis=1)
    df[f'count_0/0_{country_name}'] = df[genotype_cols].apply(lambda row: (row == '0/0').sum(), axis=1)
    df[f'count_0/1_{country_name}'] = df[genotype_cols].apply(lambda row: (row == '0/1').sum(), axis=1)
    df[f'count_missing_{country_name}'] = df[genotype_cols].apply(lambda row: (row == './.').sum(), axis=1)
    df[f'total_samples_{country_name}'] = (
        df[f'count_1/1_{country_name}']
        + df[f'count_0/0_{country_name}']
        + df[f'count_0/1_{country_name}']
        + df[f'count_missing_{country_name}']
    )
    df[f'count_derived_chromosomes_{country_name}'] = (df[f'count_1/1_{country_name}'] * 2) + df[f'count_0/1_{country_name}']
    return df


reference_path = '/work/users/g/a/gabivla/lab/SV_mosquitos/TE_paper_aedesaegytpi/reference_vcf/genotyping_reference_TEs/post_processing/merged_results_5_2_either_filtered.csv'
non_reference_path = '/work/users/g/a/gabivla/lab/SV_mosquitos/TE_paper_aedesaegytpi/nonreference_vcf/post_processing_nonreference_TEs/post_processing_overlap_df/max_min_updated_strategy/final_df_overlap_retroseq_teforest_nonreference_1000.csv'
df_reference = pd.read_csv(reference_path)
df_nonref = pd.read_csv(non_reference_path)

location_path = '/nas/longleaf/home/gabivla/lab/SV_mosquitos/location_information/mcclintock_output_locations.csv'
location_df = pd.read_csv(location_path, usecols=['Genome_code', 'country'])
location_dict = dict(zip(location_df['Genome_code'], location_df['country']))

base_cols_ref = ['chromosome', 'position_start', 'position_end', 'TE_name', 'strand', 'score']
base_cols_nonref = ['chromosome', 'position_start', 'position_end', 'TE_name', 'strand']

genome_cols_nonref = [col for col in df_nonref.columns if col not in base_cols_nonref]

rename_dict_nonref = {col: f"{col}_{location_dict.get(col, '')}" for col in genome_cols_nonref}
df_nonref = df_nonref.loc[:, ~df_nonref.columns.str.contains("Unnamed")]
df_nonref = df_nonref.rename(columns=rename_dict_nonref)

countries_nonreference = set()
for original_col in genome_cols_nonref:
    new_col = rename_dict_nonref[original_col]
    country = new_col.split('_')[-1]
    countries_nonreference.add(country)

genome_cols_ref = [col for col in df_reference.columns if col not in base_cols_ref]

# update the genome columns names
for genome in genome_cols_ref:
    if genome.startswith('SRR'):
        genome_updated = genome.replace('_1', '')
    elif genome.startswith('FEMALE') or genome.startswith('MALE'):
        genome_updated = genome.replace('_L001_L002_R1_001', '')
    # update the column name in the dataframe
    df_reference = df_reference.rename(columns={genome: genome_updated})

genome_cols_ref = [col for col in df_reference.columns if col not in base_cols_ref]

rename_dict_ref = {col: f"{col}_{location_dict.get(col, '')}" for col in genome_cols_ref}
df_reference = df_reference.loc[:, ~df_reference.columns.str.contains("Unnamed")]
df_reference = df_reference.rename(columns=rename_dict_ref)

countries_reference = set()
for original_col in genome_cols_ref:
    new_col = rename_dict_ref[original_col]
    country = new_col.split('_')[-1]
    countries_reference.add(country)

# running function for reference and nonreference TEs
# NOTE the added 'ploidy_scale' flag:
#   - reference TEs have true diploid genotype calls (0/0, 0/1, 1/1) -> derived count
#     scale runs 0..2N (chromosome count), so ploidy_scale = 2
#   - non-reference TEs (RetroSeq/TEForest) only ever call 0/0 / 1/1 (no heterozygotes
#     are called), so the "chromosome count" can only ever land on even numbers. That
#     structurally empties every odd bin of the joint SFS, which is what produced the
#     white/blank striping. Building the SFS on the individual-count scale (0..N)
#     instead removes those impossible bins entirely.
dict_dfs = {
    'reference': {'df': df_reference, 'base_columns': base_cols_ref, 'countries': countries_reference, 'ploidy_scale': 2},
    'nonreference': {'df': df_nonref, 'base_columns': base_cols_nonref, 'countries': countries_nonreference, 'ploidy_scale': 1},
}

population_sizes = {'Brazil': 19, 'Gabon': 13, 'Colombia': 24, 'USA': 28, 'Senegal': 20, 'Kenya': 19}

for key, value in dict_dfs.items():
    print(f'running it for {key} TEs')
    df = value['df']
    base_cols = value['base_columns']
    countries = value['countries']
    ploidy_scale = value['ploidy_scale']

    df_countries = {}

    for i, country in enumerate(countries):
        country_cols = [col for col in df.columns
                    if col not in base_cols and col.endswith('_' + country)]

        if not country_cols:
            continue

        country_df = df[base_cols + country_cols]
        country_df = country_df.dropna(subset=country_cols, how='all')

        print(f"\n--- Data for {country} ---")

        frequency_TEs = get_count_derived(country_df, country, base_cols)

        if ploidy_scale == 1:
            # non-reference: no heterozygotes are ever called, so collapse the
            # derived count onto the individual scale (0..N) instead of the
            # chromosome scale (0..2N) to avoid structurally-empty odd bins.
            frequency_TEs[f'count_derived_chromosomes_{country}'] = frequency_TEs[f'count_1/1_{country}']

        df_countries[country] = frequency_TEs

    country_list = list(df_countries.keys())
    country_pairs = list(itertools.combinations(country_list, 2))

    num_pairs = len(country_pairs)
    cols = min(5, num_pairs)
    rows = (num_pairs + cols - 1) // cols

    fig, axes = plt.subplots(rows, cols, figsize=(4 * cols, 4 * rows), constrained_layout=True)
    axes = np.array(axes).flatten()

    for i, (country_1, country_2) in enumerate(country_pairs):
        df_country_1 = df_countries[country_1]
        df_country_2 = df_countries[country_2]

        df_country_1['TE_information'] = df_country_1['chromosome'].astype(str) + '_' + df_country_1['position_start'].astype(str) + '_' + df_country_1['position_end'].astype(str) + '_' + df_country_1['TE_name'] + '_' + df_country_1['strand']
        df_country_2['TE_information'] = df_country_2['chromosome'].astype(str) + '_' + df_country_2['position_start'].astype(str) + '_' + df_country_2['position_end'].astype(str) + '_' + df_country_2['TE_name'] + '_' + df_country_2['strand']

        df_country_1 = df_country_1[['TE_information', f'count_derived_chromosomes_{country_1}']]
        df_country_2 = df_country_2[['TE_information', f'count_derived_chromosomes_{country_2}']]

        merged_df = df_country_1.merge(df_country_2, on="TE_information", how="outer", suffixes=("_country1", "_country2"))
        merged_df = merged_df[(merged_df[f'count_derived_chromosomes_{country_1}'] != 0) | (merged_df[f'count_derived_chromosomes_{country_2}'] != 0)]
        dac1 = merged_df[f'count_derived_chromosomes_{country_1}'].to_list()
        dac2 = merged_df[f'count_derived_chromosomes_{country_2}'].to_list()

        # n1/n2 must match the scale used for dac1/dac2: doubled population size
        # (chromosome count) for reference, raw population size (individual
        # count) for nonreference.
        n1 = population_sizes[country_1] * ploidy_scale
        n2 = population_sizes[country_2] * ploidy_scale

        # calculate joint SFS using allel
        joint_sfs = allel.joint_sfs(dac1, dac2, n1, n2)
        print(key)
        print(f'Joint SFS for {country_1} vs {country_2}:')
        print(joint_sfs)
        # plotting the joint SFS
        allel.plot_joint_sfs(joint_sfs, ax=axes[i], imshow_kwargs={'cmap': 'viridis'})
        # add the colorbar referencing the axes of the current subplot
        cbar = plt.colorbar(axes[i].images[0], ax=axes[i])
        cbar.set_label('Number of TEs')
        # setting the title and labels
        axes[i].set_title(f'Joint SFS: {country_1} vs {country_2}')
        axis_unit = 'chromosomes' if ploidy_scale == 2 else 'individuals'
        axes[i].set_xlabel(f'Derived TE count ({axis_unit}) in {country_1}')
        axes[i].set_ylabel(f'Derived TE count ({axis_unit}) in {country_2}')

    plt.savefig(f'joint_sfs_{key}.png', dpi=300)
    plt.close()