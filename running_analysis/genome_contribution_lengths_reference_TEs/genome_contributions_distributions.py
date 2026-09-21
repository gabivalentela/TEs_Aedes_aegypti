import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import gzip

# Colors matched to the violin plot
TE_COLORS = {
    'LINE':     '#E74C3C',   # red
    'SINE':     '#3498DB',   # blue
    'LTR':      '#2ECC71',   # green
    'MITE':     '#F1C40F',   # yellow
    'TIR':      '#FA8072',   # salmon
    'HELITRON': '#9B59B6',   # purple
    'PENELOPE': '#E67E22',   # orange
    'UNKNOWN':  '#888888',   # grey fallback
    'OTHER':    '#AAAAAA',   # light grey fallback
    'SNP':      '#000000',   # black for the SNP box
}

# sizes of chromosomes - NCBI
CHROM_SIZES = {
    'AaegL5_1': 310_827_022,
    'AaegL5_2': 474_425_716,
    'AaegL5_3': 409_777_670,
}

base_columns = ['chromosome', 'position_start', 'position_end', 'TE_name', 'strand', 'score']

def define_main_class(df):
    df['upper_te'] = df['TE_name'].str.upper()
    main_classes = ['MITE', 'TIR', 'LINE', 'SINE', 'LTR', 'HELITRON', 'PENELOPE', 'UNKNOWN']
    df['main_class'] = df['upper_te'].apply(
        lambda x: next((mc for mc in main_classes if mc in x), 'OTHER')
    )
    df.drop(columns=['upper_te'], inplace=True)
    return df

def calculate_lengths(df):
    df['length'] = df['position_end'] - df['position_start']
    return df

def plot_length_distributions(df):
    unique_classes = list(df['main_class'].unique())
    bins = np.histogram_bin_edges(df['length'], bins=50)
    fig, axes = plt.subplots(len(unique_classes), 1, figsize=(12, 3 * len(unique_classes)), sharex=True)
    axes = np.atleast_1d(axes)

    for axis, main_class in zip(axes, unique_classes):
        class_lengths = df[df['main_class'] == main_class]['length']
        color = TE_COLORS.get(main_class, '#AAAAAA')  # fallback to grey if class not in dict
        axis.hist(class_lengths, bins=bins, alpha=0.8, color=color)
        axis.set_ylabel('Frequency')
        axis.set_title(main_class)

    axes[-1].set_xlabel('Length')
    fig.suptitle('Length Distributions by Main Class')
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig('length_distributions.png')
    plt.close()

def compute_base_pair_contributions_TE(df):
    # dict of contributions per genome: {genome: DataFrame[main_class, percentage]}
    contributions_dict = {}
    # genotype columns define
    genotype_columns = [col for col in df.columns if col not in base_columns + ['main_class', 'length']]
    for col in genotype_columns:
        # subset with just that genome + base columns + main_class + length
        subset_df = df[base_columns + ['main_class', 'length'] + [col]]
        # keep only the ones that are 1/1 or 0/1 - present in that genome
        subset_df = subset_df[subset_df[col].isin(['1/1', '0/1', '1|1', '0|1', '1|0', '1/0'])]
        # calculate total length of base pairs for each main class
        contributions = subset_df.groupby('main_class')['length'].sum().reset_index()
        contributions = contributions.sort_values(by='length', ascending=False)
        # calculate the total length of base pairs in the genome
        total_length = sum(CHROM_SIZES.values())
        contributions['percentage'] = (contributions['length'] / total_length) * 100
        # add only main_class and percentage columns to the contributions dict
        contributions = contributions[['main_class', 'percentage']]
        contributions_dict[col] = contributions
    return contributions_dict

def compute_base_pair_contributions_SNPs(path_SNPs):
    total_length = sum(CHROM_SIZES.values())
    sample_names = None
    snp_counts = {}  # genome -> count of SNPs present

    with gzip.open(path_SNPs, 'rt', encoding='utf-8') as f:
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                header = line.rstrip('\n').split('\t')
                sample_names = header[9:]
                snp_counts = {s: 0 for s in sample_names}
                continue

            fields = line.rstrip('\n').split('\t')
            genotype_fields = fields[9:]
            for sample, genotype in zip(sample_names, genotype_fields):
                gt = genotype.split(':')[0]
                if gt in ['0/1', '1/1', '0|1', '1|1', '1|0', '1/0']:
                    snp_counts[sample] += 1

    # dict of contributions per genome: {genome: DataFrame[main_class, percentage]}
    contributions_dict = {}
    for genome, count in snp_counts.items():
        contributions_dict[genome] = pd.DataFrame({
            'main_class': ['SNP'],
            'percentage': [(count / total_length) * 100]
        })
    return contributions_dict

def combine_contributions(*contribution_dicts):
    """Stack any number of {genome: DataFrame[main_class, percentage]} dicts
    (e.g. one for TEs, one for SNPs) into a single long-format dataframe:
    columns = genome, main_class, percentage. One row per genome per class."""
    all_dfs = []
    for contributions_dict in contribution_dicts:
        for genome, class_df in contributions_dict.items():
            d = class_df.copy()
            d['genome'] = genome
            all_dfs.append(d)
    combined = pd.concat(all_dfs, ignore_index=True)
    return combined[['genome', 'main_class', 'percentage']]

def plot_contributions_boxplot(combined_df, out_path='contributions_boxplot.png'):
    """One box per main_class (TE classes + SNP), each built from the
    per-genome percentage values for that class."""
    classes = [c for c in combined_df['main_class'].unique() if c != 'SNP']
    if 'SNP' in combined_df['main_class'].unique():
        classes.append('SNP')  # keep SNP as the last box

    data = [combined_df.loc[combined_df['main_class'] == c, 'percentage'].values for c in classes]
    colors = [TE_COLORS.get(c, '#AAAAAA') for c in classes]

    fig, ax = plt.subplots(figsize=(10, 6))
    bp = ax.boxplot(data, tick_labels=classes, patch_artist=True)
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.8)

    ax.set_ylabel('% of genome')
    ax.set_xlabel('Main class')
    ax.set_title('Per-genome genome contribution by TE class and SNPs')
    plt.xticks(rotation=45, ha='right')
    fig.tight_layout()
    plt.savefig(out_path)
    plt.close()


df = pd.read_csv(
    "/work/users/g/a/gabivla/lab/SV_mosquitos/TE_paper_aedesaegytpi/reference_vcf/genotyping_reference_TEs/post_processing/merged_results_5_2_either_filtered.csv"
)
# plotting the distributions of lengths for each main class of TEs
df = define_main_class(df)
df = calculate_lengths(df)
plot_length_distributions(df)

# contribution of each main class to the total length of base pairs, per genome
contributions_TEs_dict = compute_base_pair_contributions_TE(df)

# contribution of SNPs to the total length of base pairs, per genome
path_SNPs = '/users/g/a/gabivla/VCF_SNP_variant/whole_sample_filtered_snps_110722.vcf.gz'
contributions_SNPs_dict = compute_base_pair_contributions_SNPs(path_SNPs)

# combine both into one long dataframe: genome, main_class, percentage
combined_contributions = combine_contributions(contributions_TEs_dict, contributions_SNPs_dict)

# boxplot: one box per main_class (TE classes + SNP), each point = one genome
plot_contributions_boxplot(combined_contributions)