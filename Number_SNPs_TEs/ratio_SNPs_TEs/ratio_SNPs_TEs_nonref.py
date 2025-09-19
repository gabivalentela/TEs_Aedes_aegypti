import pandas as pd
import polars as pl
import gzip

reference_TEs = '/work/users/g/a/gabivla/lab/SV_mosquitos/VCF_genotypes_TEs/separate_chromosomes_windows_reference/merged_nonreference/VCF_nonreference_post_processing.csv'
snps_file = '/proj/dschridelab/tvkent/aedes_data/VCFs/AaegL5_full_filtered_snps_110122_rep_map_100kbgenes_centromere_filtered_sites_mask.vcf.gz'

# counting number of TEs
ref_df = pd.read_csv(reference_TEs)

# Read only necessary columns from location file
location_path = '/nas/longleaf/home/gabivla/lab/SV_mosquitos/location_information/mcclintock_output_locations.csv'
location_df = pd.read_csv(location_path, usecols=['Genome_code', 'country'])
location_dict = dict(zip(location_df['Genome_code'], location_df['country']))

# Define base columns that should not be modified for TEs
te_base_cols = ['chromosome', 'position_start', 'position_end', 'TE_name', 'strand']

# Get genome columns at once
genome_cols = [col for col in ref_df.columns if col not in te_base_cols]

# Create rename dictionary and apply it once
rename_dict = {col: f"{col}_{location_dict.get(col, '')}" for col in genome_cols}
ref_df = ref_df.loc[:, ~ref_df.columns.str.contains("Unnamed")]
ref_df = ref_df.rename(columns=rename_dict)

# Extract unique countries from the renamed genome columns only
countries = set()
for original_col in genome_cols:
    new_col = rename_dict[original_col]
    country = new_col.split('_')[-1]
    countries.add(country)

ref_df = pl.DataFrame(ref_df)
ref_df = ref_df.to_pandas()

# Initialize the dictionary to store all information
all_info_by_country = {}

# Process Reference TEs
for country in countries:
    # Get all columns for this country from the renamed genome columns
    country_cols = [col for col in ref_df.columns
                    if col not in te_base_cols and col.endswith('_' + country)]
    number_of_genomes = len(country_cols)
    
    if not country_cols:
        continue

    country_df = ref_df[te_base_cols + country_cols]
    country_df = country_df.dropna(subset=country_cols, how='all')
    # remove rows were there are only 0/0 in the genotype columns
    genotypes_df = ref_df[country_cols]
    country_df = country_df[~(genotypes_df == '0/0').all(axis=1)]
    # get count of reference TEs
    number_ref_TEs = country_df[country_df['TE_name'].notnull()].shape[0]
    
    # Initialize country entry if it doesn't exist
    if country not in all_info_by_country:
        all_info_by_country[country] = {}
    
    all_info_by_country[country]['number_ref_TEs'] = number_ref_TEs

# counting number of SNPs
def read_vcf_with_header(file_path):
    with gzip.open(file_path, 'rt') as f:  # Use gzip.open with 'rt' mode
        lines = f.readlines()

    # Find the header line (starts with single #, not ##)
    header_idx = None
    for i, line in enumerate(lines):
        if line.startswith('#') and not line.startswith('##'):
            header_idx = i
            break

    if header_idx is None:
        raise ValueError("No header line found")

    # Read from header line onwards
    return pd.read_csv(file_path, compression='gzip', sep='\t',
                       skiprows=header_idx)

snps_df = read_vcf_with_header(snps_file)
# Define base columns for SNPs
snp_base_cols = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']

# Get genome columns at once
genome_cols = [col for col in snps_df.columns if col not in snp_base_cols]

# Create rename dictionary and apply it once
rename_dict = {col: f"{col}_{location_dict.get(col, '')}" for col in genome_cols}
snps_df = snps_df.loc[:, ~snps_df.columns.str.contains("Unnamed")]
snps_df = snps_df.rename(columns=rename_dict)

# Extract unique countries from the renamed genome columns only
countries_snps = set()
for original_col in genome_cols:
    new_col = rename_dict[original_col]
    country = new_col.split('_')[-1]
    countries_snps.add(country)

# Process SNPs
for country in countries_snps:
    # Get all columns for this country from the renamed genome columns
    country_cols = [col for col in snps_df.columns
                    if col not in snp_base_cols and col.endswith('_' + country)]
    number_of_genomes = len(country_cols)
    
    if not country_cols:
        continue

    country_df = snps_df[snp_base_cols + country_cols]
    country_df = country_df.dropna(subset=country_cols, how='all')
    # remove rows were there are only 0/0 in the genotype columns
    genotypes_df = snps_df[country_cols]
    genotypes_df = genotypes_df.applymap(lambda x: x.split(':')[0] if isinstance(x, str) else x)
    country_df = country_df[~(genotypes_df == '0/0').all(axis=1)]
    number_snps = country_df[country_df['POS'].notnull()].shape[0]
    
    # Initialize country entry if it doesn't exist
    if country not in all_info_by_country:
        all_info_by_country[country] = {}
    
    # Update the existing dictionary instead of overwriting
    all_info_by_country[country]['number_snps'] = number_snps

# Print the combined results
for country, info in all_info_by_country.items():
    ref_tes = info.get('number_ref_TEs', 'N/A')
    snps = info.get('number_snps', 'N/A')
    print(f"{country}: Reference TEs = {ref_tes}, SNPs = {snps}")

print("\nFull dictionary:")
print(all_info_by_country)

for country, data in all_info_by_country.items():
    ref_data = data['number_ref_TEs']
    snp_data = data['number_snps']
    ratio = ref_data / snp_data if snp_data != 0 else 'undefined'
    print(f"{country}: TE/SNP ratio = {ratio}")
    print('nonreference results')