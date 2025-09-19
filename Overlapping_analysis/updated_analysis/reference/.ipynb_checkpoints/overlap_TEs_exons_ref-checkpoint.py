from utils import read_data_annotation
from utils import divide_by_countries
import polars as pl
import os
import pandas as pd
import sys

def find_overlaps(ref_country, exons, country_name, output_dir='/work/users/g/a/gabivla/lab/SV_mosquitos/overlap_analysis'):
    print(f'Checking for overlap with exons in {country_name}')
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Filter country data for the specific country first
    df_country = ref_country[country_name]
    
    # Collect exon data as lists to reduce repeated iterations
    exon_seqids = exons['seqid'].to_list()
    exon_starts = exons['start'].to_list()
    exon_ends = exons['end'].to_list()
    
    # Prepare lists to collect overlap results
    results = []
    
    # More efficient iteration approach
    for row in df_country.iter_rows(named=True):
        for idx, exon_seqid in enumerate(exon_seqids):
            # Quick chromosome check first
            if row['chromosome'] != exon_seqid:
                continue
            
            # Overlap calculations
            exon_start = exon_starts[idx]
            exon_end = exon_ends[idx]
            
            # Determine overlap type
            if row['position_start'] >= exon_start and row['position_start'] <= exon_end:
                overlap_type = "start_overlap"
            elif row['position_end'] >= exon_start and row['position_end'] <= exon_end:
                overlap_type = "end_overlap"
            elif row['position_start'] <= exon_start and row['position_end'] >= exon_end:
                overlap_type = "full_overlap"
            else:
                continue  # No overlap
            
            # Collect results
            results.append({
                "TE_name": row['TE_name'],
                "chromosome": row['chromosome'],
                "position_start": row['position_start'],
                "position_end": row['position_end'],
                "overlap_type": overlap_type,
                "exon_start": exon_start,
                "exon_end": exon_end,
                "exon_id": exons['id'][idx] if 'id' in exons.columns else None,
                "exon_type": exons['type'][idx] if 'type' in exons.columns else None,
                "exon_seqid": exon_seqid,
                "attributes": exons['attributes'][idx] if 'attributes' in exons.columns else None
            })
    
    # Convert results to Polars DataFrame
    overlaps = pl.DataFrame(results) if results else pl.DataFrame()
    
    # Save to CSV file
    output_path = os.path.join(output_dir, f'{country_name}_overlaps_ref.csv')
    overlaps.write_csv(output_path)
    print(f'Overlaps saved to {output_path}')
    
    return overlaps

path_annotation = '../../data_annotation/VectorBase-54_AaegyptiLVP_AGWG.gff'
df_annotation = read_data_annotation(path_annotation)
# retrieve the rows that have 'exon' in the 'type' column
exons = df_annotation[df_annotation['type'] == 'exon']
    
# Convert exons DataFrame to polars
exons = pl.from_pandas(exons)

info_loc = '/nas/longleaf/home/gabivla/lab/SV_mosquitos/location_information/mcclintock_output_locations.csv'
# Read only necessary columns from location file
location_df = pd.read_csv(info_loc, usecols=['Genome_code', 'country'])
location_dict = dict(zip(location_df['Genome_code'], location_df['country']))

# VCF reference
reference_vcf = '/work/users/g/a/gabivla/lab/SV_mosquitos/VCF_genotypes_TEs/separate_chromosomes_windows_reference/merged_reference/VCF_reference_TEs_post_processing_600.csv'
df_reference = pd.read_csv(reference_vcf)
# remove first column
df_reference = df_reference.iloc[:, 1:]
df_reference = divide_by_countries(df_reference, location_dict)

country = sys.argv[1]

# Find overlaps
find_overlaps(df_reference, exons, country)
