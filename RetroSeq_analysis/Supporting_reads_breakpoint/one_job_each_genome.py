import os
import sys
import importlib.util as il
import pysam
import pandas as pd
import re
import numpy as np
import polars as pl
import argparse

# Load bedtools module
os.system('module load bedtools')

def get_reads_from_positions(hybrid_dict, bam_file, retroseq_discovery):

    # Open the BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")
    chrmo = hybrid_dict['chrm']

    # Retrieve the length of the chromosome
    chr_length = bam.get_reference_length(chrmo)
    
    # Adjust start and end positions to ensure they are within valid ranges
    start = max(0, int(hybrid_dict['start']) - 2000)  # Ensure start is non-negative
    end = min(int(hybrid_dict['end']) + 2000, chr_length)  # Ensure end does not exceed chromosome length
    family = hybrid_dict['family']
    
    # Fetch and store only the read names within the specified region
    read_codes = [read.query_name for read in bam.fetch(chrmo, start, end)]
    
    # Read the discovery file
    df_discovery = pl.read_csv(retroseq_discovery, separator="\t", has_header=False)  
    # Assuming the first column in df_discovery contains the read codes
    read_code_column_index = 4

    # Filter the rows in df_discovery where the read code is in the read_codes list
    filtered_discovery = df_discovery.filter(pl.col("column_5").is_in(read_codes))
    
    TE1 = re.search('(.+)-(.+)-hybrid', family).group(1)
    TE2 = re.search('(.+)-(.+)-hybrid', family).group(2)
        
    # Count how many rows contain the specific string in any column
    count_TE1 = filtered_discovery.filter(pl.col("column_4").str.contains(TE1, literal=True)).height
    count_TE2 = filtered_discovery.filter(pl.col("column_4").str.contains(TE2, literal=True)).height

    # Store counts in a dictionary
    te_counts = {TE1: count_TE1, TE2: count_TE2}

    # Get the TE with the highest count
    max_te = max(te_counts, key=te_counts.get)
    print(te_counts)
    print(f"The TE with the highest count is '{max_te}' with {te_counts[max_te]} occurrences.")

    bam.close()
    return max_te

def get_correct_family_TE(info_col, chrm, bam_file, retroseq_discovery):
    if 'hybrid' in info_col:
        start = int(info_col.split(",")[1])
        end = int(info_col.split(",")[2])
        family = re.search('(.+-hybrid)', info_col).group(1)
        dict_te = {
            'chrm': chrm,
            'start': start,
            'end': end,
            'family': family
        }
        te_info = get_reads_from_positions(dict_te, bam_file, retroseq_discovery)
        return te_info
    else:
        te_info = (info_col.split(",")[0]).split("-")[0]
        return te_info
                

def load_config(config_path):
    spec = il.spec_from_file_location("config", config_path)
    config = il.module_from_spec(spec)
    sys.modules[spec.name] = config
    spec.loader.exec_module(config)
    return config

def load_mcc_modules(mcc_path):
    scripts_path = os.path.join(mcc_path, 'scripts')
    print(f"Adding {scripts_path} to sys.path")
    if scripts_path not in sys.path:
        sys.path.append(scripts_path)
    print(f"Current sys.path: {sys.path}")
    
    import mccutils
    import output 
    
    return mccutils, output

def retrieve_chromosome_names(retroseq_out):
    df = pd.read_csv(retroseq_out, sep='\t', comment='#', header=None,
                     names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'])
    chromosomes = df['CHROM'].unique()
    return chromosomes

def main(retroseq_out, reference_fasta, out_dir, ref_name, sample_name, vcf_options, config_path, mcc_path, bam_file, retroseq_discovery, breakpoint_value):
    config = load_config(config_path)
    mccutils, output = load_mcc_modules(mcc_path)

    mccutils.log("retroseq", "processing RetroSeq results")

    chromosomes = retrieve_chromosome_names(retroseq_out)
    
    insertions = read_insertions(bam_file, retroseq_discovery, retroseq_out, sample_name, chromosomes, output, 
                                 support_threshold=config.PARAMS["read_support_threshold"], 
                                 breakpoint_threshold=config.PARAMS["breakpoint_confidence_threshold"])
    
    if len(insertions) >= 1:
        insertions = output.make_redundant_bed(insertions, sample_name, out_dir, method="retroseq")
        insertions = output.make_nonredundant_bed(insertions, sample_name, out_dir, method="retroseq")
        output.write_vcf(insertions, reference_fasta, sample_name, "retroseq", out_dir, vcf_options)
    else:
        mccutils.run_command(["touch", os.path.join(out_dir, f"{sample_name}_retroseq_redundant_{breakpoint_value}.bed")])
        mccutils.run_command(["touch", os.path.join(out_dir, f"{sample_name}_retroseq_nonredundant_{breakpoint_value}.bed")])
    
    mccutils.log("retroseq", "RetroSeq post processing complete")

def read_insertions(bam_file, retroseq_discovery, retroseq_vcf, sample_name, chromosomes, output, support_threshold=0, breakpoint_threshold=6):
    insertions = []

    with open(retroseq_vcf, "r") as vcf:
        for line in vcf:
            if "#" not in line:
                insert = output.Insertion(output.Retroseq())
                line = line.replace("\n","")
                split_line = line.split("\t")
                insert.chromosome = split_line[0]

                info = {}
                split_info = split_line[7].split(";")
                for i in split_info:
                    if "=" in i:
                        info[i.split("=")[0]] = i.split("=")[1]

                insert.start = int(info['MEINFO'].split(",")[1])
                insert.end = int(info['MEINFO'].split(",")[2])
                insert.family = get_correct_family_TE(info['MEINFO'], split_line[0], bam_file, retroseq_discovery)

                format_keys = split_line[8].split(":")
                format_vals = split_line[9].split(":")
                form = {}
                for x, key in enumerate(format_keys):
                    form[key] = format_vals[x]
                
                insert.support_info.support['spanning_pairs'].value = int(form['SP'])
                insert.support_info.support['supporting_reads'].value = int(form['GQ'])
                insert.support_info.support['clip3'].value = int(form['CLIP3'])
                insert.support_info.support['clip5'].value = int(form['CLIP5'])
                insert.support_info.support['call_status'].value = int(form['FL'])
                insert.support_info.support['frequency'].value = round(int(form['GQ'])/((2 * int(form['SP'])) + int(form['GQ'])), 2)
                insert.type = "non-reference"
                insert.name = f"{insert.family}|non-reference|{insert.support_info.support['frequency'].value}|{sample_name}|retroseq|rp|"

                if (insert.support_info.support['supporting_reads'].value >= support_threshold and 
                    insert.support_info.support['call_status'].value >= breakpoint_threshold and 
                    insert.chromosome in chromosomes):
                    insertions.append(insert)
    
    return insertions

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process RetroSeq data.")
    parser.add_argument("retroseq_out", help="Path to RetroSeq output file")
    parser.add_argument("reference_fasta", help="Path to reference FASTA file")
    parser.add_argument("out_dir", help="Directory to store output files")
    parser.add_argument("ref_name", help="Reference name")
    parser.add_argument("sample_name", help="Sample name")
    parser.add_argument("vcf_options", nargs="+", help="Options for VCF file generation")
    parser.add_argument("config_path", help="Path to configuration file")
    parser.add_argument("mcc_path", help="Path to MCC modules directory")
    parser.add_argument("bam_file", help="Path to BAM file")
    parser.add_argument("retroseq_discovery", help="Path to RetroSeq discovery file")
    parser.add_argument("breakpoint_value", type=int, help="Value for breakpoint - 6 or 1")  # Ensure it's parsed as an integer
    
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    
    main(args.retroseq_out, args.reference_fasta, args.out_dir, args.ref_name, args.sample_name, args.vcf_options, args.config_path, args.mcc_path, args.bam_file, args.retroseq_discovery, args.breakpoint_value)
