import pandas as pd
import os
import sys
import pysam
import multiprocessing as mp
from collections import defaultdict
from functools import partial
import re
import numpy as np
import csv

_TE_NAME_BATCH_RE = re.compile(r'(.+?)_(AaegL\d+_\d)_(\d+)_(\d+)_([-+])')
_TE_NAME_UPDATE_RE = re.compile(r'(.+)_([A-Z].+)_(\d+)_(\d+)_([-+])_(.+)')


def find_insert_size_cutoff(path):
    if not os.path.exists(path):
        print(f"!! INSERT SIZE PATH NOT FOUND: {path}")
        return {}
    try:
        list_tables = os.listdir(path)
    except FileNotFoundError:
        print(f"!! CANNOT LIST INSERT SIZE PATH: {path}")
        return {}
    dict_cutoff_values = {}
    for csv in list_tables:
        if csv.endswith('.csv'):
            try:
                match = re.search('insert_sizes_(.+).csv', csv)
                if match:
                    genome_code = match.group(1)
                    df = pd.read_csv(os.path.join(path, csv))
                    df = df[df['Cumulative Percentage'] >= 95]
                    if not df.empty:
                        insert_size_cutoff = df.iloc[0]['Insert Size']
                        if genome_code.startswith('SRR'):
                            genome_code = f'{genome_code}_1'
                        else:
                            genome_code = f'{genome_code}_L001_L002_R1_001'
                        dict_cutoff_values[genome_code] = insert_size_cutoff
                    else:
                        print(f"!! NO CUTOFF DATA found in {csv} for genome '{genome_code}'")
            except Exception as e:
                print(f"!! ERROR reading insert size file {csv}: {e}")
                continue

    print(f"\n[find_insert_size_cutoff] Found insert size cutoffs for {len(dict_cutoff_values)} genomes")
    return dict_cutoff_values


def find_bam_files(general_path, genome):
    list_bam_files = {}
    list_files = os.listdir(general_path)
    subset_genome = [file for file in list_files if genome in file and '.final.bam' in file and file.endswith('.bam')]
    list_bam_files[genome] = f'{general_path}/{subset_genome[0]}'
    return list_bam_files

def find_absence_files(list_temp2_paths):
    list_absence_files = {}
    genomes_without_absence_files = []
    for path in list_temp2_paths:
        if not os.path.exists(path):
            print(f"!! ABSENCE FILES PATH NOT FOUND: {path}")
            continue
        try:
            list_files = os.listdir(path)
        except FileNotFoundError:
            print(f"!! CANNOT LIST ABSENCE FILES PATH: {path}")
            continue
        for file in list_files:
            new_path = os.path.join(path, file)
            try:
                files = os.listdir(new_path)
            except FileNotFoundError:
                print(f"!! SUBDIRECTORY NOT FOUND: {new_path}")
                continue
            for genome in files:
                if genome.endswith('1'):
                    temp2_path = os.path.join(new_path, genome, 'results', 'temp2')
                    if not os.path.exists(temp2_path):
                        print(f"!! TEMP2 DIRECTORY NOT FOUND for genome '{genome}': {temp2_path}")
                        genomes_without_absence_files.append(genome)
                        continue
                    try:
                        list_temp2_files = os.listdir(temp2_path)
                    except FileNotFoundError:
                        print(f"!! CANNOT LIST TEMP2 DIRECTORY for genome '{genome}': {temp2_path}")
                        genomes_without_absence_files.append(genome)
                        continue
                    absence_found = False
                    for temp2_file in list_temp2_files:
                        if temp2_file.endswith('.absent.bed'):
                            full_path = os.path.join(temp2_path, temp2_file)
                            # Keep genome name as-is — no renaming of SRR*_1
                            df = pd.read_csv(full_path, sep='\t', header=None, names=['chromosome', 'start', 'end'])
                            list_absence_files[genome] = df
                            absence_found = True
                    if not absence_found:
                        print(f"!! NO .absent.bed FILE FOUND for genome '{genome}' in {temp2_path}")
                        genomes_without_absence_files.append(genome)

    print(f"\n[find_absence_files] Found absence files for {len(list_absence_files)} genomes")
    if genomes_without_absence_files:
        print(f"[find_absence_files] Genomes WITHOUT absence files ({len(genomes_without_absence_files)}): {set(genomes_without_absence_files)}")

    return list_absence_files


def get_genome_codes_TE(df):
    base_columns = ['chromosome', 'strand', 'TE_name', 'position_start', 'position_end', 'score']
    genotype_columns = [col for col in df.columns if col not in base_columns]
    column_index = {column: idx for idx, column in enumerate(df.columns)}

    te_genotype_dict = {}
    for row in df.itertuples(index=False, name=None):
        te_name = (
            f"{row[column_index['TE_name']]}_{row[column_index['chromosome']]}_"
            f"{row[column_index['position_start']]}_{row[column_index['position_end']]}_{row[column_index['strand']]}"
        )
        sample_genotypes = {sample: row[column_index[sample]] for sample in genotype_columns}
        sample_genotypes = {sample: gt for sample, gt in sample_genotypes.items() if gt == '1/1' or gt == '0/0'}
        te_genotype_dict[te_name] = sample_genotypes

    return te_genotype_dict

def is_fr_orientation(read):
    """
    True if this read and its mate point inward toward each other: -> <-
    (i.e. proper FR / 'innie' orientation), False otherwise.
    """
    # must be on opposite strands
    if read.is_reverse == read.mate_is_reverse:
        return False

    if read.reference_start < read.next_reference_start:
        # this read is the upstream (left) mate -> it must be forward
        return not read.is_reverse
    elif read.reference_start > read.next_reference_start:
        # this read is the downstream (right) mate -> it must be reverse
        return read.is_reverse
    else:
        # same start position - can't order by position, fall back to
        # just requiring opposite strands (already checked above)
        return True

def is_qualifying_read_spanning(read):
    if read.mapping_quality == 0:
        return False
    if read.is_unmapped or read.mate_is_unmapped:
        return False
    if read.is_duplicate or read.is_secondary or read.is_supplementary:
        return False
    if read.next_reference_id != read.reference_id:
        return False
    if not is_fr_orientation(read):
        return False
    return True

def is_qualifying_read_split(read):
    if read.mapping_quality == 0:
        return False
    if read.is_unmapped:
        return False
    if read.is_duplicate or read.is_secondary:
        return False
    if not read.mate_is_unmapped and read.next_reference_id == read.reference_id:
        if not is_fr_orientation(read):
            return False
    return True

REF_CONSUMING_OPS = set("MDN=X")

def cigar_ref_length(cigar_str):
    length = 0
    for count, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar_str):
        if op in REF_CONSUMING_OPS:
            length += int(count)
    return length

def get_mate_end_position(read, bam=None):
    if read.is_unmapped or read.mate_is_unmapped:
        return None, None
    if not read.has_tag("MC"):
        print(f"{read.query_name}: no MC tag found")
        return None, None
    mate_cigar = read.get_tag("MC")
    mate_end = read.next_reference_start + cigar_ref_length(mate_cigar)

    mate_mapq = read.get_tag("MQ") if read.has_tag("MQ") else None
    if mate_mapq is None:
        print(f"{read.query_name}: no MQ tag found")
    return mate_end, mate_mapq

def spanning_support(read, chromosome, start, end, bam, insert_threshold, genome_code):
    support_read = set()
    read_length = read.query_length
    tolerance = insert_threshold - (2*read_length)

    read_start = read.reference_start
    read_end = read.reference_end
    mate_start = read.next_reference_start
    mate_end, mate_mapq = get_mate_end_position(read, bam)

    if read_end > start - tolerance and read_end < start and chromosome == read.reference_name:
        if mate_start > end and mate_start < end + tolerance and chromosome == read.next_reference_name:
            #print(f"[{genome_code}] spanning-pair support: {read.query_name} | "
            #      f"upstream={read.reference_start}-{read.reference_end}")
            support_read.add(read)

    if read_start > end and read_start < end + tolerance and chromosome == read.reference_name:
        if mate_end > start - tolerance and mate_end < start and chromosome == read.next_reference_name:
            #print(f"[{genome_code}] spanning-pair support: {read.query_name} | "
            #     f"downstream={read.reference_start}-{read.reference_end}")
            support_read.add(read)

    support_read = list(support_read)

    return support_read


def clips_at_boundary(read, boundary, boundary_tolerance=10):
    left_bound = boundary + boundary_tolerance
    right_bound = boundary - boundary_tolerance
    cigar = read.cigartuples
    if not cigar:
        return False
    left_clip = cigar[0][1] if cigar[0][0] in (4, 5) else 0
    right_clip = cigar[-1][1] if cigar[-1][0] in (4, 5) else 0
    if left_clip > 10 and read.reference_end is not None and read.reference_end < left_bound:
        return True
    if right_clip > 10 and read.reference_start is not None and read.reference_start > right_bound:
        return True
    return False

def split_support(read, start, end, chromosome, genome_code, bam, insert_threshold, boundary_tolerance=10):
    support_read_start = set()
    support_read_end = set()
    read_start = read.reference_start
    read_end = read.reference_end
    mate_start = read.next_reference_start
    mate_end, mate_mapq = get_mate_end_position(read, bam)
    read_length = read.query_length
    tolerance = insert_threshold - read_length

    if read.cigartuples is None:
        return support_read_start, support_read_end

    if read.reference_end < start and read.reference_end > start - boundary_tolerance and chromosome == read.reference_name:
        boundary = start
        boundary_side = 'start'
        mate_on_opposite_side = mate_start > end and mate_start < end + tolerance and chromosome == read.next_reference_name
    elif read.reference_start > end and read.reference_start < end + boundary_tolerance and chromosome == read.reference_name:
        boundary = end
        boundary_side = 'end'
        mate_on_opposite_side = (mate_end is not None) and (mate_end < start and mate_end > start - tolerance) and chromosome == read.next_reference_name
    else:
        return support_read_start, support_read_end

    if not mate_on_opposite_side:
        return support_read_start, support_read_end

    if clips_at_boundary(read, boundary, boundary_tolerance):
        cigar = read.cigartuples
        left_clip = cigar[0][1] if cigar[0][0] in (4, 5) else 0
        right_clip = cigar[-1][1] if cigar[-1][0] in (4, 5) else 0
        #print(f"[{genome_code}] split-support read at boundary {boundary}: "
        #      f"{read.query_name} | flag={read.flag} | "
        #      f"pos={read.reference_start}-{read.reference_end} | "
        #      f"left_clip={left_clip} right_clip={right_clip}")
        if boundary_side == 'start':
            support_read_start.add(read)
        else:
            support_read_end.add(read)

    support_read_start = list(support_read_start)
    support_read_end = list(support_read_end)

    return support_read_start, support_read_end


def evidence_for_TE_absence(reads, chromosome, start, end, insert_threshold, bam, genome_code=None, min_reads_per_side=1, sides_support_TEs='both'):
    if reads is None:
        return None

    reads_support_spanning_list = []
    reads_support_split_list_start = []
    reads_support_split_list_end = []

    for read in reads:
        if is_qualifying_read_spanning(read):
            reads_support_spanning = spanning_support(read, chromosome, start, end, bam, insert_threshold, genome_code)
            reads_support_spanning_list.extend(reads_support_spanning)
        if is_qualifying_read_split(read):
            reads_support_split_start, reads_support_split_end = split_support(read, start, end, chromosome, genome_code, bam,insert_threshold, boundary_tolerance=10)
            reads_support_split_list_start.extend(reads_support_split_start)
            reads_support_split_list_end.extend(reads_support_split_end)

    # sum of supporting reads from spanning and split reads
    #total_supporting_reads = len(reads_support_spanning_list) + len(reads_support_split_list_start) + len(reads_support_split_list_end)
    # make sure no duplicate reads are counted
    # get set of unique read names from all supporting reads
    unique_supporting_reads = set()
    unique_supporting_reads.update([r.query_name for r in reads_support_spanning_list])
    unique_supporting_reads.update([r.query_name for r in reads_support_split_list_start])
    unique_supporting_reads.update([r.query_name for r in reads_support_split_list_end])
    total_supporting_reads = len(unique_supporting_reads)

    if sides_support_TEs == 'both':
        if len(reads_support_spanning_list) > 0 and total_supporting_reads >= min_reads_per_side:
            return '0/0'
        elif (len(reads_support_split_list_start) > 0 or len(reads_support_split_list_end) > 0) and total_supporting_reads >= min_reads_per_side:
            return '0/0'
    elif sides_support_TEs == 'either':
        if len(reads_support_spanning_list) > 0 and total_supporting_reads >= min_reads_per_side:
            return '0/0'
        elif (len(reads_support_split_list_start) > 0 or len(reads_support_split_list_end) > 0) and total_supporting_reads >= min_reads_per_side:
            return '0/0'
    return None

def evidence_te_presence(reads, chromosome, start, end, insert_threshold, bam_file,
                          genome_code=None, min_overlap=20,
                          min_mapq=1, min_reads_per_side=2, sides_support_TEs='both'):
    if reads is None:
        print(f"[{genome_code}] !! READS NOT FOUND for genome '{genome_code}' - cannot check TE presence")
        return './.', {}

    te_key = f'{chromosome}:{start}_{end}'
    type_of_evidence = {te_key: {}}  # te_key -> {read_name: {..., 'evidence': [ {...}, {...} ]}}

    reads_presence_start = set()
    reads_presence_end = set()

    for read in reads:
        if (read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate
                or read.mapping_quality < min_mapq):
            continue

        read_start = read.reference_start
        read_end = read.reference_end
        if read_end is None:
            continue

        mate_start = read.next_reference_start
        mate_end, mate_mapq = get_mate_end_position(read, bam_file)
        has_valid_mate = (read.is_paired and not read.mate_is_unmapped and read.is_proper_pair
                   and mate_start is not None and mate_start >= 0 and mate_end is not None
                   and chromosome == read.next_reference_name and mate_mapq is not None and mate_mapq >= min_mapq)

        # read-side overlap calcs (always available)
        right_most_pos_read = max(read_start, read_end)
        leftmost_pos_read = min(read_start, read_end)
        overlap_inside_TE_read = max(min(right_most_pos_read, end) - max(leftmost_pos_read, start), 0)
        outside_of_TE_read_start = max(0, min(right_most_pos_read, start) - leftmost_pos_read)
        outside_of_TE_read_end   = max(0, right_most_pos_read - max(leftmost_pos_read, end))

        # mate-side overlap calcs (only meaningful/trusted if has_valid_mate)
        if has_valid_mate:
            right_most_pos_mate = max(mate_start, mate_end)
            leftmost_pos_mate = min(mate_start, mate_end)
            overlap_inside_TE_mate = max(min(right_most_pos_mate, end) - max(leftmost_pos_mate, start), 0)
            outside_of_TE_mate_start = max(0, min(right_most_pos_mate, start) - leftmost_pos_mate)
            outside_of_TE_mate_end   = max(0, right_most_pos_mate - max(leftmost_pos_mate, end))
        else:
            overlap_inside_TE_mate = outside_of_TE_mate_start = outside_of_TE_mate_end = 0

        def record_evidence(side, evidence_type, used_mate):
            """Append evidence for this read/TE, keeping full history (mapq etc.) for later filtering."""
            entry = type_of_evidence[te_key].setdefault(
                read.query_name,
                {
                    'mapping_quality': read.mapping_quality,
                    'mate_mapping_quality': mate_mapq,
                    'has_valid_mate': has_valid_mate,
                    'evidence': [],
                }
            )
            entry['evidence'].append({
                'side': side,
                'evidence_type': evidence_type,
                'used_mate': used_mate,
            })

        # --- just_read: read ALONE can support this; mate (if valid) is an alternate/independent path ---
        if overlap_inside_TE_read >= min_overlap and outside_of_TE_read_start >= min_overlap:
            reads_presence_start.add(read.query_name)
            record_evidence('start', 'just_read_start', used_mate=False)
        if has_valid_mate and overlap_inside_TE_mate >= min_overlap and outside_of_TE_mate_start >= min_overlap:
            reads_presence_start.add(read.query_name)
            record_evidence('start', 'just_read_start_via_mate', used_mate=True)

        if overlap_inside_TE_read >= min_overlap and outside_of_TE_read_end >= min_overlap:
            reads_presence_end.add(read.query_name)
            record_evidence('end', 'just_read_end', used_mate=False)
        if has_valid_mate and overlap_inside_TE_mate >= min_overlap and outside_of_TE_mate_end >= min_overlap:
            reads_presence_end.add(read.query_name)
            record_evidence('end', 'just_read_end_via_mate', used_mate=True)

        # --- read_and_mate: strictly requires a valid mate, since it combines read + mate overlap ---
        if has_valid_mate:
            if (overlap_inside_TE_read >= min_overlap and outside_of_TE_mate_start >= min_overlap) or \
               (overlap_inside_TE_mate >= min_overlap and outside_of_TE_read_start >= min_overlap):
                reads_presence_start.add(read.query_name)
                record_evidence('start', 'read_and_mate_start', used_mate=True)

            if (overlap_inside_TE_read >= min_overlap and outside_of_TE_mate_end >= min_overlap) or \
               (overlap_inside_TE_mate >= min_overlap and outside_of_TE_read_end >= min_overlap):
                reads_presence_end.add(read.query_name)
                record_evidence('end', 'read_and_mate_end', used_mate=True)

    total_supporting_reads = len(reads_presence_start | reads_presence_end)

    genotype = './.'
    if sides_support_TEs == 'both':
        if len(reads_presence_start) > 0 and len(reads_presence_end) > 0 and total_supporting_reads >= min_reads_per_side:
            genotype = '1/1'
    elif sides_support_TEs == 'either':
        if (len(reads_presence_start) > 0 or len(reads_presence_end) > 0) and total_supporting_reads >= min_reads_per_side:
            genotype = '1/1'

    return genotype, type_of_evidence

def evidence_for_TE_absence_file(df_absence, chromosome, start, end, genome):
    if df_absence is None:
        return './.'
    overlap = df_absence[
        (df_absence['chromosome'] == chromosome) &
        (df_absence['start'] == start) &
        (df_absence['end'] == end)
    ]
    if not overlap.empty:
        #print(f"[{genome}] Evidence for TE absence from absence file: {chromosome}:{start}-{end} overlaps with "
        #      f"{overlap.iloc[0]['chromosome']}:{overlap.iloc[0]['start']}-{overlap.iloc[0]['end']}")
        return '0/0'
    return './.'



def process_te_batch(te_batch, bam_files, cutoff_values, full_absence_list, number_supporting_reads_absence, number_supporting_reads_presence, sides_support_TEs):
    results = {}
    
    if not te_batch:
        return results

    genome = te_batch[0][1]
    bam_file = bam_files.get(genome, 'Not Found')
    bam = None if bam_file == 'Not Found' else pysam.AlignmentFile(bam_file, "rb")

    absence_file = full_absence_list.get(genome, None)
    df_absence = None
    if absence_file is not None:
        df_absence = absence_file

    try:
        for te_info in te_batch:
            TE_name, genome, chromosome, start, end = te_info
            reads = list(bam.fetch(chromosome, max(0, start - 1000), end + 1000))

            insert_threshold = int(cutoff_values.get(genome, 1000))  # Default to 1000 if not found

            if bam_file == 'Not Found':
                print(f"[{genome}] !! BAM FILE NOT FOUND for genome '{genome}'")
            if genome not in cutoff_values:
                print(f"[{genome}] !! INSERT SIZE CUTOFF NOT FOUND for genome '{genome}' - using default 1000bp")

            te_name_genome = f'{TE_name}_{genome}'
            te_key = f'{chromosome}:{start}_{end}'

            genotype_presence, presence_evidence = evidence_te_presence(
                reads, chromosome, start, end, insert_threshold, bam, genome_code=genome, min_reads_per_side=number_supporting_reads_presence,
                sides_support_TEs=sides_support_TEs
            )
            genotype_absence_reads = evidence_for_TE_absence(
                reads, chromosome, start, end, insert_threshold, bam, genome_code=genome, min_reads_per_side=number_supporting_reads_absence, sides_support_TEs=sides_support_TEs
            )
            genotype_absence_file = evidence_for_TE_absence_file(
                df_absence, chromosome, start, end, genome
            )

            # Normalize None to './.'
            if genotype_presence is None:
                genotype_presence = './.'
            if genotype_absence_reads is None:
                genotype_absence_reads = './.'
            if genotype_absence_file is None:
                genotype_absence_file = './.'

            # Absence file is ground truth — checked first
            if genotype_absence_file == '0/0':
                final_genotype = '0/1' if genotype_presence == '1/1' else '0/0'

            # No file evidence — fall back to read-based evidence
            else:  # genotype_absence_file == './.'
                if genotype_presence == '1/1' and genotype_absence_reads == '0/0':
                    final_genotype = '0/1'
                elif genotype_presence == '1/1' and genotype_absence_reads == './.':
                    final_genotype = '1/1'
                elif genotype_presence == './.' and genotype_absence_reads == '0/0':
                    final_genotype = '0/0'
                else:  # presence './.' and absence_reads './.'
                    final_genotype = './.'

            results[te_name_genome] = final_genotype

    finally:
        if bam is not None:
            bam.close()

    return results


def check_bam_files(dict_absent_tes, bam_files, cutoff_values, strictness, full_absence_list,
                     number_supporting_reads_absence, number_supporting_reads_presence, sides_support_TEs, num_processes=None):
    if num_processes is None:
        num_processes = int(os.environ.get('SLURM_CPUS_PER_TASK', mp.cpu_count()))

    print(f"\n[check_bam_files] Using {num_processes} processes")

    genome_te_batches = defaultdict(list)

    for TE_name, absent_genomes in dict_absent_tes.items():
        match = _TE_NAME_BATCH_RE.search(TE_name)
        if not match:
            print(f'No match found for {TE_name}')
            continue

        chromosome, start, end = match.group(2), int(match.group(3)), int(match.group(4))
        for genome in absent_genomes:
            genome_te_batches[genome].append((TE_name, genome, chromosome, start, end))

    all_batches = []
    for genome, te_list in genome_te_batches.items():
        chrom_batches = defaultdict(list)
        for te_info in te_list:
            chrom_batches[te_info[2]].append(te_info)
        for chrom_batch in chrom_batches.values():
            all_batches.append(chrom_batch)

    results = {}
    all_mapq_records = []
    with mp.Pool(processes=num_processes) as pool:
        if strictness == 'strict':
            batch_results = pool.map(
                partial(process_te_batch, bam_files=bam_files, cutoff_values=cutoff_values, full_absence_list=full_absence_list, number_supporting_reads_absence=number_supporting_reads_absence, number_supporting_reads_presence=number_supporting_reads_presence, sides_support_TEs=sides_support_TEs),
                all_batches
            )
        else:
            pass

        for batch_result in batch_results:
            results.update(batch_result)

    return results

def update_vcf(vcf_file, updated_genotypes, strictness_label):
    print(f"\n[update_vcf] Received {len(updated_genotypes)} keys to update {strictness_label}")
    for te_info, genotype in updated_genotypes.items():
        #print(f"\n[update_vcf] Processing key: '{te_info}' → genotype='{genotype}'")

        match = _TE_NAME_UPDATE_RE.search(te_info)
        if not match:
            print(f"[update_vcf] !! REGEX FAILED for: '{te_info}'")
            continue

        te_name, chromosome, start, end, strand, genome_code = match.groups()
        start, end = int(start), int(end)
        #print(f"[update_vcf] Parsed → te_name='{te_name}' | chrom='{chromosome}' | "
        #      f"start={start} | end={end} | strand='{strand}' | genome='{genome_code}'")

        row_match = vcf_file[
            (vcf_file['chromosome'] == chromosome) &
            (vcf_file['position_start'] == start) &
            (vcf_file['position_end'] == end) &
            (vcf_file['TE_name'] == te_name) &
            (vcf_file['strand'] == strand)
        ]

        if row_match.empty:
            #print(f"[update_vcf] !! NO ROW MATCH for te_name='{te_name}' | chrom='{chromosome}' | "
            #      f"start={start} | end={end} | strand='{strand}'")
            nearby = vcf_file[vcf_file['chromosome'] == chromosome][
                ['chromosome', 'TE_name', 'position_start', 'position_end', 'strand']
            ]
            #print(f"[update_vcf] Rows in vcf_file for chrom '{chromosome}':\n{nearby.to_string()}")
            continue

        row_index = row_match.index[0]
        if genome_code not in vcf_file.columns:
            print(f"[update_vcf] !! Column '{genome_code}' not found in vcf_file")
            continue

        #print(f"[update_vcf] ✓ Updating row {row_index}, column '{genome_code}' to '{genotype}'")
        vcf_file.at[row_index, genome_code] = genotype

    return vcf_file



def main_check_bams(df_has_ref, number_supporting_reads_absence, number_supporting_reads_presence, sides_support_TEs, mapq_output_file=None):
    df = df_has_ref.copy()
    # getting evidence for presence check - for the ones that were 1/1 initially
    genome_codes_no_TE = get_genome_codes_TE(df)
    updated_genotypes_strict = check_bam_files(
        genome_codes_no_TE, list_bam_files, cutoff_values, 'strict', full_absence_list,
        number_supporting_reads_absence, number_supporting_reads_presence, sides_support_TEs
    )
    updated_vcf_strict = update_vcf(df.copy(), updated_genotypes_strict, 'strict')
    return updated_vcf_strict


# ── Entry point ──────────────────────────────────────────────────────────────
chromosome = sys.argv[1]
genome = sys.argv[2]
number_supporting_reads_presence = int(sys.argv[3]) # 1 or 2 or 5 
number_supporting_reads_absence = int(sys.argv[4]) # 1 or 2 or 5 
sides_support_TEs = str(sys.argv[5]) # both or either

list_temp2_paths = [
    '/work/users/g/a/gabivla/lab/SV_mosquitos/new_library_pipeline_reference_rerun/outputs',
    '/users/g/a/gabivla/new_library_pipeline_reference/outputs'
]
bam_paths = '/users/g/a/gabivla/BAM_files_reference_genotyping'
full_absence_list_all = find_absence_files(list_temp2_paths)
full_absence_list = {genome: full_absence_list_all.get(genome)}
del full_absence_list_all
list_bam_files = find_bam_files(bam_paths, genome)
cutoff_values = find_insert_size_cutoff('/work/users/g/a/gabivla/lab/SV_mosquitos/insert_sizes_bam/')
print(f'done finding absence files, bam files and cutoff values')

# making directory for outputs
os.makedirs(f'./results_{number_supporting_reads_absence}_{number_supporting_reads_presence}_{sides_support_TEs}', exist_ok=True)

path_file = '/work/users/g/a/gabivla/lab/SV_mosquitos/TEs_Aedes_aegypti/running_pipeline_new_library/new_run_reference_mcclintock/making_final_vcf/merge_chromosomes/combined_reference_tes.vcf'
# keep only base columns and the specific genome column
columns_to_keep = ['chrom', 'start', 'end', 'strand', 'TE_name', 'score', genome]
df_merged_onecode = pd.read_csv(path_file, sep='\t', usecols=columns_to_keep)
df_merged_onecode = df_merged_onecode[df_merged_onecode['chrom'] == chromosome]
new_column_names = {
    'chrom': 'chromosome',
    'strand': 'strand',
    'TE_name': 'TE_name',
    'start': 'position_start',
    'end': 'position_end'
}
df_merged_onecode.rename(columns=new_column_names, inplace=True)
print(f'done making df')
df_merged_onecode_strict = main_check_bams(df_merged_onecode, number_supporting_reads_absence, number_supporting_reads_presence, sides_support_TEs)
print(f'done checking bams and updating vcf')
df_merged_onecode_strict.to_csv(f'./results_{number_supporting_reads_absence}_{number_supporting_reads_presence}_{sides_support_TEs}/df_merged_onecode_strict_subset_{chromosome}_{genome}.csv', index=False)