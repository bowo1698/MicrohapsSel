#!/usr/bin/env python3
"""
Convert CSV genotype to VCF format
Modified from Dzianis Prakapenka to accept CSV input and optimised memory with parallel processing

Author: Agus Wibowo
Date: 2025-01-26

Geneotype data expected:
ID,snp_1,snp_2,snp_3,...,generation,population_type,individual_id
ind_1,0,1,2,...,5,reference,ref_gen5_ind1
ind_2,1,2,0,...,5,reference,ref_gen5_ind2
ind_3,2,0,1,...,5,reference,ref_gen5_ind3

Map file:
SNPID    Chr    Position
snp_1    1      12345
snp_2    1      23456
snp_3    2      34567
snp_4    2      45678

Example:
1. Single VCF file (sequential)
python convert-to-vcf.py \
  -i genotype_reference.csv \
  -m snp_map.txt \
  -o salmon_genotypes.vcf

2. Per-Chromosome VCF (Sequential)
python convert-to-vcf.py \
  -i genotype_reference.csv \
  -m snp_map.txt \
  -o salmon_genotypes.vcf \
  --split-by-chr

3. Parallel Processing and store all vcf file per chr
python convert-to-vcf.py \
  -i genotype_reference.csv \
  -m snp_map.txt \
  -o salmon_genotypes.vcf \
  --split-by-chr \
  --parallel \
  --ncores 8

4. Parallel + Merge (Delete Chr Files)
python convert-to-vcf.py \
  -i genotype_reference.csv \
  -m snp_map.txt \
  -o salmon_genotypes.vcf \
  --split-by-chr \
  --parallel \
  --ncores 25 \
  --merge-after-split

5. Parallel + Merge (Keep Both)
python convert-to-vcf.py \
  -i genotype_reference.csv \
  -m snp_map.txt \
  -o salmon_genotypes.vcf \
  --split-by-chr \
  --parallel \
  --ncores 25 \
  --merge-after-split \
  --keep-chr-files
"""
from __future__ import print_function
import argparse
import os
import re
from datetime import datetime
import numpy as np
import pandas as pd
from multiprocessing import Pool, cpu_count

def make_arg_parser():
    app_name = "convert-csv-to-vcf.py"
    description = "Converts AlphaSimR genotype CSV to VCF format"
    parser = argparse.ArgumentParser(
        prog=app_name,
        description=description,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-i", "--input",
                       required=True,
                       help="Path to genotype CSV file (e.g., genotype_reference.csv)")
    parser.add_argument("-m", "--map",
                       required=True,
                       help="Path to map file (format: 'SNPID Chr Position')")
    parser.add_argument("--ref",
                       default="T",
                       help="Reference allele")
    parser.add_argument("--alt",
                       default="C",
                       help="Alternative allele")
    parser.add_argument("--qual",
                       default=".",
                       help="VCF QUAL column")
    parser.add_argument("--filter",
                       default="PASS",
                       help="VCF FILTER column")
    parser.add_argument("--info",
                       default=".",
                       help="VCF INFO column")
    parser.add_argument("--format-field",
                       default="GT",
                       help="VCF FORMAT column")
    parser.add_argument("-o", "--output",
                       default="output.vcf",
                       help="Output VCF file name")
    parser.add_argument("--parallel", action="store_true",
                       help="Process chromosomes in parallel (only with --split-by-chr)")
    parser.add_argument("--ncores", type=int, default=cpu_count(),
                       help="Number of cores for parallel processing")
    parser.add_argument("--split-by-chr",
                   action="store_true",
                   help="Create separate VCF files per chromosome")
    parser.add_argument("--merge-after-split",
                    action="store_true",
                    help="Merge chromosome VCFs into single file after parallel processing (requires --split-by-chr --parallel)")
    parser.add_argument("--keep-chr-files",
                    action="store_true",
                    help="Keep individual chromosome VCF files after merging (only with --merge-after-split)")
    parser.add_argument("-v", "--verbose",
                    action="store_true",
                    default=False,
                    help="Verbose output")
    return parser

def natural_sort_key(s):
    """
    Sort strings with numbers naturally (e.g., chr1, chr2, chr10 not chr1, chr10, chr2)
    
    Args:
        s: String to convert to sort key
    
    Returns:
        list: Sort key with integers parsed
    
    Example:
        >>> sorted(['chr10', 'chr2', 'chr1'], key=natural_sort_key)
        ['chr1', 'chr2', 'chr10']
    """
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split('([0-9]+)', s)]

def read_map_file(map_file, verbose=False):
    """
    Read map file and create SNP dictionary
    
    Expected format:
        SNPID    Chr    Position
        snp_1    1      12345
        snp_2    1      23456
    
    Args:
        map_file: Path to map file
        verbose: Print progress messages
    
    Returns:
        dict: {snp_id: [chr_id, position, snp_id]}
    """
    if verbose:
        print(f"Reading map file: {map_file}")
    
    if not os.path.exists(map_file):
        raise FileNotFoundError(f"Map file not found: {map_file}")
    
    snp_dict = {}
    
    with open(map_file, 'r') as f:
        # Skip header line
        _ = f.readline()
        
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 3:
                snp_id = parts[0]
                chr_id = parts[1]
                position = parts[2]
                snp_dict[snp_id] = [chr_id, position, snp_id]
    
    if verbose:
        print(f"  Loaded {len(snp_dict)} SNPs")
    
    return snp_dict


def read_genotype_metadata(csv_file, verbose=False):
    """
    Read metadata from AlphaSimR genotype CSV (header and individual IDs only)
    
    Expected format:
        ID,snp_1,snp_2,...,generation,population_type,individual_id
        ind_1,0,1,2,...,5,reference,reference_gen5_ind1
    
    Args:
        csv_file: Path to genotype CSV file
        verbose: Print progress messages
    
    Returns:
        tuple: (snp_cols, individual_ids)
            - snp_cols: List of SNP column names
            - individual_ids: List of individual IDs
    """
    if verbose:
        print(f"Reading genotype CSV: {csv_file}")
    
    if not os.path.exists(csv_file):
        raise FileNotFoundError(f"Genotype file not found: {csv_file}")
    
    df_header = pd.read_csv(csv_file, nrows=0)
    metadata_cols = ['ID', 'individual', 'generation', 'population_type', 'individual_id']
    snp_cols = [col for col in df_header.columns if col not in metadata_cols]
    
    # Read only ID column
    df_ids = pd.read_csv(csv_file, usecols=['ID'] if 'ID' in df_header.columns else ['individual'])
    individual_ids = df_ids.iloc[:, 0].astype(str).tolist()
    
    return snp_cols, individual_ids

def read_genotype_chunk(csv_file, snp_cols, chunksize=1000):
    """
    Generator that yields genotype data in chunks (memory-efficient streaming)
    
    Args:
        csv_file: Path to genotype CSV file
        snp_cols: List of SNP column names to extract
        chunksize: Number of individuals to process per chunk
    
    Yields:
        tuple: (ind_ids, geno_matrix)
            - ind_ids: numpy array of individual IDs (shape: n_individuals)
            - geno_matrix: numpy int8 array of genotypes (shape: n_individuals × n_snps)
                          Values: 0=homozygous ref, 1=heterozygous, 2=homozygous alt, -1=missing
    """
    use_cols = ['ID'] + snp_cols if 'ID' in pd.read_csv(csv_file, nrows=0).columns else ['individual'] + snp_cols
    
    for chunk in pd.read_csv(csv_file, usecols=use_cols, chunksize=chunksize, dtype=str):
        ind_ids = chunk.iloc[:, 0].values
        # More explicit: convert to numeric first
        geno_df = chunk.iloc[:, 1:].apply(pd.to_numeric, errors='coerce')
        geno_matrix = geno_df.fillna(-1).astype(np.int8).values
        yield ind_ids, geno_matrix

def encode_genotype_vectorized(geno_array):
    """
    Convert numeric genotype array to VCF GT format strings (vectorized)
    
    Args:
        geno_array: numpy array of integer genotypes
    
    Returns:
        numpy array of GT strings
    
    Encoding:
        0  -> '0/0' (homozygous reference)
        1  -> '0/1' (heterozygous)
        2  -> '1/1' (homozygous alternate)
        -1 -> './.' (missing)
    """
    mapping = {0: '0/0', 1: '0/1', 2: '1/1', -1: './.'}
    return np.vectorize(mapping.get)(geno_array)


def write_vcf_header(out_file, individual_ids, verbose=False):
    """
    Write VCF v4.2 header to file
    
    Args:
        out_file: File object to write to
        individual_ids: List of individual/sample IDs for header
        verbose: Print progress messages
    """
    now = datetime.now()
    
    comments = [
        "##fileformat=VCFv4.2",
        f"##fileDate={now.strftime('%Y%m%d')}",
        "##source=AlphaSimR-to-VCF-converter",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
    ]
    
    # Column header
    header = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    header.extend(individual_ids)
    
    # Write
    for comment in comments:
        print(comment, file=out_file)
    print('\t'.join(header), file=out_file)
    
    if verbose:
        print(f"  Header written with {len(individual_ids)} individuals")

def write_vcf_body_streaming(out_file, csv_file, snp_cols, snp_dict, args, chr_filter=None, chunksize=1000, verbose=False):
    """
    Write VCF variant lines using chunked streaming (memory-efficient)
    
    Args:
        out_file: File object to write to
        csv_file: Path to genotype CSV file
        snp_cols: List of all SNP column names
        snp_dict: Dictionary mapping SNP ID to [chr, pos, snp_id]
        args: argparse.Namespace with VCF format parameters (ref, alt, qual, etc.)
        chr_filter: If provided, only write SNPs from this chromosome
        chunksize: Number of individuals to process per chunk
        verbose: Print progress messages
    
    Returns:
        int: Number of variants written
    """
    variants_written = 0
    
    # Filter SNPs by chromosome
    if chr_filter is not None:
        snp_indices = [i for i, snp in enumerate(snp_cols) 
                      if snp in snp_dict and snp_dict[snp][0] == str(chr_filter)]
        filtered_snps = [snp_cols[i] for i in snp_indices]
    else:
        snp_indices = [i for i, snp in enumerate(snp_cols) if snp in snp_dict]
        filtered_snps = [snp_cols[i] for i in snp_indices]
    
    if not filtered_snps:
        return 0
    
    # Process in chunks to manage memory
    all_genotypes = []
    all_ind_ids = []
    
    for ind_ids, geno_matrix in read_genotype_chunk(csv_file, snp_cols, chunksize):
        all_ind_ids.extend(ind_ids)
        all_genotypes.append(geno_matrix[:, snp_indices])
    
    # Concatenate chunks
    full_genotype_matrix = np.vstack(all_genotypes)
    
    # Write variants (transpose: SNPs × Individuals)
    for idx, snp_col in enumerate(filtered_snps):
        chr_id, position, snp_id = snp_dict[snp_col]
        
        variant_line = [chr_id, position, snp_id, args.ref, args.alt,
                       args.qual, args.filter, args.info, args.format_field]
        
        # Vectorized encoding
        genotypes_for_snp = full_genotype_matrix[:, idx]
        gt_strings = encode_genotype_vectorized(genotypes_for_snp)
        variant_line.extend(gt_strings)
        
        print('\t'.join(variant_line), file=out_file)
        variants_written += 1
    
    if verbose:
        print(f"  Wrote {variants_written} variants")
    
    return variants_written

def process_chromosome(args_tuple):
    """
    Worker function for parallel chromosome processing (multiprocessing.Pool)
    
    Args:
        args_tuple: Tuple of (chr_id, csv_file, snp_cols, snp_dict, args, output_file)
    
    Returns:
        tuple: (chr_id, variants_written, output_file_path)
    """
    chr_id, csv_file, snp_cols, snp_dict, args, output_file = args_tuple
    
    # Read metadata
    _, individual_ids = read_genotype_metadata(csv_file)
    
    with open(output_file, 'w') as out_file:
        write_vcf_header(out_file, individual_ids, verbose=False)
        variants = write_vcf_body_streaming(out_file, csv_file, snp_cols, 
                                           snp_dict, args, chr_filter=chr_id, 
                                           chunksize=1000)
    
    return chr_id, variants, output_file

def merge_vcf_files(chr_files, output_file, keep_chr_files=False, verbose=False):
    """
    Merge multiple chromosome-specific VCF files into one
    
    Args:
        chr_files: List of tuples (chr_id, vcf_filepath)
        output_file: Path to merged output VCF
        keep_chr_files: If False, delete individual chr VCF files after merge
        verbose: Print detailed progress
    """
    if verbose:
        print(f"\nMerging {len(chr_files)} chromosome VCF files...")
    
    # Sort by chromosome number
    chr_files_sorted = sorted(chr_files, key=lambda x: int(x[0]) if x[0].isdigit() else x[0])
    
    with open(output_file, 'w') as out_vcf:
        header_written = False
        total_variants = 0
        
        for chr_id, vcf_file in chr_files_sorted:
            if verbose:
                print(f"  Processing chr{chr_id} from {vcf_file}")
            
            with open(vcf_file, 'r') as in_vcf:
                variant_count = 0
                
                for line in in_vcf:
                    # Write header only once (from first file)
                    if line.startswith('##'):
                        if not header_written:
                            out_vcf.write(line)
                    elif line.startswith('#CHROM'):
                        if not header_written:
                            out_vcf.write(line)
                            header_written = True
                    else:
                        # Write variant lines
                        out_vcf.write(line)
                        variant_count += 1
                
                total_variants += variant_count
                if verbose:
                    print(f"    Added {variant_count} variants")
            
            # Delete chromosome file if not keeping
            if not keep_chr_files:
                os.remove(vcf_file)
                if verbose:
                    print(f"    Deleted {vcf_file}")
    
    if verbose:
        print(f"\n  ✓ Merged VCF: {output_file}")
        print(f"  Total variants: {total_variants}")
    
    return total_variants

def main():
    parser = make_arg_parser()
    
    args = parser.parse_args()

    # Validate parameter combinations
    if args.merge_after_split and not (args.split_by_chr and args.parallel):
        parser.error("--merge-after-split requires both --split-by-chr and --parallel")

    if args.keep_chr_files and not args.merge_after_split:
        parser.error("--keep-chr-files only works with --merge-after-split")
    
    print("\n" + "="*60)
    print("AlphaSimR CSV to VCF Converter (Optimized)")
    print("="*60)
    
    snp_dict = read_map_file(args.map, verbose=args.verbose)
    snp_cols, individual_ids = read_genotype_metadata(args.input)
    
    if args.verbose:
        print(f"\nProcessing: {len(snp_cols)} SNPs, {len(individual_ids)} individuals")
    
    if args.split_by_chr:
        chromosomes = set(info[0] for info in snp_dict.values())
        chromosomes = sorted(chromosomes, key=lambda x: int(x) if x.isdigit() else x)
        
        print(f"\nWriting {len(chromosomes)} chromosome-specific VCF files...")
        
        if args.parallel and len(chromosomes) > 1:
            ncores = min(args.ncores, len(chromosomes))
            print(f"Using {ncores} cores for parallel processing")
            
            tasks = []
            for chr_id in chromosomes:
                output_file = args.output.replace('.vcf', f'_chr{chr_id}.vcf')
                tasks.append((chr_id, args.input, snp_cols, snp_dict, args, output_file))
            
            with Pool(processes=ncores) as pool:
                results = pool.map(process_chromosome, tasks)
            
            # Collect results
            chr_files = []
            for chr_id, variants, output_file in results:
                print(f"  ✓ Chr{chr_id}: {variants} variants -> {output_file}")
                chr_files.append((chr_id, output_file))
            
            # Merge if requested
            if args.merge_after_split:
                print("\n" + "="*60)
                print("Merging chromosome VCF files...")
                print("="*60)
                
                total_variants = merge_vcf_files(
                    chr_files, 
                    args.output, 
                    keep_chr_files=args.keep_chr_files,
                    verbose=args.verbose
                )
                
                if args.keep_chr_files:
                    print(f"\n  ✓ Merged file: {args.output} ({total_variants:,} total variants)")
                    print(f"  ✓ Individual chromosome files retained")
                else:
                    print(f"\n  ✓ Merged file: {args.output} ({total_variants:,} total variants)")
                    print(f"  ✓ Individual chromosome files deleted")
        else:
            # Sequential processing
            for chr_id in chromosomes:
                output_file = args.output.replace('.vcf', f'_chr{chr_id}.vcf')
                
                with open(output_file, 'w') as out_file:
                    write_vcf_header(out_file, individual_ids, verbose=args.verbose)
                    variants = write_vcf_body_streaming(out_file, args.input, snp_cols, 
                                                       snp_dict, args, chr_filter=chr_id)
                
                print(f"  ✓ {output_file} ({variants} variants)")
    else:
        # Single VCF file
        print(f"\nWriting single VCF file: {args.output}")
        
        with open(args.output, 'w') as out_file:
            write_vcf_header(out_file, individual_ids, verbose=args.verbose)
            variants = write_vcf_body_streaming(out_file, args.input, snp_cols, 
                                               snp_dict, args, chr_filter=None)
        
        print(f"  ✓ {args.output} ({variants} variants)")
    
    print("\n" + "="*60)
    print("✓ Conversion completed successfully!")
    print("="*60 + "\n")

if __name__ == '__main__':
    main()