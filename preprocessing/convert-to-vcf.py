#!/usr/bin/env python3
"""
Convert AlphaSimR CSV genotype to VCF format
Modified from convert-to-vcf.py to accept CSV input

Author: Modified for AlphaSimR pipeline
Date: 2025-01-26
"""
from __future__ import print_function
import argparse
import csv
import re
from datetime import datetime

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
    parser.add_argument("--split-by-chr",
                       action="store_true",
                       help="Create separate VCF files per chromosome")
    parser.add_argument("-v", "--verbose",
                       action="store_true",
                       default=False,
                       help="Verbose output")
    return parser

def natural_sort_key(s):
    """Sort strings with numbers naturally"""
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split('([0-9]+)', s)]

def read_map_file(map_file, verbose=False):
    """
    Read map file and create SNP dictionary
    
    Returns:
        dict: {snp_id: [chr, position, snp_id]}
    """
    if verbose:
        print(f"Reading map file: {map_file}")
    
    snp_dict = {}
    
    with open(map_file, 'r') as f:
        # Skip header
        header = f.readline()
        
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


def read_genotype_csv(csv_file, verbose=False):
    """
    Read AlphaSimR genotype CSV file
    
    Expected format:
    individual,snp_1,snp_2,...,generation,population_type,individual_id
    ind_1,0,1,2,...,5,reference,reference_gen5_ind1
    
    Returns:
        list: SNP column names (e.g., ['snp_1', 'snp_2', ...])
        dict: {individual_id: [geno_snp1, geno_snp2, ...]}
    """
    if verbose:
        print(f"Reading genotype CSV: {csv_file}")
    
    snp_cols = []
    genotypes = {}
    
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        
        # Get SNP column names (exclude metadata columns)
        all_cols = reader.fieldnames
        metadata_cols = ['ID', 'individual', 'generation', 'population_type', 'individual_id']
        snp_cols = [col for col in all_cols if col not in metadata_cols]
        
        if verbose:
            print(f"  Found {len(snp_cols)} SNP columns")
        
        # Read genotypes
        for row in reader:
            # Use individual_id if available, otherwise use individual
            ind_id = row.get('ID', row.get('individual_id', row.get('individual', 'unknown')))
            
            # Extract genotypes for SNP columns only
            geno_values = [row[snp] for snp in snp_cols]
            genotypes[ind_id] = geno_values
    
    if verbose:
        print(f"  Loaded genotypes for {len(genotypes)} individuals")
    
    return snp_cols, genotypes


def encode_genotype(geno_str):
    """
    Convert numeric genotype (0/1/2) to VCF GT format
    
    0 -> 0/0 (homozygous reference)
    1 -> 0/1 (heterozygous)
    2 -> 1/1 (homozygous alternate)
    NA/missing -> ./.
    """
    code = {'0': '0/0', '1': '0/1', '2': '1/1'}
    return code.get(str(geno_str).strip(), './.')


def write_vcf_header(out_file, individual_ids, verbose=False):
    """Write VCF header lines"""
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


def write_vcf_body(out_file, snp_cols, snp_dict, genotypes, args, chr_filter=None, verbose=False):
    """
    Write VCF body (variant lines)
    
    Args:
        chr_filter: If provided, only write SNPs from this chromosome
    """
    variants_written = 0
    
    individual_ids = sorted(genotypes.keys(), key=natural_sort_key)
    
    for snp_col in snp_cols:
        # Get SNP info from map
        if snp_col not in snp_dict:
            if verbose:
                print(f"  Warning: {snp_col} not in map file, skipping")
            continue
        
        chr_id, position, snp_id = snp_dict[snp_col]
        
        # Filter by chromosome if specified
        if chr_filter is not None and chr_id != str(chr_filter):
            continue
        
        # Build variant line
        variant_line = [
            chr_id,                  # CHROM
            position,                # POS
            snp_id,                  # ID
            args.ref,                # REF
            args.alt,                # ALT
            args.qual,               # QUAL
            args.filter,             # FILTER
            args.info,               # INFO
            args.format_field        # FORMAT
        ]
        
        # Get SNP index in genotype data
        snp_idx = snp_cols.index(snp_col)
        
        # Add genotypes for all individuals
        for ind_id in individual_ids:
            geno_value = genotypes[ind_id][snp_idx]
            gt_format = encode_genotype(geno_value)
            variant_line.append(gt_format)
        
        # Write variant line
        print('\t'.join(variant_line), file=out_file)
        variants_written += 1
    
    if verbose:
        print(f"  Wrote {variants_written} variants")
    
    return variants_written


def main():
    parser = make_arg_parser()
    args = parser.parse_args()
    
    print("\n" + "="*60)
    print("AlphaSimR CSV to VCF Converter")
    print("="*60)
    
    # Step 1: Read map file
    snp_dict = read_map_file(args.map, verbose=args.verbose)
    
    # Step 2: Read genotype CSV
    snp_cols, genotypes = read_genotype_csv(args.input, verbose=args.verbose)
    
    # Get individual IDs (sorted for consistent order)
    individual_ids = sorted(genotypes.keys(), key=natural_sort_key)
    
    if args.verbose:
        print(f"\nProcessing:")
        print(f"  {len(snp_cols)} SNPs")
        print(f"  {len(individual_ids)} individuals")
    
    # Step 3: Write VCF file(s)
    if args.split_by_chr:
        # Get unique chromosomes
        chromosomes = set(info[0] for info in snp_dict.values())
        chromosomes = sorted(chromosomes, key=lambda x: int(x) if x.isdigit() else x)
        
        print(f"\nWriting {len(chromosomes)} chromosome-specific VCF files...")
        
        for chr_id in chromosomes:
            output_file = args.output.replace('.vcf', f'_chr{chr_id}.vcf')
            
            if args.verbose:
                print(f"\nChromosome {chr_id} -> {output_file}")
            
            with open(output_file, 'w') as out_file:
                write_vcf_header(out_file, individual_ids, verbose=args.verbose)
                write_vcf_body(out_file, snp_cols, snp_dict, genotypes, args, 
                             chr_filter=chr_id, verbose=args.verbose)
            
            print(f"  ✓ {output_file}")
    
    else:
        # Single VCF file for all chromosomes
        print(f"\nWriting VCF file: {args.output}")
        
        with open(args.output, 'w') as out_file:
            write_vcf_header(out_file, individual_ids, verbose=args.verbose)
            write_vcf_body(out_file, snp_cols, snp_dict, genotypes, args, 
                         chr_filter=None, verbose=args.verbose)
        
        print(f"  ✓ {args.output}")
    
    print("\n" + "="*60)
    print("✓ Conversion completed successfully!")
    print("="*60 + "\n")


if __name__ == '__main__':
    main()