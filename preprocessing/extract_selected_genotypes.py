#!/usr/bin/env python3
"""
Extract selected SNPs from full genotype file based on microhaplotype selection

Input:
  - Full genotype file (CSV/TSV with SNP columns)
  - SNP selection file (snp_selection_detailed.csv)
  
Output:
  - Subset genotype file with only selected SNPs
  - Metadata file with reduction statistics

Author: Agus Wibowo
"""

import pandas as pd
import argparse
import os

def extract_selected_genotypes(genotype_file, snp_selection_file, 
                               output_file, sep=',', 
                               meta_cols=None,
                               verbose=True):
    """
    Extract genotypes for selected SNPs only
    
    Args:
        genotype_file: Full genotype file path
        snp_selection_file: SNP selection CSV from microhaplotype analysis
        output_file: Output genotype file path
        sep: Separator for genotype file (default: ',')
        meta_cols: Metadata columns to preserve (None = auto-detect from genotype file)
        verbose: Print progress
    """
    
    if verbose:
        print("\n" + "="*70)
        print("EXTRACTING SELECTED SNP GENOTYPES")
        print("="*70)
    
    # Step 1: Load SNP selection
    if verbose:
        print(f"\n[1/5] Loading SNP selection file...")
    
    snp_selection = pd.read_csv(snp_selection_file)
    selected_snps = snp_selection['snp_id'].unique()
    n_selected = len(selected_snps)
    
    if verbose:
        print(f"  Total entries: {len(snp_selection):,}")
        print(f"  Unique SNPs selected: {n_selected:,}")
        print(f"  Chromosomes: {snp_selection['chr'].nunique()}")
    
    # Step 2: Load genotype file (header only first)
    if verbose:
        print(f"\n[2/5] Reading genotype file header...")
    
    geno_header = pd.read_csv(genotype_file, sep=sep, nrows=0)
    all_cols = geno_header.columns.tolist()
    
    # FIX: Proper metadata detection
    if meta_cols is None:
        common_meta = ['ID', 'individual_id', 'generation', 'population_type', 
                       'population', 'sex', 'farm_id', 'sire', 'dam',
                       'father_id', 'mother_id']
        
        meta_cols = [col for col in all_cols if col in common_meta]
        
        if len(meta_cols) == 0:
            meta_cols = [col for col in all_cols 
                         if not col.startswith(('AX-', 'rs', 'SNP', 'chr', 'snp'))]
            if len(meta_cols) == 0:
                meta_cols = all_cols[:3]
        
        if verbose:
            print(f"  Auto-detected metadata: {meta_cols}")
    
    # Identify SNP columns (all non-metadata)
    snp_cols = [col for col in all_cols if col not in meta_cols]
    n_original = len(snp_cols)
    
    if verbose:
        print(f"  Total columns: {len(all_cols):,}")
        print(f"  Metadata columns: {len(meta_cols)}")
        print(f"  SNP columns: {n_original:,}")
    
    # Step 3: Validate SNP presence
    if verbose:
        print(f"\n[3/5] Validating SNP presence...")
    
    available_snps = set(snp_cols)
    selected_set = set(selected_snps)
    
    valid_snps = selected_set & available_snps
    missing_snps = selected_set - available_snps
    
    if verbose:
        print(f"  SNPs found: {len(valid_snps):,}/{n_selected:,}")
        if len(missing_snps) > 0:
            print(f"  Missing SNPs: {len(missing_snps)}")
            if len(missing_snps) <= 10:
                print(f"    {list(missing_snps)}")
    
    if len(valid_snps) == 0:
        raise ValueError("No selected SNPs found in genotype file!")
    
    # Step 4: Sort SNPs by chromosome and position
    if verbose:
        print(f"\n[4/5] Sorting SNPs by genomic order...")
    
    snp_order_df = snp_selection[['snp_id', 'chr', 'position']].drop_duplicates()
    snp_order_df = snp_order_df[snp_order_df['snp_id'].isin(valid_snps)]
    snp_order_df = snp_order_df.sort_values(['chr', 'position'])
    
    ordered_snps = snp_order_df['snp_id'].tolist()
    
    if verbose:
        print(f"  Ordered SNPs: {len(ordered_snps):,}")
    
    # Step 5: Extract genotypes
    if verbose:
        print(f"\n[5/5] Extracting genotypes...")
    
    geno_df = pd.read_csv(genotype_file, sep=sep)
    cols_to_keep = meta_cols + ordered_snps
    geno_subset = geno_df[cols_to_keep].copy()
    
    geno_subset.to_csv(output_file, index=False)
    
    # Generate metadata
    reduction_pct = (1 - len(ordered_snps) / n_original) * 100
    
    # Validation warning
    if len(ordered_snps) >= n_original * 0.95:
        print("\n  WARNING: Minimal SNP reduction detected!")
        print(f"   Reduction only {reduction_pct:.1f}%")
        print(f"   Detected metadata: {meta_cols}")
        print(f"   Re-run with: --meta-cols ID generation population_type")
    
    if verbose:
        print(f"  Input shape:  {geno_df.shape}")
        print(f"  Output shape: {geno_subset.shape}")
        print(f"  Saved to: {output_file}")

    metadata = {
        'genotype_file': os.path.basename(genotype_file),
        'snp_selection_file': os.path.basename(snp_selection_file),
        'n_individuals': len(geno_subset),
        'n_snps_original': n_original,
        'n_snps_selected': len(ordered_snps),
        'n_snps_missing': len(missing_snps),
        'reduction_pct': reduction_pct,
        'file_size_mb': os.path.getsize(output_file) / (1024**2)
    }
    
    meta_file = output_file.replace('.csv', '_metadata.csv')
    pd.DataFrame([metadata]).to_csv(meta_file, index=False)
    
    if verbose:
        print(f"\n" + "="*70)
        print("SUMMARY")
        print("="*70)
        print(f"Original SNPs:     {n_original:,}")
        print(f"Selected SNPs:     {len(ordered_snps):,}")
        print(f"Reduction:         {reduction_pct:.1f}%")
        print(f"Output file size:  {metadata['file_size_mb']:.2f} MB")
        print(f"Metadata saved:    {meta_file}")
        print("="*70)
    
    return geno_subset, metadata

def main():
    parser = argparse.ArgumentParser(
        description="Extract selected SNPs from full genotype file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument("-g", "--genotype",
                       required=True,
                       help="Full genotype file (CSV/TSV)")
    
    parser.add_argument("-s", "--selection",
                       required=True,
                       help="SNP selection file (snp_selection_detailed.csv)")
    
    parser.add_argument("-o", "--output",
                       required=True,
                       help="Output genotype file")
    
    parser.add_argument("--sep",
                       default=",",
                       help="Separator for genotype file")
    
    parser.add_argument("--meta-cols",
                       nargs='+',
                       default=None,
                       help="Metadata column names (default: auto-detect non-SNP columns)")
    
    parser.add_argument("-v", "--verbose",
                       action="store_true",
                       help="Verbose output")
    
    args = parser.parse_args()
    
    extract_selected_genotypes(
        genotype_file=args.genotype,
        snp_selection_file=args.selection,
        output_file=args.output,
        sep=args.sep,
        meta_cols=args.meta_cols,
        verbose=args.verbose
    )
    
    print("\n✓ Genotype extraction completed!")

if __name__ == '__main__':
    main()