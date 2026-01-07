#!/usr/bin/env python3
"""
Haplotype block definition using hybrid approach for PHASED data:
- Input: Phased haplotype files (output from convert-from-vcf.py)
- Three methods available:
  1. snp_count: Fixed/sliding window by SNP count (with optional LD filtering)
  2. physical: Sliding window by physical distance (bp) (with optional LD filtering)
  3. ld_haploblock: LD-based haploblock discovery (Jonas 2017) with D' threshold
     - pure mode: Select best N-SNP haplotype per haploblock
     - micro mode: Split haploblocks into fixed physical windows and create microhaplotype
- Output format compatible with get-hap-geno.py

Output block info format:
blk1 0 3
blk2 4 7
blk3 8 11

Output genotype format:
ID	hap_1_1	hap_1_1	hap_1_2	hap_1_2	hap_1_3	hap_1_3	hap_1_4	hap_1_4
ind001	1	2	3	4	1	1	2	5
ind002	2	2	4	4	1	2	5	6
ind003	1	3	3	5	2	2	2	2
ind004	3	4	5	6	1	3	6	7
ind005	2	1	4	3	2	1	5	2

Author: Agus Wibowo
"""

import argparse
import pandas as pd
import numpy as np
import os


def read_haplotype_file(hap_file, noheader=False, verbose=False):
    """
    Read phased haplotype file - SEPARATED FORMAT
    
    Format: ID hap1_snp1 hap2_snp1 hap1_snp2 hap2_snp2 hap1_snp3 hap2_snp3 ...
    
    Each individual has 2*n_snps columns (2 haplotypes per SNP)
    
    Returns:
        numpy array: (n_individuals * 2) × n_snps
        Rows are haplotypes (hap1_ind1, hap2_ind1, hap1_ind2, hap2_ind2, ...)
        Columns are SNPs
    """
    if verbose:
        print(f"    Reading: {hap_file}")
    
    try:
        with open(hap_file, 'r') as f:
            lines = f.readlines()
        
        if len(lines) == 0:
            raise ValueError("File is empty")
        
        # Skip header if present
        start_idx = 1 if not noheader else 0
        
        if start_idx >= len(lines):
            raise ValueError("No data after header")
        
        # Parse first line to get dimensions
        first_parts = lines[start_idx].strip().split()
        n_hap_cols = len(first_parts) - 1  # Exclude ID column
        n_snps = n_hap_cols // 2  # Each SNP has 2 haplotype columns
        
        if verbose:
            print(f"    Total haplotype columns: {n_hap_cols}")
            print(f"    Inferred n_snps: {n_snps}")
            print(f"    First 10 values: {first_parts[1:11]}")
        
        # Collect all haplotypes
        all_hap1 = []  # First haplotype of each individual
        all_hap2 = []  # Second haplotype of each individual
        
        for line_num, line in enumerate(lines[start_idx:], start=start_idx+1):
            parts = line.strip().split()
            
            if len(parts) < 2:
                continue
            
            # Skip ID column
            alleles = parts[1:]
            
            if len(alleles) != n_hap_cols:
                if verbose:
                    print(f"    Warning: Line {line_num} has {len(alleles)} cols, expected {n_hap_cols}")
                continue
            
            # Separate into two haplotypes per individual
            # Format: hap1_snp1, hap2_snp1, hap1_snp2, hap2_snp2, ...
            hap1 = []
            hap2 = []
            
            for i in range(0, len(alleles), 2):
                if i+1 < len(alleles):
                    try:
                        a1 = int(alleles[i])
                        a2 = int(alleles[i+1])
                        
                        # Convert to 0/1 binary if using 1/2 coding
                        # Assumes: 1=ref(0), 2=alt(1)
                        hap1.append(0 if a1 == 1 else 1)
                        hap2.append(0 if a2 == 1 else 1)
                    except ValueError:
                        hap1.append(np.nan)
                        hap2.append(np.nan)
            
            all_hap1.append(hap1)
            all_hap2.append(hap2)
        
        if len(all_hap1) == 0:
            raise ValueError("No valid haplotype data found")
        
        # Combine: interleave hap1 and hap2
        # Result: [hap1_ind1, hap2_ind1, hap1_ind2, hap2_ind2, ...]
        all_haplotypes = []
        for h1, h2 in zip(all_hap1, all_hap2):
            all_haplotypes.append(h1)
            all_haplotypes.append(h2)
        
        # Convert to numpy array
        hap_array = np.array(all_haplotypes, dtype=float)
        
        if verbose:
            print(f"    Final shape: {hap_array.shape} (n_haplotypes × n_snps)")
            unique_vals = np.unique(hap_array[~np.isnan(hap_array)])
            print(f"    Unique values: {unique_vals}")
            
            # Validation check
            expected_snps = n_snps
            actual_snps = hap_array.shape[1]
            if expected_snps != actual_snps:
                print(f"    WARNING: Expected {expected_snps} SNPs, got {actual_snps}")
        
        return hap_array
        
    except Exception as e:
        if verbose:
            print(f"    ERROR: {str(e)}")
            import traceback
            traceback.print_exc()
        return None

def calculate_haplotype_ld(hap_matrix):
    """
    Vectorized calculation of mean pairwise LD (r²) from haplotype data
    
    Args:
        hap_matrix: numpy array (n_haplotypes × n_snps)
                    Each row is a haplotype (0 or 1 at each SNP)
    
    Returns:
        float: Mean r² across all SNP pairs
    """
    n_snps = hap_matrix.shape[1]
    
    if n_snps < 2:
        return 0.0
    
    # Remove rows with any missing data
    mask = ~np.any(np.isnan(hap_matrix), axis=1)
    clean_haps = hap_matrix[mask]
    
    # Need at least 10 individuals for reliable LD
    if clean_haps.shape[0] < 10:
        return 0.0
    
    # Check for invariant SNPs (all 0 or all 1)
    snp_variances = np.var(clean_haps, axis=0)
    variant_snps = snp_variances > 0
    
    # Need at least 2 variant SNPs
    if np.sum(variant_snps) < 2:
        return 0.0
    
    # Keep only variant SNPs
    clean_haps = clean_haps[:, variant_snps]
    n_variant_snps = clean_haps.shape[1]
    
    try:
        # Vectorized correlation matrix computation
        # Transpose: n_snps × n_haplotypes (corrcoef expects features as rows)
        corr_matrix = np.corrcoef(clean_haps.T)
        
        # Extract upper triangle (all pairwise correlations)
        upper_indices = np.triu_indices(n_variant_snps, k=1)
        pairwise_r = corr_matrix[upper_indices]
        
        # Square to get r²
        r_squared = pairwise_r ** 2
        
        # Return mean LD
        return np.mean(r_squared)
        
    except:
        # Fallback if correlation fails (e.g., numerical issues)
        return 0.0


def get_chromosome_number(filename):
    """Extract chromosome number from filename"""
    import re
    match = re.search(r'chr(\d+)', filename)
    if match:
        return int(match.group(1))
    
    # Try to extract any number
    nums = ''.join(filter(str.isdigit, filename))
    return int(nums) if nums else 0

def calculate_d_prime(snp1, snp2):
    """
    Calculate D' (normalized linkage disequilibrium) between two SNPs
    D' = |D| / D_max
    where D = p_AB - p_A * p_B
    
    Args:
        snp1, snp2: arrays of 0/1 genotypes
    Returns:
        float: D' value (0-1)
    """
    # Remove missing data
    mask = ~(np.isnan(snp1) | np.isnan(snp2))
    s1 = snp1[mask]
    s2 = snp2[mask]
    
    if len(s1) < 10:
        return 0.0
    
    # Allele frequencies
    p_A = np.mean(s1)
    p_B = np.mean(s2)
    
    # Haplotype frequency (both = 1)
    p_AB = np.mean(s1 * s2)
    
    # Linkage disequilibrium
    D = p_AB - (p_A * p_B)
    
    # D_max for normalization
    if D >= 0:
        D_max = min(p_A * (1 - p_B), (1 - p_A) * p_B)
    else:
        D_max = min(p_A * p_B, (1 - p_A) * (1 - p_B))
    
    if D_max == 0:
        return 0.0
    
    return abs(D / D_max)

def build_ld_haploblocks(hap_matrix, snp_positions=None, d_prime_threshold=0.45, 
                         min_block_snps=2, verbose=False):
    """
    Build haploblocks based on consecutive D' >= threshold
    Following Jonas et al. (2017) with optional physical distance constraint
    
    Args:
        hap_matrix: numpy array (n_haplotypes × n_snps)
        snp_positions: array of SNP physical positions (bp),
        d_prime_threshold: minimum D' between consecutive SNPs (default 0.45)
        min_block_snps: minimum SNPs per haploblock (default 2)
        verbose: print progress
    
    Returns:
        list of dicts: haploblock info
    """
    n_snps = hap_matrix.shape[1]
    
    if verbose:
        print(f"    Building haploblocks with D' threshold {d_prime_threshold}")
    
    # Calculate D' for all consecutive pairs
    d_primes = []
    for i in range(n_snps - 1):
        d_prime = calculate_d_prime(hap_matrix[:, i], hap_matrix[:, i + 1])
        d_primes.append(d_prime)
    
    if verbose:
        print(f"    Mean consecutive D': {np.mean(d_primes):.3f}")
        print(f"    Pairs ≥ threshold: {np.sum(np.array(d_primes) >= d_prime_threshold)}/{len(d_primes)}")
    
    # Build haploblocks where D' >= threshold consecutively
    haploblocks = []
    current_block = [0]  # Start with first SNP
    
    for i, d_prime in enumerate(d_primes):
        
        if d_prime >= d_prime_threshold:
            # Strong LD → extend current block
            current_block.append(i + 1)
        else:
            # Weak LD or exceeded physical constraint → end current block, start new
            if len(current_block) >= min_block_snps:
                block_start = current_block[0]
                block_end = current_block[-1]
                
                block_info = {
                    'start_idx': block_start,
                    'end_idx': block_end,
                    'snp_indices': current_block.copy(),
                    'n_snps': len(current_block)
                }
                
                # Add physical span if positions available
                if snp_positions is not None:
                    block_info['physical_span'] = snp_positions[block_end] - snp_positions[block_start]
                
                haploblocks.append(block_info)
            
            current_block = [i + 1]
    
    # Add last block
    if len(current_block) >= min_block_snps:
        block_start = current_block[0]
        block_end = current_block[-1]
        
        block_info = {
            'start_idx': block_start,
            'end_idx': block_end,
            'snp_indices': current_block.copy(),
            'n_snps': len(current_block)
        }
        
        if snp_positions is not None:
            block_info['physical_span'] = snp_positions[block_end] - snp_positions[block_start]
        
        haploblocks.append(block_info)
    
    if verbose:
        sizes = [b['n_snps'] for b in haploblocks]
        print(f"    Built {len(haploblocks)} haploblocks")
        if len(sizes) > 0:
            print(f"    Block size - mean: {np.mean(sizes):.1f}, median: {np.median(sizes):.0f}, range: [{min(sizes)}, {max(sizes)}]")
            
            if snp_positions is not None and 'physical_span' in haploblocks[0]:
                spans = [b['physical_span'] for b in haploblocks]
                print(f"    Physical span (bp) - mean: {np.mean(spans):.0f}, median: {np.median(spans):.0f}, range: [{min(spans):.0f}, {max(spans):.0f}]")
    
    return haploblocks

def criterion_b_score(block_haps, aft=0.08, md=0.10):
    """
    Calculate Criterion-B score (Jonas et al. 2016)
    Lower score = better balance between allele frequency and diversity
    
    Args:
        block_haps: numpy array (n_haplotypes × n_snps)
        aft: Allele Frequency Threshold (predictable if freq >= aft)
        md: Maximum Deviation parameter
    
    Returns:
        float: criterion-B score (minimize)
    """
    n_snps = block_haps.shape[1]
    HS = 2 ** n_snps  # Theoretical maximum alleles
    
    # Convert to unique haplotype patterns
    unique_haps = []
    for hap in block_haps:
        if not np.any(np.isnan(hap)):
            unique_haps.append(tuple(hap.astype(int)))
    
    if len(unique_haps) == 0:
        return np.inf
    
    # Count frequencies
    from collections import Counter
    counts = Counter(unique_haps)
    freqs = np.array([c / len(unique_haps) for c in counts.values()])
    
    # Predictable alleles (freq >= AFT)
    n_predictable = np.sum(freqs >= aft)
    
    if n_predictable == 0:
        return np.inf
    
    pred_freqs = freqs[freqs >= aft]
    
    # Term 1: Squared deviation from ideal frequency
    ideal_freq = 1 / HS
    sd_term = np.sum((pred_freqs - ideal_freq) ** 2)
    
    # Term 2: Weighted predictable alleles
    if n_predictable == 1:
        w = 0
    else:
        w = (md * n_predictable) / (HS * (n_predictable - 1))
    
    score = sd_term - (w * n_predictable)
    
    return score

def select_best_haplotype_from_haploblock(haploblock, hap_matrix, 
                                           haplotype_size=4, 
                                           aft=0.08, md=0.10,
                                           verbose=False):
    """
    Select best 4-SNP haplotype from haploblock using Criterion-B
    
    Args:
        haploblock: dict with 'snp_indices'
        hap_matrix: full haplotype matrix
        haplotype_size: target haplotype size (default 4)
        aft, md: Criterion-B parameters
    
    Returns:
        dict: selected haplotype info
    """
    snp_indices = haploblock['snp_indices']
    n_snps = len(snp_indices)
    
    # Case 1: Block too short → extend with flanking SNPs
    if n_snps < haplotype_size:
        start = snp_indices[0]
        end = snp_indices[-1]
        
        # How many SNPs to add
        needed = haplotype_size - n_snps
        left = needed // 2
        right = needed - left
        
        # Extend boundaries
        new_start = max(0, start - left)
        new_end = min(hap_matrix.shape[1] - 1, end + right)
        
        # Get exact haplotype_size SNPs
        extended = list(range(new_start, new_end + 1))
        if len(extended) > haplotype_size:
            extended = extended[:haplotype_size]
        elif len(extended) < haplotype_size:
            # Extend further if needed
            while len(extended) < haplotype_size and extended[-1] < hap_matrix.shape[1] - 1:
                extended.append(extended[-1] + 1)
        
        selected_indices = extended[:haplotype_size]
        block_haps = hap_matrix[:, selected_indices]
        score = criterion_b_score(block_haps, aft, md)
        
        return {
            'snp_indices': selected_indices,
            'n_snps': len(selected_indices),
            'criterion_b_score': score,
            'original_block_size': n_snps,
            'selection_type': 'extended'
        }
    
    # Case 2: Block = haplotype_size → use directly
    elif n_snps == haplotype_size:
        block_haps = hap_matrix[:, snp_indices]
        score = criterion_b_score(block_haps, aft, md)
        
        return {
            'snp_indices': snp_indices,
            'n_snps': n_snps,
            'criterion_b_score': score,
            'original_block_size': n_snps,
            'selection_type': 'exact'
        }
    
    # Case 3: Block longer → try all 4-SNP windows, select best
    else:
        best_score = np.inf
        best_indices = None
        
        for start_pos in range(n_snps - haplotype_size + 1):
            window_indices = snp_indices[start_pos:start_pos + haplotype_size]
            block_haps = hap_matrix[:, window_indices]
            
            score = criterion_b_score(block_haps, aft, md)
            
            if score < best_score:
                best_score = score
                best_indices = window_indices
        
        return {
            'snp_indices': best_indices,
            'n_snps': len(best_indices),
            'criterion_b_score': best_score,
            'original_block_size': n_snps,
            'selection_type': 'selected'
        }

def split_ld_block_by_physical_window(haploblock, hap_matrix, chr_snps, 
                                       window_bp=125, min_snps=2,
                                       aft=0.08, md=0.10, verbose=False):
    """
    Split LD haploblock into non-overlapping physical windows
    Then select best microhaplotype using Criterion-B for each window
    
    Args:
        haploblock: dict with 'snp_indices'
        hap_matrix: full haplotype matrix
        chr_snps: DataFrame with SNP positions
        window_bp: physical window size in bp
        min_snps: minimum SNPs per window
        aft, md: Criterion-B parameters
    
    Returns:
        list of dicts: microhaplotype info for each valid window
    """
    snp_indices = haploblock['snp_indices']
    
    # Get positions for SNPs in this block
    positions = chr_snps.iloc[snp_indices]['Position'].values
    start_pos = positions[0]
    end_pos = positions[-1]
    block_span = end_pos - start_pos
    
    # If block shorter than window, return as single microhaplotype
    if block_span <= window_bp:
        if len(snp_indices) >= min_snps:
            block_haps = hap_matrix[:, snp_indices]
            score = criterion_b_score(block_haps, aft, md)
            
            return [{
                'snp_indices': snp_indices,
                'n_snps': len(snp_indices),
                'criterion_b_score': score,
                'physical_span': block_span,
                'split_type': 'single'
            }]
        else:
            return []
    
    # Split into non-overlapping windows
    microhaplotypes = []
    current_start_idx = 0
    
    while current_start_idx < len(snp_indices):
        window_start_pos = positions[current_start_idx]
        window_end_pos = window_start_pos + window_bp
        
        # Collect SNPs within this window
        window_snp_indices = []
        for i in range(current_start_idx, len(snp_indices)):
            if positions[i] < window_end_pos:
                window_snp_indices.append(snp_indices[i])
            else:
                break
        
        # Check if window meets min_snps criterion
        if len(window_snp_indices) >= min_snps:
            # Extract haplotypes for this window
            block_haps = hap_matrix[:, window_snp_indices]
            
            # Calculate Criterion-B score
            score = criterion_b_score(block_haps, aft, md)
            
            actual_span = positions[current_start_idx + len(window_snp_indices) - 1] - window_start_pos
            
            microhaplotypes.append({
                'snp_indices': window_snp_indices,
                'n_snps': len(window_snp_indices),
                'criterion_b_score': score,
                'physical_span': actual_span,
                'split_type': 'split'
            })
        
        # Move to next window (non-overlapping)
        current_start_idx += len(window_snp_indices)
        
        # Safety: if no SNPs added, move forward to avoid infinite loop
        if len(window_snp_indices) == 0:
            current_start_idx += 1
    
    return microhaplotypes

def define_microhaplotype_blocks_from_haplotypes(hap_files, map_file,
                                                  method='snp_count',
                                                  window_snps=4, step_snps=None,
                                                  window_bp=125, step_bp=None,
                                                  min_snps=2, max_snps=4,
                                                  min_ld_r2=None,
                                                  haplotype_type='micro',
                                                  d_prime_threshold=0.45,
                                                  aft=0.08,
                                                  md=0.10,
                                                  noheader=False,
                                                  verbose=False):
    """
    Define microhaplotype blocks from phased haplotype files
    
    Args:
        hap_files: List of phased haplotype files (per chromosome)
        map_file: SNP map file (SNPID Chr Position)
        method: Block definition method ('snp_count' or 'physical')
        window_snps: Window size in number of SNPs (for method='snp_count')
        step_snps: Step size in SNPs (for method='snp_count', None = non-overlapping)
        window_bp: Window size in base pairs (for method='physical')
        step_bp: Step size in bp (for method='physical', None = non-overlapping)
        min_snps: Minimum SNPs per block
        max_snps: Maximum SNPs per block
        min_ld_r2: Minimum mean LD threshold (None = no filter)
        d_prime_threshold: D' threshold for ld_haploblock method (default: 0.45)
        aft: Allele Frequency Threshold for Criterion-B (default: 0.08)
        md: Maximum Deviation for Criterion-B (default: 0.10)
        noheader: If True, haplotype files have no header
        verbose: Print detailed progress
    """
    
    # Validate method
    if method not in ['snp_count', 'physical', 'ld_haploblock', 'snp_count_simple']:
        raise ValueError(f"method must be 'snp_count', 'physical', 'ld_haploblock', or 'snp_count_simple', got '{method}'")
    
    # Read map file
    snp_map = pd.read_csv(map_file, sep="\t")
    if len(snp_map.columns) == 3:
        snp_map.columns = ['SNPID', 'Chr', 'Position']
    else:
        snp_map.columns = ['SNPID', 'Chr', 'Position', 'GeneticPos']
    
    # Set step size based on method
    if method == 'snp_count':
        if step_snps is None:
            step_snps = window_snps  # Non-overlapping by default
    elif method == 'physical':
        if step_bp is None:
            step_bp = window_bp  # Non-overlapping by default
    
    # Sort haplotype files by chromosome
    hap_files_sorted = sorted(hap_files, key=get_chromosome_number)
    
    all_blocks = {}
    
    for hap_file in hap_files_sorted:
        # Extract chromosome number
        chr_num = get_chromosome_number(hap_file)
        
        if verbose:
            print(f"\nProcessing chromosome {chr_num}: {hap_file}")
        
        # Get SNPs for this chromosome
        chr_snps = snp_map[snp_map['Chr'] == chr_num].reset_index(drop=True)
        chr_snps = chr_snps.sort_values('Position').reset_index(drop=True)
        
        if len(chr_snps) == 0:
            if verbose:
                print(f"  Warning: No SNPs found for chromosome {chr_num}")
            continue
        
        # Load haplotype data if LD filtering requested
        hap_matrix = None
        if min_ld_r2 is not None:
            try:
                hap_matrix = read_haplotype_file(hap_file, noheader=noheader, verbose=verbose)
                
                if hap_matrix is None:
                    if verbose:
                        print(f"  Warning: Could not read haplotypes, skipping LD filter for this chr")
                else:
                    if verbose:
                        print(f"  Loaded: {hap_matrix.shape[0]} haplotypes × {hap_matrix.shape[1]} SNPs")
                    
                    # Validate dimensions
                    if hap_matrix.shape[1] != len(chr_snps):
                        if verbose:
                            print(f"  WARNING: Dimension mismatch!")
                            print(f"    Haplotype file has {hap_matrix.shape[1]} SNPs")
                            print(f"    Map file has {len(chr_snps)} SNPs for chr {chr_num}")
                            print(f"  → Using intersection: {min(hap_matrix.shape[1], len(chr_snps))} SNPs")
                        
                        # Use intersection
                        n_common = min(hap_matrix.shape[1], len(chr_snps))
                        hap_matrix = hap_matrix[:, :n_common]
                        chr_snps = chr_snps.iloc[:n_common].reset_index(drop=True)
                        
                        if verbose:
                            print(f"  Adjusted to: {len(chr_snps)} SNPs")
            
            except Exception as e:
                if verbose:
                    print(f"  Error loading haplotypes: {str(e)}")
                hap_matrix = None
        
        if verbose:
            print(f"  Map SNPs: {len(chr_snps)}")
        
        blocks = []
        
        # === METHOD 1: SNP COUNT BASED ===
        if method == 'snp_count':
            for start_snp_idx in range(0, len(chr_snps), step_snps):
                end_snp_idx = min(start_snp_idx + window_snps - 1, len(chr_snps) - 1)
                
                # Get SNPs in this window
                window_snps_df = chr_snps.iloc[start_snp_idx:end_snp_idx+1]
                n_snps = len(window_snps_df)
                
                # Check SNP count criteria
                if min_snps <= n_snps <= max_snps:
                    
                    # Get physical positions
                    start_pos = window_snps_df['Position'].min()
                    end_pos = window_snps_df['Position'].max()
                    
                    # LD filtering if requested
                    pass_ld_filter = True
                    mean_ld = None
                    
                    if min_ld_r2 is not None and hap_matrix is not None:
                        # Extract haplotypes for SNPs in block
                        snp_indices = list(range(start_snp_idx, end_snp_idx + 1))
                        block_haps = hap_matrix[:, snp_indices]
                        
                        # Calculate mean LD
                        mean_ld = calculate_haplotype_ld(block_haps)
                        
                        # DEBUG PRINT
                        if len(blocks) < 5:
                            print(f"    Block candidate: SNP {start_snp_idx}-{end_snp_idx}, n_snps={n_snps}, LD={mean_ld:.3f}, threshold={min_ld_r2}")
                        
                        pass_ld_filter = (mean_ld >= min_ld_r2)
                    
                    if pass_ld_filter:
                        block_info = {
                            'chr': chr_num,
                            'start_pos': start_pos,
                            'end_pos': end_pos,
                            'start_idx': start_snp_idx,
                            'end_idx': end_snp_idx,
                            'n_snps': n_snps
                        }
                        
                        if mean_ld is not None:
                            block_info['mean_ld_r2'] = mean_ld
                        
                        blocks.append(block_info)
        
        # === METHOD 2: PHYSICAL DISTANCE BASED ===
        elif method == 'physical':
            start_pos = chr_snps['Position'].min()
            end_pos = chr_snps['Position'].max()
            
            # Vectorized window assignment 
            chr_snps = chr_snps.copy()  # Avoid modifying original
            chr_snps['window_id'] = ((chr_snps['Position'] - start_pos) // step_bp).astype(int)
            
            # Group SNPs by window 
            for window_id, window_snps_df in chr_snps.groupby('window_id'):
                
                # Calculate actual window boundaries
                window_start = start_pos + (window_id * step_bp)
                window_end = window_start + window_bp
                
                # Filter SNPs strictly within window
                window_snps_df = window_snps_df[
                    (window_snps_df['Position'] >= window_start) &
                    (window_snps_df['Position'] < window_end)
                ]
                
                n_snps = len(window_snps_df)
                
                # Check SNP count criteria
                if min_snps <= n_snps <= max_snps:
                    
                    # Get SNP indices 
                    start_idx = window_snps_df.index[0]
                    end_idx = window_snps_df.index[-1]
                    
                    # LD filtering if requested
                    pass_ld_filter = True
                    mean_ld = None
                    
                    if min_ld_r2 is not None and hap_matrix is not None:
                        # Extract haplotypes for SNPs in block
                        snp_indices = window_snps_df.index.tolist()
                        block_haps = hap_matrix[:, snp_indices]
                        
                        # Calculate mean LD
                        mean_ld = calculate_haplotype_ld(block_haps)
                        
                        # DEBUG PRINT
                        if len(blocks) < 5:
                            print(f"    Block candidate: pos {window_start}-{window_end}, n_snps={n_snps}, LD={mean_ld:.3f}, threshold={min_ld_r2}")
                        
                        pass_ld_filter = (mean_ld >= min_ld_r2)
                    
                    if pass_ld_filter:
                        block_info = {
                            'chr': chr_num,
                            'start_pos': window_start,
                            'end_pos': window_end,
                            'start_idx': start_idx,
                            'end_idx': end_idx,
                            'n_snps': n_snps
                        }
                        
                        if mean_ld is not None:
                            block_info['mean_ld_r2'] = mean_ld
                        
                        blocks.append(block_info)
        
        # === METHOD 3: SNP COUNT SIMPLE (block-by-snp.py logic) ===
        elif method == 'snp_count_simple':
            # Simple fixed-size blocks like block-by-snp.py
            # No LD filtering, no overlap
            n_snps = len(chr_snps)
            sd = window_snps  # block size
            rem = n_snps % sd
            blks = (n_snps - rem) // sd
            
            # Create full blocks
            for i in range(0, blks * sd, sd):
                block_info = {
                    'chr': chr_num,
                    'start_pos': chr_snps.iloc[i]['Position'],
                    'end_pos': chr_snps.iloc[i + sd - 1]['Position'],
                    'start_idx': i,
                    'end_idx': i + sd - 1,
                    'n_snps': sd
                }
                blocks.append(block_info)
            
            # Last partial block if remainder exists
            if rem != 0:
                start_idx = n_snps - rem
                block_info = {
                    'chr': chr_num,
                    'start_pos': chr_snps.iloc[start_idx]['Position'],
                    'end_pos': chr_snps.iloc[n_snps - 1]['Position'],
                    'start_idx': start_idx,
                    'end_idx': n_snps - 1,
                    'n_snps': rem
                }
                blocks.append(block_info)
        
        # === METHOD 4: LD-BASED HAPLOBLOCKS (Jonas 2017) ===
        elif method == 'ld_haploblock':
            # Must load haplotype data
            if hap_matrix is None:
                hap_matrix = read_haplotype_file(hap_file, noheader=noheader, verbose=verbose)
            
            if hap_matrix is None:
                if verbose:
                    print(f"  ERROR: Cannot read haplotypes for LD haploblock method")
                continue

            # Get SNP positions for physical constraint
            snp_positions = chr_snps['Position'].values

            # STAGE 1: Build haploblocks based on consecutive D'
            if verbose:
                print(f"  STAGE 1: Building LD haploblocks...")
            
            # Build haploblocks based on consecutive D'
            haploblocks = build_ld_haploblocks(
                hap_matrix,
                snp_positions=snp_positions,
                d_prime_threshold=d_prime_threshold,  # USE PARAMETER
                min_block_snps=min_snps,
                verbose=verbose
            )
            
            if verbose:
                print(f"  Built {len(haploblocks)} LD haploblocks")
                if haplotype_type == 'pure':
                    print(f"  STAGE 2: Selecting best {window_snps}-SNP haplotypes from each block...")
                else:
                    print(f"  STAGE 2: Splitting blocks into {window_bp}bp windows (microhaplotypes)...")

            
            # STAGE 2: Split each LD block into physical windows + select via Criterion-B
            n_output = 0
            for haploblock in haploblocks:
                # Choice 1: Pure haplotype (Jonas 2017 original)
                if haplotype_type == 'pure':
                    selected = select_best_haplotype_from_haploblock(
                        haploblock,
                        hap_matrix,
                        haplotype_size=window_snps,
                        aft=aft,
                        md=md,
                        verbose=False
                    )
                    
                    sel_indices = selected['snp_indices']
                    
                    block_info = {
                        'chr': chr_num,
                        'start_pos': chr_snps.iloc[sel_indices[0]]['Position'],
                        'end_pos': chr_snps.iloc[sel_indices[-1]]['Position'],
                        'start_idx': sel_indices[0],
                        'end_idx': sel_indices[-1],
                        'n_snps': selected['n_snps'],
                        'criterion_b_score': selected['criterion_b_score'],
                        'original_block_size': selected['original_block_size'],
                        'selection_type': selected['selection_type']
                    }
                    
                    blocks.append(block_info)
                    n_output += 1
                
                # Choice 2: Microhaplotype (physical window split)
                else:  # haplotype_type == 'micro'
                    microhaplotypes = split_ld_block_by_physical_window(
                        haploblock,
                        hap_matrix,
                        chr_snps,
                        window_bp=window_bp,
                        min_snps=min_snps,
                        aft=aft,
                        md=md,
                        verbose=False
                    )
                    
                    for micro in microhaplotypes:
                        sel_indices = micro['snp_indices']
                        
                        block_info = {
                            'chr': chr_num,
                            'start_pos': chr_snps.iloc[sel_indices[0]]['Position'],
                            'end_pos': chr_snps.iloc[sel_indices[-1]]['Position'],
                            'start_idx': sel_indices[0],
                            'end_idx': sel_indices[-1],
                            'n_snps': micro['n_snps'],
                            'criterion_b_score': micro['criterion_b_score'],
                            'physical_span': micro['physical_span'],
                            'split_type': micro['split_type']
                        }
                        
                        blocks.append(block_info)
                        n_output += 1

            if verbose:
                output_label = "haplotypes" if haplotype_type == 'pure' else "microhaplotypes"
                print(f"  Generated {n_output} {output_label} from {len(haploblocks)} LD blocks")
                
                # Statistics
                scores = [b['criterion_b_score'] for b in blocks if not np.isinf(b['criterion_b_score'])]
                if len(scores) > 0:
                    print(f"  Mean Criterion-B score: {np.mean(scores):.6f}")
                
                if haplotype_type == 'pure':
                    sel_types = [b['selection_type'] for b in blocks]
                    from collections import Counter
                    type_counts = Counter(sel_types)
                    print(f"  Selection types: {dict(type_counts)}")
                else:
                    spans = [b['physical_span'] for b in blocks]
                    print(f"  Physical span - mean: {np.mean(spans):.1f} bp, median: {np.median(spans):.0f} bp")
                    
                    split_types = [b['split_type'] for b in blocks]
                    from collections import Counter
                    type_counts = Counter(split_types)
                    print(f"  Split types: {dict(type_counts)}")
                    
        if verbose:
            print(f"  Defined blocks: {len(blocks)}")
            if min_ld_r2 is not None and len(blocks) > 0:
                avg_ld = np.mean([b['mean_ld_r2'] for b in blocks if 'mean_ld_r2' in b])
                print(f"  Mean LD (r²): {avg_ld:.3f}")
        
        all_blocks[chr_num] = blocks
    
    return all_blocks

def deduplicate_blocks(all_blocks, verbose=False):
    """
    Remove duplicate blocks with identical SNP ranges
    Following Jonas et al. (2017): "When adjacent haploblocks selected 
    the same 4-SNP window, only one was retained"
    
    Args:
        all_blocks: dict of chromosome -> list of blocks
        verbose: print deduplication stats
    
    Returns:
        dict: deduplicated blocks
    """
    total_removed = 0
    
    for chr_num in all_blocks:
        blocks = all_blocks[chr_num]
        
        if len(blocks) == 0:
            continue
        
        seen_ranges = set()
        unique_blocks = []
        duplicates = 0
        
        for block in blocks:
            # Create tuple of SNP indices as unique identifier
            snp_range = tuple(range(block['start_idx'], block['end_idx'] + 1))
            
            if snp_range not in seen_ranges:
                seen_ranges.add(snp_range)
                unique_blocks.append(block)
            else:
                duplicates += 1
        
        all_blocks[chr_num] = unique_blocks
        total_removed += duplicates
        
        if verbose and duplicates > 0:
            print(f"  Chr {chr_num}: removed {duplicates} duplicate block(s) "
                  f"({len(blocks)} → {len(unique_blocks)})")
    
    if verbose and total_removed > 0:
        print(f"\nTotal duplicates removed: {total_removed}")
    
    return all_blocks

def calculate_genome_coverage(all_blocks, map_file, verbose=False):
    """
    Calculate detailed genome coverage statistics
    
    Returns dict with coverage metrics
    """
    # Read map file
    snp_map = pd.read_csv(map_file, sep="\t")
    if len(snp_map.columns) == 3:
        snp_map.columns = ['SNPID', 'Chr', 'Position']
    
    coverage_stats = {
        'total_chromosomes': len(all_blocks),
        'total_blocks': 0,
        'total_snps_in_blocks': 0,
        'total_genome_length_bp': 0,
        'covered_length_bp': 0,
        'coverage_pct': 0.0,
        'mean_block_spacing_kb': 0.0,
        'chromosomes': []
    }
    
    all_spacings = []
    
    for chr_num in sorted(all_blocks.keys()):
        blocks = all_blocks[chr_num]
        
        if len(blocks) == 0:
            continue
        
        # Get chromosome info from map
        chr_snps = snp_map[snp_map['Chr'] == chr_num]
        chr_length = chr_snps['Position'].max() - chr_snps['Position'].min()
        
        # Calculate covered length (sum of block spans)
        covered = sum(b['end_pos'] - b['start_pos'] for b in blocks)
        
        # Calculate block spacing
        if len(blocks) > 1:
            positions = sorted([b['start_pos'] for b in blocks])
            spacings = [positions[i+1] - positions[i] for i in range(len(positions)-1)]
            all_spacings.extend(spacings)
            mean_spacing = np.mean(spacings) / 1000  # Convert to kb
        else:
            mean_spacing = 0
        
        chr_coverage = (covered / chr_length * 100) if chr_length > 0 else 0
        
        coverage_stats['chromosomes'].append({
            'chr': chr_num,
            'n_blocks': len(blocks),
            'length_bp': chr_length,
            'covered_bp': covered,
            'coverage_pct': chr_coverage,
            'mean_spacing_kb': mean_spacing
        })
        
        # Accumulate totals
        coverage_stats['total_blocks'] += len(blocks)
        coverage_stats['total_snps_in_blocks'] += sum(b['n_snps'] for b in blocks)
        coverage_stats['total_genome_length_bp'] += chr_length
        coverage_stats['covered_length_bp'] += covered
    
    # Calculate overall statistics
    if coverage_stats['total_genome_length_bp'] > 0:
        coverage_stats['coverage_pct'] = (
            coverage_stats['covered_length_bp'] / 
            coverage_stats['total_genome_length_bp'] * 100
        )
    
    if len(all_spacings) > 0:
        coverage_stats['mean_block_spacing_kb'] = np.mean(all_spacings) / 1000
        coverage_stats['median_block_spacing_kb'] = np.median(all_spacings) / 1000
        coverage_stats['min_block_spacing_kb'] = np.min(all_spacings) / 1000
        coverage_stats['max_block_spacing_kb'] = np.max(all_spacings) / 1000
    
    # Print summary if verbose
    if verbose:
        print(f"\n" + "="*70)
        print("GENOME COVERAGE ANALYSIS")
        print("="*70)
        print(f"Total genome length:     {coverage_stats['total_genome_length_bp']:,} bp")
        print(f"Covered by blocks:       {coverage_stats['covered_length_bp']:,} bp")
        print(f"Coverage:                {coverage_stats['coverage_pct']:.2f}%")
        print(f"Total blocks:            {coverage_stats['total_blocks']:,}")
        print(f"Total SNPs in blocks:    {coverage_stats['total_snps_in_blocks']:,}")
        
        if coverage_stats['mean_block_spacing_kb'] > 0:
            print(f"\nBlock spacing statistics:")
            print(f"  Mean:                  {coverage_stats['mean_block_spacing_kb']:.1f} kb")
            print(f"  Median:                {coverage_stats['median_block_spacing_kb']:.1f} kb")
            print(f"  Range:                 [{coverage_stats['min_block_spacing_kb']:.1f}, {coverage_stats['max_block_spacing_kb']:.1f}] kb")
        
        print(f"\nPer-chromosome coverage:")
        print(f"  {'Chr':<5} {'Blocks':<8} {'Length (Mb)':<12} {'Covered (kb)':<14} {'Coverage %':<12}")
        print(f"  {'-'*5} {'-'*8} {'-'*12} {'-'*14} {'-'*12}")
        
        for chr_stat in coverage_stats['chromosomes']:
            print(f"  {chr_stat['chr']:<5} "
                  f"{chr_stat['n_blocks']:<8} "
                  f"{chr_stat['length_bp']/1e6:<12.2f} "
                  f"{chr_stat['covered_bp']/1e3:<14.1f} "
                  f"{chr_stat['coverage_pct']:<12.2f}")
        
        print("="*70)
    
    return coverage_stats

def calculate_allele_metrics(hap_matrix, block_info):
    """
    Calculate allele frequency and diversity metrics for a block
    
    Returns dict with metrics
    """
    snp_indices = list(range(block_info['start_idx'], block_info['end_idx'] + 1))
    block_haps = hap_matrix[:, snp_indices]
    
    # Get unique haplotype patterns
    unique_haps = []
    for hap in block_haps:
        if not np.any(np.isnan(hap)):
            unique_haps.append(tuple(hap.astype(int)))
    
    if len(unique_haps) == 0:
        return None
    
    from collections import Counter
    counts = Counter(unique_haps)
    freqs = np.array([c / len(unique_haps) for c in counts.values()])
    
    # Metrics
    n_alleles = len(freqs)
    n_predictable = np.sum(freqs >= 0.08)  # AFT threshold
    min_freq = np.min(freqs)
    max_freq = np.max(freqs)
    mean_freq = np.mean(freqs)
    
    # Shannon entropy
    entropy = -np.sum(freqs * np.log2(freqs + 1e-10))
    
    # Expected heterozygosity (He)
    he = 1 - np.sum(freqs ** 2)
    
    # PIC (Polymorphic Information Content)
    pic = 1 - np.sum(freqs ** 2) - np.sum([2 * freqs[i]**2 * freqs[j]**2 
                                             for i in range(len(freqs)) 
                                             for j in range(i+1, len(freqs))])
    
    # Effective number of alleles
    ne = 1 / np.sum(freqs ** 2)
    
    # Rare alleles proportion (freq < 0.05)
    rare_prop = np.sum(freqs < 0.05) / n_alleles
    
    return {
        'n_alleles': n_alleles,
        'n_predictable': n_predictable,
        'min_freq': min_freq,
        'max_freq': max_freq,
        'mean_freq': mean_freq,
        'freq_entropy': entropy,
        'expected_heterozygosity': he,
        'pic': pic,
        'effective_alleles': ne,
        'rare_alleles_prop': rare_prop
    }

def generate_haplotype_genotypes(all_blocks, hap_files, output_dir,
                                  missing_code=None, missing_value="-9999",
                                  noheader=False, verbose=False):
    """
    Generate haplotype genotype files from block definitions
    Compatible with masbayes input format
    
    Args:
        all_blocks: dict of chromosome -> list of blocks
        hap_files: list of haplotype files
        output_dir: output directory
        missing_code: coding for missing alleles (e.g., "9")
        missing_value: value for blocks with missing alleles
        noheader: if True, haplotype files have no header
        verbose: print progress
    """
    os.makedirs(output_dir, exist_ok=True)
    
    if verbose:
        print(f"\nGenerating haplotype genotypes...")
    
    # Sort haplotype files by chromosome
    hap_files_sorted = sorted(hap_files, key=get_chromosome_number)
    
    for hap_file in hap_files_sorted:
        chr_num = get_chromosome_number(hap_file)
        
        if chr_num not in all_blocks or len(all_blocks[chr_num]) == 0:
            if verbose:
                print(f"  Chr {chr_num}: No blocks defined, skipping")
            continue
        
        blocks = all_blocks[chr_num]
        
        if verbose:
            print(f"  Chr {chr_num}: Processing {len(blocks)} blocks")
        
        # Read haplotype file
        try:
            with open(hap_file, 'r') as f:
                lines = f.readlines()
            
            # Skip header if present
            start_idx = 1 if not noheader else 0
            
            # Parse haplotypes
            haps = []
            for line in lines[start_idx:]:
                parts = line.strip().split()
                if len(parts) >= 2:
                    haps.append(parts)
            
        except Exception as e:
            if verbose:
                print(f"  Error reading {hap_file}: {e}")
            continue
        
        # Encode haplotypes per block
        codes = []
        for block in blocks:
            label = 1
            hap_dict = {}
            
            for hh in haps:
                # Reconstruct full haplotype string
                if len(hh) == 2:
                    h = hh[1]
                else:
                    h = ''.join(hh[1:])
                
                # Extract substring for this block
                # Format: hap1_snp1 hap2_snp1 hap1_snp2 hap2_snp2 ...
                # So SNP index i corresponds to positions [i*2, i*2+1]
                start_pos = block['start_idx'] * 2
                end_pos = block['end_idx'] * 2 + 2
                h_part = h[start_pos:end_pos]
                
                # Split into two haplotypes (even/odd positions)
                parents = [h_part[::2], h_part[1::2]]
                
                for p in parents:
                    # Check for missing alleles
                    if missing_code is not None and missing_code in p:
                        hap_dict[p] = missing_value
                    elif p not in hap_dict:
                        hap_dict[p] = str(label)
                        label += 1
            
            codes.append(hap_dict)
        
        # Write genotype file
        outfile = os.path.join(output_dir, f'hap_geno_{chr_num}')
        
        try:
            with open(outfile, 'w') as o:
                # Header
                header = ['ID'] + [f"hap_{chr_num}_{b}\thap_{chr_num}_{b}" 
                                   for b in range(1, len(blocks) + 1)]
                o.write('\t'.join(header) + '\n')
                
                # Data rows
                for hh in haps:
                    if len(hh) == 2:
                        ind, h = hh
                    else:
                        ind = hh[0]
                        h = ''.join(hh[1:])
                    
                    output = [ind]
                    
                    for c, block in enumerate(blocks):
                        start_pos = block['start_idx'] * 2
                        end_pos = block['end_idx'] * 2 + 2
                        h_part = h[start_pos:end_pos]
                        
                        parents = [h_part[::2], h_part[1::2]]
                        
                        for p in parents:
                            output.append(codes[c][p])
                    
                    o.write('\t'.join(output) + '\n')
            
            if verbose:
                print(f"    Written: {outfile}")
        
        except Exception as e:
            if verbose:
                print(f"  Error writing {outfile}: {e}")

def write_snp_selection_file(all_blocks, map_file, output_dir, verbose=False):
    """
    Generate detailed SNP selection file showing exact SNPs in each block
    
    Output: {output_dir}/stats/snp_selection_detailed.csv
    """
    stats_dir = os.path.join(output_dir, 'stats')
    os.makedirs(stats_dir, exist_ok=True)
    
    # Read map file
    snp_map = pd.read_csv(map_file, sep="\t")
    if len(snp_map.columns) == 3:
        snp_map.columns = ['SNPID', 'Chr', 'Position']
    else:
        snp_map.columns = ['SNPID', 'Chr', 'Position', 'GeneticPos']
    
    snp_data = []
    block_counter = 1
    
    for chr_num in sorted(all_blocks.keys()):
        blocks = all_blocks[chr_num]
        chr_snps = snp_map[snp_map['Chr'] == chr_num].reset_index(drop=True)
        
        for block in blocks:
            block_id = f"blk{block_counter}"
            
            # Get all SNP indices in this block
            for idx in range(block['start_idx'], block['end_idx'] + 1):
                snp_info = chr_snps.iloc[idx]
                
                snp_data.append({
                    'block_id': block_id,
                    'chr': chr_num,
                    'snp_index': idx,
                    'snp_id': snp_info['SNPID'],
                    'position': snp_info['Position']
                })
            
            block_counter += 1
    
    snp_df = pd.DataFrame(snp_data)
    snp_file = os.path.join(stats_dir, 'snp_selection_detailed.csv')
    snp_df.to_csv(snp_file, index=False)
    
    if verbose:
        print(f"  ✓ SNP selection: {snp_file}")
        print(f"    Total SNPs selected: {len(snp_df):,}")
        print(f"    Unique SNPs: {snp_df['snp_id'].nunique():,}")
    
    return snp_file

def write_csv_outputs(all_blocks, hap_files, map_file, output_dir, 
                      method_name, noheader=False, verbose=False):
    """
    Generate comprehensive CSV outputs for analysis and publication
    
    Outputs saved to {output_dir}/stats/
    """
    
    stats_dir = os.path.join(output_dir, 'stats')
    os.makedirs(stats_dir, exist_ok=True)
    
    if verbose:
        print(f"\nGenerating analysis CSV files...")
        print(f"  Output directory: {stats_dir}/")
    
    # Read map file
    snp_map = pd.read_csv(map_file, sep="\t")
    if len(snp_map.columns) == 3:
        snp_map.columns = ['SNPID', 'Chr', 'Position']
    
    # === 1. MICROHAPLOTYPE COORDINATES CSV ===
    coords_data = []
    block_counter = 1
    
    for chr_num in sorted(all_blocks.keys()):
        blocks = all_blocks[chr_num]
        
        for block in blocks:
            row = {
                'block_id': f"blk{block_counter}",
                'chr': block['chr'],
                'start_pos': block['start_pos'],
                'end_pos': block['end_pos'],
                'n_snps': block['n_snps'],
                'physical_span_bp': block['end_pos'] - block['start_pos']
            }
            
            # Optional fields
            if 'mean_ld_r2' in block:
                row['mean_ld_r2'] = block['mean_ld_r2']
            if 'criterion_b_score' in block:
                row['criterion_b_score'] = block['criterion_b_score']
            if 'selection_type' in block:
                row['selection_type'] = block['selection_type']
            if 'original_block_size' in block:
                row['original_block_size'] = block['original_block_size']
            
            coords_data.append(row)
            block_counter += 1
    
    coords_df = pd.DataFrame(coords_data)
    coords_file = os.path.join(stats_dir, f'microhaplotype_coordinates.csv')
    coords_df.to_csv(coords_file, index=False)
    if verbose:
        print(f"  ✓ Coordinates: {coords_file}")
    
    # === 2. ALLELE FREQUENCY & DIVERSITY METRICS ===
    # Need to reload haplotype data
    allele_data = []
    
    for hap_file in sorted(hap_files, key=get_chromosome_number):
        chr_num = get_chromosome_number(hap_file)
        
        if chr_num not in all_blocks or len(all_blocks[chr_num]) == 0:
            continue
        
        hap_matrix = read_haplotype_file(hap_file, noheader=noheader, verbose=False)
        
        if hap_matrix is None:
            continue
        
        for i, block in enumerate(all_blocks[chr_num]):
            metrics = calculate_allele_metrics(hap_matrix, block)
            
            if metrics is not None:
                row = {
                    'block_id': coords_data[sum(len(all_blocks[c]) for c in sorted(all_blocks.keys()) if c < chr_num) + i]['block_id'],
                    'chr': chr_num,
                    **metrics
                }
                allele_data.append(row)
    
    allele_df = pd.DataFrame(allele_data)
    allele_file = os.path.join(stats_dir, f'allele_frequency_metrics.csv')
    allele_df.to_csv(allele_file, index=False)
    if verbose:
        print(f"  ✓ Allele metrics: {allele_file}")
    
    # === 3. CHROMOSOME SUMMARY STATS ===
    chr_summary = []
    
    for chr_num in sorted(all_blocks.keys()):
        blocks = all_blocks[chr_num]
        
        if len(blocks) == 0:
            continue
        
        chr_snps = snp_map[snp_map['Chr'] == chr_num]
        total_snps = len(chr_snps)
        chr_length = chr_snps['Position'].max() - chr_snps['Position'].min()
        
        # Calculate coverage
        covered_bp = sum(b['end_pos'] - b['start_pos'] for b in blocks)
        coverage_pct = (covered_bp / chr_length * 100) if chr_length > 0 else 0
        
        # Calculate spacing
        positions = sorted([b['start_pos'] for b in blocks])
        if len(positions) > 1:
            spacings = [positions[i+1] - positions[i] for i in range(len(positions)-1)]
            mean_spacing_kb = np.mean(spacings) / 1000
        else:
            mean_spacing_kb = 0
        
        blocks_per_mb = len(blocks) / (chr_length / 1e6) if chr_length > 0 else 0
        
        chr_summary.append({
            'chr': chr_num,
            'n_blocks': len(blocks),
            'total_snps': total_snps,
            'chr_length_bp': chr_length,
            'coverage_bp': covered_bp,
            'coverage_pct': coverage_pct,
            'mean_spacing_kb': mean_spacing_kb,
            'blocks_per_mb': blocks_per_mb
        })
    
    chr_df = pd.DataFrame(chr_summary)
    chr_file = os.path.join(stats_dir, f'chromosome_summary.csv')
    chr_df.to_csv(chr_file, index=False)
    if verbose:
        print(f"  ✓ Chromosome summary: {chr_file}")
    
    # === 4. ALLELE FREQUENCY SPECTRUM ===
    if len(allele_data) > 0:
        bins = [0, 0.01, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 1.0]
        bin_labels = ['0.00-0.01', '0.01-0.05', '0.05-0.10', '0.10-0.15', 
                    '0.15-0.20', '0.20-0.25', '0.25-0.30', '0.30-0.40', 
                    '0.40-0.50', '0.50-1.00']
        
        # Collect ALL allele frequencies from all blocks
        all_freqs = []
        
        for hap_file in sorted(hap_files, key=get_chromosome_number):
            chr_num = get_chromosome_number(hap_file)
            
            if chr_num not in all_blocks or len(all_blocks[chr_num]) == 0:
                continue
            
            hap_matrix = read_haplotype_file(hap_file, noheader=noheader, verbose=False)
            
            if hap_matrix is None:
                continue
            
            for block in all_blocks[chr_num]:
                snp_indices = list(range(block['start_idx'], block['end_idx'] + 1))
                block_haps = hap_matrix[:, snp_indices]
                
                # Get unique haplotype patterns
                unique_haps = []
                for hap in block_haps:
                    if not np.any(np.isnan(hap)):
                        unique_haps.append(tuple(hap.astype(int)))
                
                if len(unique_haps) > 0:
                    from collections import Counter
                    counts = Counter(unique_haps)
                    freqs = [c / len(unique_haps) for c in counts.values()]
                    all_freqs.extend(freqs)  # ADD TO LIST
        
        # Bin the frequencies
        hist, _ = np.histogram(all_freqs, bins=bins)
        
        spectrum_data = []
        for i, label in enumerate(bin_labels):
            spectrum_data.append({
                'freq_bin': label,
                'method': method_name,
                'count': int(hist[i])
            })
        
        spectrum_df = pd.DataFrame(spectrum_data)
        spectrum_file = os.path.join(stats_dir, f'allele_frequency_spectrum.csv')
        spectrum_df.to_csv(spectrum_file, index=False)
        if verbose:
            print(f"  ✓ Frequency spectrum: {spectrum_file}")
    
    # === 5. OVERALL SUMMARY ===
    summary = {
        'method': method_name,
        'total_chromosomes': len(all_blocks),
        'total_blocks': sum(len(blocks) for blocks in all_blocks.values()),
        'total_snps_covered': sum(b['n_snps'] for blocks in all_blocks.values() for b in blocks),
        'mean_snps_per_block': np.mean([b['n_snps'] for blocks in all_blocks.values() for b in blocks]),
        'mean_physical_span': np.mean([b['end_pos'] - b['start_pos'] for blocks in all_blocks.values() for b in blocks]),
        'genome_coverage_pct': np.mean([s['coverage_pct'] for s in chr_summary])
    }
    
    if len(allele_data) > 0:
        summary['mean_alleles_per_block'] = allele_df['n_alleles'].mean()
        summary['mean_expected_het'] = allele_df['expected_heterozygosity'].mean()
        summary['mean_pic'] = allele_df['pic'].mean()
    
    summary_df = pd.DataFrame([summary])
    summary_file = os.path.join(stats_dir, f'method_summary.csv')
    summary_df.to_csv(summary_file, index=False)
    if verbose:
        print(f"  ✓ Method summary: {summary_file}")
    
    if verbose:
        print(f"\n  Generated {5 if len(allele_data) > 0 else 4} CSV files in {stats_dir}/")
    
    return stats_dir

def write_output_files(all_blocks, output_dir, verbose=False):
    """
    Write hap_block and hap_info files per chromosome
    Compatible with get-hap-geno.py format
    """
    
    os.makedirs(output_dir, exist_ok=True)
    
    block_counter = 1
    total_blocks = 0
    
    for chrom in sorted(all_blocks.keys()):
        blocks = all_blocks[chrom]
        
        if len(blocks) == 0:
            continue
        
        # Output files
        hap_block_file = os.path.join(output_dir, f"hap_block_{chrom}")
        hap_info_file = os.path.join(output_dir, f"hap_block_info_{chrom}")
        
        with open(hap_block_file, 'w') as hb, open(hap_info_file, 'w') as hi:
            for block in blocks:
                block_id = f"blk{block_counter}"
                start_idx = block['start_idx']
                end_idx = block['end_idx']
                
                # hap_block: blk1 start_idx end_idx
                hb.write(f"{block_id}\t {start_idx} {end_idx}\n")
                
                # hap_info: blk1 idx1 idx2 idx3 ...
                indices = '\t'.join(str(i) for i in range(start_idx, end_idx + 1))
                hi.write(f"{block_id}\t {indices}\n")
                
                block_counter += 1
        
        total_blocks += len(blocks)
        
        if verbose:
            print(f"  Chr {chrom}: {len(blocks)} blocks → {os.path.basename(hap_block_file)}")
    
    return total_blocks


def main():
    parser = argparse.ArgumentParser(
        description="Microhaplotype block definition from PHASED haplotypes - Hybrid approach",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument("-i", "--input",
                       nargs='+',
                       required=True,
                       help="Phased haplotype files (e.g., output/ref/haplo/chr*)")
    
    parser.add_argument("-m", "--map",
                       required=True,
                       help="SNP map file (SNPID Chr Position)")
    
    parser.add_argument("--method",
                       choices=['snp_count', 'physical', 'ld_haploblock', 'snp_count_simple'],
                       default='snp_count',
                       help="Block definition method: "
                           "'snp_count' = fixed sliding window with LD filter, "
                           "'snp_count_simple' = fixed blocks (block-by-snp.py), "
                           "'physical' = physical distance window, "
                           "'ld_haploblock' = LD-based haploblocks (Jonas 2017)")
    
    parser.add_argument("--window-bp",
                       type=int,
                       default=125,
                       help="Window size in base pairs (for method='physical' or 'ld_haploblock')")
    
    parser.add_argument("--haplotype-type",
                       choices=['pure', 'micro'],
                       default='pure',
                       help="For ld_haploblock method - "
                        "    pure: Select best N-SNP haplotype per block (Jonas 2017)\n"
                        "    | micro: Split blocks into physical windows (microhaplotype)")
    
    parser.add_argument("--step-bp",
                       type=int,
                       default=None,
                       help="Step size in bp (for method='physical')")
    
    parser.add_argument("-w", "--window",
                       type=int,
                       default=4,
                       help="Window size in SNPs (for method='snp_count') or bp (deprecated, use --window-bp)")
    
    parser.add_argument("-s", "--step",
                       type=int,
                       default=None,
                       help="Step size in number of SNPs (default: same as window = non-overlapping)")
    
    parser.add_argument("--min-snps",
                       type=int,
                       default=2,
                       help="Minimum SNPs per block")
    
    parser.add_argument("--max-snps",
                       type=int,
                       default=4,
                       help="Maximum SNPs per block")
    
    parser.add_argument("--d-prime-threshold",
                       type=float,
                       default=0.45,
                       help="D' threshold for ld_haploblock method (Jonas 2017: 0.45)")
    
    parser.add_argument("--aft",
                       type=float,
                       default=0.08,
                       help="Allele Frequency Threshold for Criterion-B (default: 0.08)")
    
    parser.add_argument("--md",
                       type=float,
                       default=0.10,
                       help="Maximum Deviation for Criterion-B (default: 0.10)")
    
    parser.add_argument("--no-dedup",
                      action="store_true",
                      help="Skip deduplication of identical SNP ranges (not recommended)")
    
    parser.add_argument("--min-ld",
                       type=float,
                       default=None,
                       help="Minimum mean LD (r²) threshold (e.g., 0.3)")
    
    parser.add_argument("--noheader",
                       action="store_true",
                       help="Haplotype files have no header")
    
    parser.add_argument("-o", "--output",
                       default="hap_info_microhap",
                       help="Output directory name")
    
    parser.add_argument("--generate-genotypes",
                       type=str,
                       default=None,
                       metavar="PREFIX",
                       help="Generate haplotype genotype files with specified output prefix (e.g., 'mh_geno_ld')")
    
    parser.add_argument("--missing",
                       type=str,
                       default=None,
                       help="Coding for missing alleles (e.g., '9'), blocks with this will be set to missing-value")
    
    parser.add_argument("--missing-value",
                       type=str,
                       default="-9999",
                       help="Value to use for blocks with missing alleles")
    
    parser.add_argument("-v", "--verbose",
                       action="store_true",
                       help="Verbose output")
    
    args = parser.parse_args()
    
    # Print configuration
    print("\n" + "="*70)
    print("MICROHAPLOTYPE BLOCK DEFINITION - HYBRID (PHASED DATA)")
    print("="*70)
    print(f"Haplotype files:   {len(args.input)} files")
    print(f"Map file:          {args.map}")

    if args.method == 'ld_haploblock':
        print(f"Method:            LD-based haploblocks (Jonas 2017)")
        print(f"D' threshold:      {args.d_prime_threshold}")
        if args.haplotype_type == 'pure':
            print(f"Haplotype size:    {args.window} SNPs (best selection)")
        else:
            print(f"Physical window:   {args.window_bp} bp")
        print(f"Criterion-B AFT:   {args.aft}")
        print(f"Criterion-B MD:    {args.md}")
    else:
        print(f"Method:            {args.method}")
        print(f"Window size:       {args.window} SNPs")
        step_desc = f"{args.step} SNPs (overlapping)" if args.step else f"{args.window} SNPs (non-overlapping)"
        print(f"Step size:         {step_desc}")
        print(f"SNP range:         {args.min_snps} - {args.max_snps}")
        if args.min_ld:
            print(f"Min LD (r²):       {args.min_ld}")

    print(f"Output directory:  {args.output}")
    print("="*70 + "\n")
    
    # Define blocks
    all_blocks = define_microhaplotype_blocks_from_haplotypes(
        hap_files=args.input,
        map_file=args.map,
        method=args.method,
        window_snps=args.window,
        step_snps=args.step,
        window_bp=args.window_bp,
        step_bp=args.step_bp,
        min_snps=args.min_snps,
        max_snps=args.max_snps,
        min_ld_r2=args.min_ld,
        d_prime_threshold=args.d_prime_threshold,
        haplotype_type=args.haplotype_type,
        aft=args.aft, 
        md=args.md,
        noheader=args.noheader,
        verbose=args.verbose
    )
    
    if not args.no_dedup:
        if args.verbose:
            print("\nDeduplicating blocks...")
        all_blocks = deduplicate_blocks(all_blocks, verbose=args.verbose)
    elif args.verbose:
        print("\nSkipping deduplication (--no-dedup flag set)")

    # Generate CSV statistics and analysis files
    method_label = f"{args.method}"
    if args.method == 'ld_haploblock':
        method_label = f"ld_haploblock_{args.haplotype_type}"

    write_csv_outputs(
        all_blocks=all_blocks,
        hap_files=args.input,
        map_file=args.map,
        output_dir=args.output,
        method_name=method_label,
        noheader=args.noheader,
        verbose=args.verbose
    )

    write_snp_selection_file(
        all_blocks=all_blocks,
        map_file=args.map,
        output_dir=args.output,
        verbose=args.verbose
    )

    # Write output files
    print("\nWriting output files:")
    total_blocks = write_output_files(all_blocks, args.output, verbose=args.verbose)

    if args.generate_genotypes is not None:
        # Determine output directory for genotypes
        # Get parent directory of args.output
        parent_dir = os.path.dirname(args.output) if os.path.dirname(args.output) else '.'
        # Auto-ensure different folder by appending suffix if same
        if args.generate_genotypes == os.path.basename(args.output):
            geno_output_dir = os.path.join(parent_dir, args.generate_genotypes + "_geno")
        else:
            geno_output_dir = os.path.join(parent_dir, args.generate_genotypes)
        
        if args.verbose:
            print(f"\nGenotype output directory: {geno_output_dir}")
        
        generate_haplotype_genotypes(
            all_blocks=all_blocks,
            hap_files=args.input,
            output_dir=geno_output_dir,
            missing_code=args.missing,
            missing_value=args.missing_value,
            noheader=args.noheader,
            verbose=args.verbose
        )

    # Calculate genome coverage statistics
    coverage_stats = calculate_genome_coverage(
        all_blocks=all_blocks,
        map_file=args.map,
        verbose=args.verbose
    )

    # Save coverage stats to CSV
    stats_dir = os.path.join(args.output, 'stats')
    os.makedirs(stats_dir, exist_ok=True)

    coverage_df = pd.DataFrame([{
        'method': method_label,
        'total_blocks': coverage_stats['total_blocks'],
        'genome_length_mb': coverage_stats['total_genome_length_bp'] / 1e6,
        'covered_length_kb': coverage_stats['covered_length_bp'] / 1e3,
        'coverage_pct': coverage_stats['coverage_pct'],
        'mean_spacing_kb': coverage_stats.get('mean_block_spacing_kb', 0)
    }])

    coverage_file = os.path.join(stats_dir, 'genome_coverage.csv')
    coverage_df.to_csv(coverage_file, index=False)

    if args.verbose:
        print(f"  ✓ Coverage stats saved: {coverage_file}")
    
    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"Chromosomes processed:  {len(all_blocks)}")
    print(f"Total blocks defined:   {total_blocks}")
    print(f"Output directory:       {args.output}/")
    print("="*70)
    
    print("\n✓ Block definition completed!")

if __name__ == '__main__':
    main()