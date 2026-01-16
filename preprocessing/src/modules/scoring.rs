// src/modules/scoring.rs
use std::collections::HashMap;
use ndarray::{Array2, Axis};

use super::types::{SnpInfo, Haploblock, SelectedHaplotype, Microhaplotype, AlleleMetrics};

pub fn criterion_b_score(block_haps: &Array2<f64>, aft: f64, md: f64) -> f64 {
    let n_snps = block_haps.ncols();
    let hs = 2_f64.powi(n_snps as i32);

    let mut unique_haps = HashMap::new();
    for row in block_haps.rows() {
        let has_nan = row.iter().any(|v| v.is_nan());
        if !has_nan {
            let hap_tuple: Vec<i32> = row.iter().map(|v| *v as i32).collect();
            *unique_haps.entry(hap_tuple).or_insert(0) += 1;
        }
    }

    if unique_haps.is_empty() {
        return f64::INFINITY;
    }

    let total = unique_haps.values().sum::<i32>() as f64;
    let freqs: Vec<f64> = unique_haps.values().map(|&c| c as f64 / total).collect();

    let n_predictable = freqs.iter().filter(|&&f| f >= aft).count();
    
    if n_predictable == 0 {
        return f64::INFINITY;
    }

    let pred_freqs: Vec<f64> = freqs.iter().filter(|&&f| f >= aft).copied().collect();
    let ideal_freq = 1.0 / hs;
    let sd_term: f64 = pred_freqs.iter().map(|f| (f - ideal_freq).powi(2)).sum();

    let w = if n_predictable == 1 {
        0.0
    } else {
        (md * n_predictable as f64) / (hs * (n_predictable as f64 - 1.0))
    };

    sd_term - (w * n_predictable as f64)
}

pub fn select_best_haplotype_from_haploblock(
    haploblock: &Haploblock,
    hap_matrix: &Array2<f64>,
    haplotype_size: usize,
    aft: f64,
    md: f64,
) -> SelectedHaplotype {
    let snp_indices = &haploblock.snp_indices;
    let n_snps = snp_indices.len();

    if n_snps < haplotype_size {
        let start = snp_indices[0];
        let end = *snp_indices.last().unwrap();
        let needed = haplotype_size - n_snps;
        let left = needed / 2;
        let right = needed - left;

        let new_start = start.saturating_sub(left);
        let new_end = (end + right).min(hap_matrix.ncols() - 1);

        let mut extended: Vec<usize> = (new_start..=new_end).collect();
        if extended.len() > haplotype_size {
            extended.truncate(haplotype_size);
        } else {
            while extended.len() < haplotype_size && *extended.last().unwrap() < hap_matrix.ncols() - 1 {
                extended.push(extended.last().unwrap() + 1);
            }
        }

        let selected_indices: Vec<usize> = extended.iter().take(haplotype_size).copied().collect();
        let block_haps = hap_matrix.select(Axis(1), &selected_indices);
        let score = criterion_b_score(&block_haps, aft, md);

        SelectedHaplotype {
            snp_indices: selected_indices,
            n_snps: haplotype_size,
            criterion_b_score: score,
            original_block_size: n_snps,
            selection_type: "extended".to_string(),
        }
    } else if n_snps == haplotype_size {
        let block_haps = hap_matrix.select(Axis(1), snp_indices);
        let score = criterion_b_score(&block_haps, aft, md);

        SelectedHaplotype {
            snp_indices: snp_indices.clone(),
            n_snps,
            criterion_b_score: score,
            original_block_size: n_snps,
            selection_type: "exact".to_string(),
        }
    } else {
        let mut best_score = f64::INFINITY;
        let mut best_indices = Vec::new();

        for start_pos in 0..=(n_snps - haplotype_size) {
            let window_indices: Vec<usize> = snp_indices[start_pos..start_pos + haplotype_size].to_vec();
            let block_haps = hap_matrix.select(Axis(1), &window_indices);
            let score = criterion_b_score(&block_haps, aft, md);

            if score < best_score {
                best_score = score;
                best_indices = window_indices;
            }
        }

        SelectedHaplotype {
            snp_indices: best_indices,
            n_snps: haplotype_size,
            criterion_b_score: best_score,
            original_block_size: n_snps,
            selection_type: "selected".to_string(),
        }
    }
}

pub fn split_ld_block_by_physical_window(
    haploblock: &Haploblock,
    hap_matrix: &Array2<f64>,
    chr_snps: &[SnpInfo],
    window_bp: i32,
    min_snps: usize,
    aft: f64,
    md: f64,
) -> Vec<Microhaplotype> {
    let snp_indices = &haploblock.snp_indices;
    let positions: Vec<i64> = snp_indices.iter().map(|&i| chr_snps[i].position).collect();
    let start_pos = positions[0];
    let end_pos = *positions.last().unwrap();
    let block_span = end_pos - start_pos;

    if block_span <= window_bp as i64 {
        if snp_indices.len() >= min_snps {
            let block_haps = hap_matrix.select(Axis(1), snp_indices);
            let score = criterion_b_score(&block_haps, aft, md);

            return vec![Microhaplotype {
                snp_indices: snp_indices.clone(),
                n_snps: snp_indices.len(),
                criterion_b_score: score,
                physical_span: block_span,
                split_type: "single".to_string(),
            }];
        } else {
            return vec![];
        }
    }

    let mut microhaplotypes = Vec::new();
    let mut current_start_idx = 0;

    while current_start_idx < snp_indices.len() {
        let window_start_pos = positions[current_start_idx];
        let window_end_pos = window_start_pos + window_bp as i64;

        let mut window_snp_indices = Vec::new();
        for i in current_start_idx..snp_indices.len() {
            if positions[i] < window_end_pos {
                window_snp_indices.push(snp_indices[i]);
            } else {
                break;
            }
        }

        if window_snp_indices.len() >= min_snps {
            let block_haps = hap_matrix.select(Axis(1), &window_snp_indices);
            let score = criterion_b_score(&block_haps, aft, md);
            let actual_span = positions[current_start_idx + window_snp_indices.len() - 1] - window_start_pos;

            microhaplotypes.push(Microhaplotype {
                snp_indices: window_snp_indices.clone(),
                n_snps: window_snp_indices.len(),
                criterion_b_score: score,
                physical_span: actual_span,
                split_type: "split".to_string(),
            });
        }

        current_start_idx += window_snp_indices.len();
        if window_snp_indices.is_empty() {
            current_start_idx += 1;
        }
    }

    microhaplotypes
}

pub fn calculate_allele_metrics(hap_matrix: &Array2<f64>, block_start_idx: usize, block_end_idx: usize) -> Option<AlleleMetrics> {
    let snp_indices: Vec<usize> = (block_start_idx..=block_end_idx).collect();
    let block_haps = hap_matrix.select(Axis(1), &snp_indices);
    
    // Get unique haplotype patterns
    let mut unique_haps = Vec::new();
    for row in block_haps.rows() {
        let has_nan = row.iter().any(|v| v.is_nan());
        if !has_nan {
            let hap_tuple: Vec<i32> = row.iter().map(|v| *v as i32).collect();
            unique_haps.push(hap_tuple);
        }
    }
    
    if unique_haps.is_empty() {
        return None;
    }
    
    // Count frequencies
    let mut counts = HashMap::new();
    for hap in &unique_haps {
        *counts.entry(hap.clone()).or_insert(0) += 1;
    }
    
    let total = unique_haps.len() as f64;
    let freqs: Vec<f64> = counts.values().map(|&c| c as f64 / total).collect();
    
    // Basic metrics
    let n_alleles = freqs.len();
    let n_predictable = freqs.iter().filter(|&&f| f >= 0.08).count();
    let min_freq = freqs.iter().copied().fold(f64::INFINITY, f64::min);
    let max_freq = freqs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mean_freq = freqs.iter().sum::<f64>() / n_alleles as f64;
    
    // Shannon entropy
    let entropy = -freqs.iter()
        .map(|&f| if f > 0.0 { f * f.log2() } else { 0.0 })
        .sum::<f64>();
    
    // Expected heterozygosity (He)
    let he = 1.0 - freqs.iter().map(|f| f.powi(2)).sum::<f64>();
    
    // Observed heterozygosity (Ho)
    let n_individuals = block_haps.nrows() / 2;
    let mut heterozygous_count = 0;
    
    for i in (0..block_haps.nrows()).step_by(2) {
        if i + 1 < block_haps.nrows() {
            let hap1: Vec<i32> = block_haps.row(i)
                .iter()
                .filter(|v| !v.is_nan())
                .map(|v| *v as i32)
                .collect();
            let hap2: Vec<i32> = block_haps.row(i + 1)
                .iter()
                .filter(|v| !v.is_nan())
                .map(|v| *v as i32)
                .collect();
            
            if !hap1.is_empty() && !hap2.is_empty() && hap1 != hap2 {
                heterozygous_count += 1;
            }
        }
    }
    
    let ho = if n_individuals > 0 {
        heterozygous_count as f64 / n_individuals as f64
    } else {
        0.0
    };
    
    // PIC (Polymorphic Information Content)
    let sum_squared = freqs.iter().map(|f| f.powi(2)).sum::<f64>();
    let mut cross_term = 0.0;
    for i in 0..freqs.len() {
        for j in (i + 1)..freqs.len() {
            cross_term += 2.0 * freqs[i].powi(2) * freqs[j].powi(2);
        }
    }
    let pic = 1.0 - sum_squared - cross_term;
    
    // Effective number of alleles
    let ne = 1.0 / freqs.iter().map(|f| f.powi(2)).sum::<f64>();
    
    // Rare alleles proportion (freq < 0.05)
    let rare_prop = freqs.iter().filter(|&&f| f < 0.05).count() as f64 / n_alleles as f64;
    
    Some(AlleleMetrics {
        n_alleles,
        n_predictable,
        min_freq,
        max_freq,
        mean_freq,
        freq_entropy: entropy,
        expected_heterozygosity: he,
        observed_heterozygosity: ho,
        pic,
        effective_alleles: ne,
        rare_alleles_prop: rare_prop,
    })
}