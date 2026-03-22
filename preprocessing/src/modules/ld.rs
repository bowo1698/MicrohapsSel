// src/modules/ld.rs
use ndarray::Array2;
use super::types::Haploblock;

pub fn calculate_d_prime(snp1: &[f64], snp2: &[f64]) -> f64 {
    let pairs: Vec<(f64, f64)> = snp1
        .iter()
        .zip(snp2.iter())
        .filter(|(a, b)| !a.is_nan() && !b.is_nan())
        .map(|(a, b)| (*a, *b))
        .collect();

    if pairs.len() < 10 {
        return 0.0;
    }

    let n = pairs.len() as f64;
    let p_a: f64 = pairs.iter().map(|(a, _)| a).sum::<f64>() / n;
    let p_b: f64 = pairs.iter().map(|(_, b)| b).sum::<f64>() / n;
    let p_ab: f64 = pairs.iter().map(|(a, b)| a * b).sum::<f64>() / n;

    let d = p_ab - (p_a * p_b);

    let d_max = if d >= 0.0 {
        f64::min(p_a * (1.0 - p_b), (1.0 - p_a) * p_b)
    } else {
        f64::min(p_a * p_b, (1.0 - p_a) * (1.0 - p_b))
    };

    if d_max == 0.0 {
        0.0
    } else {
        d.abs() / d_max
    }
}

pub fn build_ld_haploblocks(
    hap_matrix: &Array2<f64>,
    snp_positions: Option<&[i64]>,
    d_prime_threshold: f64,
    min_block_snps: usize,
    verbose: bool,
) -> Vec<Haploblock> {
    let n_snps = hap_matrix.ncols();

    if verbose {
        println!("    Building haploblocks with D' threshold {}", d_prime_threshold);
    }

    let mut d_primes = Vec::new();
    let mut prev_col = hap_matrix.column(0).to_vec();
    for i in 0..n_snps - 1 {
        let next_col = hap_matrix.column(i + 1).to_vec();
        let d_prime = calculate_d_prime(&prev_col, &next_col);
        d_primes.push(d_prime);
        prev_col = next_col;  // reuse
    }

    if verbose {
        let mean_d = d_primes.iter().sum::<f64>() / d_primes.len() as f64;
        let above_thresh = d_primes.iter().filter(|&&d| d >= d_prime_threshold).count();
        println!("    Mean consecutive D': {:.3}", mean_d);
        println!("    Pairs ≥ threshold: {}/{}", above_thresh, d_primes.len());
    }

    let mut haploblocks = Vec::new();
    let mut current_block = vec![0];

    for (i, &d_prime) in d_primes.iter().enumerate() {
        if d_prime >= d_prime_threshold {
            current_block.push(i + 1);
        } else {
            if current_block.len() >= min_block_snps {
                let block_start = current_block[0];
                let block_end = *current_block.last().unwrap();

                let physical_span = if let Some(pos) = snp_positions {
                    Some(pos[block_end] - pos[block_start])
                } else {
                    None
                };

                haploblocks.push(Haploblock {
                    start_idx: block_start,
                    end_idx: block_end,
                    snp_indices: current_block.clone(),
                    n_snps: current_block.len(),
                    physical_span,
                });
            }
            current_block = vec![i + 1];
        }
    }

    if current_block.len() >= min_block_snps {
        let block_start = current_block[0];
        let block_end = *current_block.last().unwrap();

        let physical_span = if let Some(pos) = snp_positions {
            Some(pos[block_end] - pos[block_start])
        } else {
            None
        };

        haploblocks.push(Haploblock {
            start_idx: block_start,
            end_idx: block_end,
            snp_indices: current_block.clone(),
            n_snps: current_block.len(),
            physical_span,
        });
    }

    if verbose {
        let sizes: Vec<usize> = haploblocks.iter().map(|b| b.n_snps).collect();
        println!("    Built {} haploblocks", haploblocks.len());
        if !sizes.is_empty() {
            let mean = sizes.iter().sum::<usize>() as f64 / sizes.len() as f64;
            let mut sorted = sizes.clone();
            sorted.sort();
            let median = sorted[sorted.len() / 2];
            println!("    Block size - mean: {:.1}, median: {}, range: [{}, {}]",
                     mean, median, sizes.iter().min().unwrap(), sizes.iter().max().unwrap());

            if let Some(block) = haploblocks.first() {
                if let Some(_) = block.physical_span {
                    let spans: Vec<i64> = haploblocks.iter()
                        .filter_map(|b| b.physical_span)
                        .collect();
                    let mean_span = spans.iter().sum::<i64>() as f64 / spans.len() as f64;
                    let mut sorted_spans = spans.clone();
                    sorted_spans.sort();
                    let median_span = sorted_spans[sorted_spans.len() / 2];
                    println!("    Physical span (bp) - mean: {:.0}, median: {:.0}, range: [{:.0}, {:.0}]",
                             mean_span, median_span, 
                             spans.iter().min().unwrap(), 
                             spans.iter().max().unwrap());
                }
            }
        }
    }

    haploblocks
}