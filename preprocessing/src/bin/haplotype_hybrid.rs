// src/bin/haplotype_hybrid.rs
use clap::{Parser, ValueEnum};
use std::path::Path;
use anyhow::Result;
use microhaplotype_preprocessing::modules::*;

#[derive(Parser, Debug)]
#[command(name = "haplotype_hybrid")]
#[command(about = "Haplotype block definition using hybrid approach for PHASED data")]
struct Args {
    #[arg(short = 'i', long = "input", required = true, num_args = 1..)]
    input: Vec<String>,

    #[arg(short = 'm', long = "map", required = true)]
    map: String,

    #[arg(long = "method", default_value = "ld-haploblock")]
    method: MethodArg,

    #[arg(long = "window-bp", default_value_t = 125)]
    window_bp: i32,

    #[arg(long = "haplotype-type", default_value = "pure")]
    haplotype_type: HaplotypeTypeArg,

    #[arg(short = 'w', long = "window", default_value_t = 4)]
    window: usize,

    #[arg(long = "min-snps", default_value_t = 2)]
    min_snps: usize,

    #[arg(long = "max-snps", default_value_t = usize::MAX)]
    max_snps: usize,

    #[arg(long = "d-prime-threshold", default_value_t = 0.45)]
    d_prime_threshold: f64,

    #[arg(long = "aft", default_value_t = 0.08)]
    aft: f64,

    #[arg(long = "md", default_value_t = 0.10)]
    md: f64,

    #[arg(long = "no-dedup")]
    no_dedup: bool,

    #[arg(long = "ld-prune")]
    ld_prune: bool,

    #[arg(long = "ld-prune-r2", default_value_t = 0.8)]
    ld_prune_r2: f64,

    #[arg(long = "ld-prune-window-bp", default_value_t = 500_000)]
    ld_prune_window_bp: i64,

    #[arg(long = "no-rare-priority")]
    no_rare_priority: bool,

    #[arg(long = "min-pic")]
    min_pic: Option<f64>,

    #[arg(long = "top-k")]
    top_k: Option<usize>,

    #[arg(long = "min-effective-alleles")]
    min_effective_alleles: Option<f64>,

    #[arg(long = "max-rare-alleles-prop")]
    max_rare_alleles_prop: Option<f64>,

    #[arg(long = "min-ld")]
    min_ld: Option<f64>,

    #[arg(long = "noheader", default_value_t = true)]
    noheader: bool,

    #[arg(short = 'o', long = "output", default_value = "hap_info_microhap")]
    output: String,

    #[arg(long = "generate-genotypes")]
    generate_genotypes: Option<String>,

    #[arg(long = "missing")]
    missing: Option<String>,

    #[arg(long = "missing-value", default_value = "-9999")]
    missing_value: String,

    #[arg(short = 'v', long = "verbose")]
    verbose: bool,
}

#[derive(Debug, Clone, ValueEnum)]
enum MethodArg {
    LdHaploblock,
    SnpCountSimple,
}

#[derive(Debug, Clone, ValueEnum)]
enum HaplotypeTypeArg {
    Pure,
    Micro,
}

fn main() -> Result<()> {
    let args = Args::parse();

    // Print configuration
    println!("\n{}", "=".repeat(70));
    println!("MICROHAPLOTYPE BLOCK DEFINITION - HYBRID (PHASED DATA)");
    println!("{}", "=".repeat(70));
    println!("Haplotype files:   {} files", args.input.len());
    println!("Map file:          {}", args.map);

    match args.method {
        MethodArg::LdHaploblock => {
            println!("Method:            LD-based haploblocks (Jonas 2017)");
            println!("D' threshold:      {}", args.d_prime_threshold);
            match args.haplotype_type {
                HaplotypeTypeArg::Pure => {
                    println!("Haplotype size:    {} SNPs (best selection)", args.window);
                }
                HaplotypeTypeArg::Micro => {
                    println!("Physical window:   {} bp", args.window_bp);
                }
            }
            println!("Criterion-B AFT:   {}", args.aft);
            println!("Criterion-B MD:    {}", args.md);
            if args.ld_prune {
                println!("LD pruning:        r²>{}, window={}bp", args.ld_prune_r2, args.ld_prune_window_bp);
            }
        }
        MethodArg::SnpCountSimple => {
            println!("Method:            snp_count_simple");
            println!("Window size:       {} SNPs", args.window);
            println!("SNP range:         {} - {}", args.min_snps, args.max_snps);
            if let Some(ld) = args.min_ld {
                println!("Min LD (r²):       {}", ld);
            }
        }
    }

    println!("Output directory:  {}", args.output);
    if args.min_pic.is_some() || args.top_k.is_some() {
        if let Some(p) = args.min_pic { println!("Min PIC filter:    {}", p); }
        if let Some(k) = args.top_k   { println!("Top-K blocks:      {}", k); }
    }
    println!("{}\n", "=".repeat(70));

    if let Some(v) = args.min_effective_alleles { println!("Min eff. alleles:  {}", v); }
    if let Some(v) = args.max_rare_alleles_prop { println!("Max rare prop:     {}", v); }
    println!("{}\n", "=".repeat(70));

    // Convert CLI args to config
    let config = BlockDefinitionConfig {
        method: match args.method {
            MethodArg::LdHaploblock => Method::LdHaploblock,
            MethodArg::SnpCountSimple => Method::SnpCountSimple,
        },
        haplotype_type: match args.haplotype_type {
            HaplotypeTypeArg::Pure => HaplotypeType::Pure,
            HaplotypeTypeArg::Micro => HaplotypeType::Micro,
        },
        window_bp: args.window_bp,
        window_snps: args.window,
        min_snps: args.min_snps,
        max_snps: args.max_snps,
        d_prime_threshold: args.d_prime_threshold,
        aft: args.aft,
        md: args.md,
        min_ld: args.min_ld,
        noheader: args.noheader,
        verbose: args.verbose,
    };

    // Read map file
    let snp_map = read_map_file(&args.map)?;

    // Expand glob patterns and sort by chromosome
    let hap_files = expand_and_sort_files(&args.input)?;

    // Define blocks
    let all_blocks = define_microhaplotype_blocks(&hap_files, &snp_map, &config)?;

    // Deduplication
    let all_blocks = if !args.no_dedup {
        if args.verbose {
            println!("\nDeduplicating blocks...");
        }
        deduplicate_blocks(all_blocks, args.verbose)
    } else {
        if args.verbose {
            println!("\nSkipping deduplication (--no-dedup flag set)");
        }
        all_blocks
    };

    let all_blocks = if args.ld_prune {
        if args.verbose {
            println!("\nLD pruning blocks (r²>{}, window={}bp, rare_priority={})...",
                args.ld_prune_r2, args.ld_prune_window_bp, !args.no_rare_priority);
        }
        let prune_config = LdPruneConfig {
            r2_threshold: args.ld_prune_r2,
            window_bp: args.ld_prune_window_bp,
            prioritize_rare: !args.no_rare_priority,
            rare_freq_threshold: 0.05,
        };
        ld_prune_blocks(all_blocks, &hap_files, &prune_config, args.noheader, args.verbose)
    } else {
        all_blocks
    };

    let all_blocks = if args.min_effective_alleles.is_some() || args.max_rare_alleles_prop.is_some() {
        if args.verbose {
            println!("\nFiltering blocks by allele quality...");
        }
        filter_blocks_by_allele_quality(
            all_blocks,
            &hap_files,
            args.min_effective_alleles,
            args.max_rare_alleles_prop,
            args.noheader,
            args.verbose,
        )
    } else {
        all_blocks
    };

    // Ranking 
    let all_blocks = if args.min_pic.is_some() || args.top_k.is_some() {
        if args.verbose {
            println!("\nRanking/filtering blocks by PIC...");
        }
        rank_and_filter_blocks(
            all_blocks,
            &hap_files,
            args.min_pic,
            args.top_k,
            args.noheader,
            args.verbose,
        )
    } else {
        all_blocks
    };

    // Generate CSV statistics
    let method_label = match config.method {
        Method::LdHaploblock => match config.haplotype_type {
            HaplotypeType::Pure => "ld_haploblock_pure",
            HaplotypeType::Micro => "ld_haploblock_micro",
        },
        Method::SnpCountSimple => "snp_count_simple",
    };

    write_csv_outputs(
        &all_blocks,
        &hap_files,
        &snp_map,
        &args.output,
        method_label,
        args.noheader,
        args.verbose,
    )?;

    write_snp_selection_file(&all_blocks, &snp_map, &args.output, args.verbose)?;

    // Write output files
    println!("\nWriting output files:");
    let total_blocks = write_output_files(&all_blocks, &args.output, args.verbose)?;

    // Generate genotypes if requested
    if let Some(ref geno_prefix) = args.generate_genotypes {
        let parent_dir = Path::new(&args.output)
            .parent()
            .unwrap_or(Path::new("."));
        
        let geno_output_dir = if geno_prefix == Path::new(&args.output).file_name().unwrap().to_str().unwrap() {
            parent_dir.join(format!("{}_geno", geno_prefix))
        } else {
            parent_dir.join(geno_prefix)
        };

        if args.verbose {
            println!("\nGenotype output directory: {}", geno_output_dir.display());
        }

        generate_haplotype_genotypes(
            &all_blocks,
            &hap_files,
            &geno_output_dir,
            args.missing.as_deref(),
            &args.missing_value,
            args.noheader,
            args.verbose,
        )?;
    }

    // Calculate genome coverage
    let coverage_stats = calculate_genome_coverage(&all_blocks, &snp_map, args.verbose)?;

    // Save coverage stats
    save_coverage_stats(&coverage_stats, &args.output, method_label, args.verbose)?;

    // Summary
    println!("\n{}", "=".repeat(70));
    println!("SUMMARY");
    println!("{}", "=".repeat(70));
    println!("Chromosomes processed:  {}", all_blocks.len());
    println!("Total blocks defined:   {}", total_blocks);
    println!("Output directory:       {}/", args.output);
    println!("{}", "=".repeat(70));
    println!("\n✓ Block definition completed!");

    Ok(())
}