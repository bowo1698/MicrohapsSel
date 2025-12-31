suppressPackageStartupMessages({
  library(AlphaSimR)
  library(tidyverse)
  library(vcfR)
  library(data.table)
})

# ============================================================================
# SETUP
# ============================================================================

output_dir <- "output/trout_gens"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

vcf_file <- "output/trout/Trout_phased.vcf.gz"
map_file <- "output/trout/Trout_CD_map.txt"

# Simulation parameters
n_iterations <- 5
n_generations <- 3
n_sires <- 40
n_dams_per_sire <- 10
n_offspring_per_family <- 20
n_families <- n_sires * n_dams_per_sire  # 400
n_offspring_per_gen <- n_families * n_offspring_per_family  # 8,000
h2_target <- 0.4
n_qtl_per_chr <- 100
top_pct_selection <- 0.4
organism <- "trout"

# ============================================================================
# LOAD & PROCESS VCF (FOUNDERS)
# ============================================================================

cat("Reading VCF file...\n")
vcf <- read.vcfR(vcf_file, verbose = FALSE)

cat("Parsing phased genotypes...\n")
gt <- extract.gt(vcf, element = "GT")
n_founders <- ncol(gt)
n_snps <- nrow(gt)

cat(sprintf("  Loaded: %d SNPs, %d founders\n", n_snps, n_founders))

# Parse haplotypes
hap_matrix <- matrix(0, nrow = 2 * n_founders, ncol = n_snps)
for(i in 1:n_founders) {
  for(j in 1:n_snps) {
    alleles <- strsplit(gt[j,i], "\\|")[[1]]
    hap_matrix[2*i-1, j] <- as.numeric(alleles[1])
    hap_matrix[2*i, j] <- as.numeric(alleles[2])
  }
}

# Filter monomorphic
cat("Filtering monomorphic SNPs...\n")
snp_var <- apply(hap_matrix, 2, var, na.rm = TRUE)
polymorphic <- snp_var > 0 & !is.na(snp_var)
cat(sprintf("  Monomorphic SNPs removed: %d\n", sum(!polymorphic)))

hap_matrix <- hap_matrix[, polymorphic]
vcf_fix_filtered <- vcf@fix[polymorphic, ]
cat(sprintf("  Remaining polymorphic SNPs: %d\n", ncol(hap_matrix)))

# Create genetic map by chromosome
cat("Creating genetic map by chromosome...\n")
chr_numbers <- as.numeric(vcf_fix_filtered[,"CHROM"])
chr_list <- sort(unique(chr_numbers))
genMap <- vector("list", length(chr_list))
haplotypes <- vector("list", length(chr_list))

for(i in seq_along(chr_list)) {
  chr <- chr_list[i]
  chr_idx <- which(chr_numbers == chr)
  positions <- as.numeric(vcf_fix_filtered[chr_idx, "POS"])
  
  if(any(duplicated(positions))) {
    cat(sprintf("  WARNING: Chr %d has %d duplicate positions, removing...\n", 
                chr, sum(duplicated(positions))))
    keep_idx <- !duplicated(positions)
    chr_idx <- chr_idx[keep_idx]
    positions <- positions[keep_idx]
  }
  
  # Convert physical position (bp) to genetic position (Morgans)
  genetic_pos <- (positions - min(positions)) / 1e8
  genMap[[i]] <- genetic_pos
  haplotypes[[i]] <- hap_matrix[, chr_idx, drop = FALSE]

  # Check
  if(length(genMap[[i]]) != ncol(haplotypes[[i]])) {
    stop(sprintf("Mismatch at chr %d: genMap=%d, hap=%d", 
                 chr, length(genMap[[i]]), ncol(haplotypes[[i]])))
  }
  
  cat(sprintf("  Chr %d: %d SNPs, genetic length %.3f M\n", 
              chr, length(genMap[[i]]), max(genetic_pos)))
}

cat(sprintf("  Total chromosomes: %d\n", length(genMap)))

# Final check
cat("Final validation...\n")
total_genmap <- sum(sapply(genMap, length))
total_hap <- sum(sapply(haplotypes, ncol))
cat(sprintf("  Total SNPs in genMap: %d\n", total_genmap))
cat(sprintf("  Total SNPs in haplotypes: %d\n", total_hap))

if(total_genmap != total_hap) {
  stop("FATAL: Total mismatch!")
}

# Create founder population
cat("Creating founder population...\n")
founderPop <- newMapPop(genMap = genMap, haplotypes = haplotypes, ploidy = 2L)
cat("  SUCCESS: Founder population created\n")

SP_base <- SimParam$new(founderPop)
founders <- newPop(founderPop, simParam = SP_base)

# Founders genotype & pedigree (NO phenotype)
founders_geno <- pullSegSiteGeno(founders, simParam = SP_base)
snp_names_base <- sprintf("SNP_%d", 1:ncol(founders_geno))
colnames(founders_geno) <- snp_names_base

founders_ped <- tibble(
  individual_id = paste0("founder_", 1:nInd(founders)),
  generation = 0,
  population_type = "founder",
  mother_id = NA,
  father_id = NA,
  iteration = NA,
  family_id = NA
)

cat("Creating SNP map file...\n")

# Build chromosome and position vectors
snp_chr <- integer(ncol(founders_geno))
snp_pos <- numeric(ncol(founders_geno))
snp_idx <- 1

n_chr <- length(genMap)

for(chr_i in 1:n_chr) {
  n_snp_chr <- ncol(haplotypes[[chr_i]])
  snp_chr[snp_idx:(snp_idx + n_snp_chr - 1)] <- chr_list[chr_i]
  snp_pos[snp_idx:(snp_idx + n_snp_chr - 1)] <- genMap[[chr_i]] * 1e8  # Morgan to bp
  snp_idx <- snp_idx + n_snp_chr
}

# Create map dataframe
map_df <- tibble(
  SNPID = sprintf("AX-%08d", 1:ncol(founders_geno)),
  Chr = snp_chr,
  Position = round(snp_pos)
)

write_tsv(map_df, file.path(output_dir, sprintf("%s_snp_map.txt", organism)))
cat(sprintf("  SNP map saved: %d SNPs\n", nrow(map_df)))

cat(sprintf("Founders: %d individuals, %d SNPs\n", nInd(founders), ncol(founders_geno)))

# ============================================================================
# SAMPLE QTL POSITIONS (ONCE, BEFORE ITERATION LOOP)
# ============================================================================
total_qtl <- n_chr * n_qtl_per_chr

set.seed(456)
qtl_idx <- integer(total_qtl)
idx <- 1
chr_cumsum <- c(0, cumsum(sapply(haplotypes, ncol)))

for(chr_i in 1:n_chr) {
  n_seg_chr <- ncol(haplotypes[[chr_i]])
  qtl_chr <- sort(sample(1:n_seg_chr, n_qtl_per_chr, replace = FALSE))
  qtl_idx[idx:(idx + n_qtl_per_chr - 1)] <- chr_cumsum[chr_i] + qtl_chr
  idx <- idx + n_qtl_per_chr
}

cat(sprintf("QTL positions sampled: %d QTLs\n", total_qtl))

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

select_parents <- function(phenotypes, pop_object, n_sires, n_dams, top_pct = 0.4) {
  # Rank by TBV
  ranked <- phenotypes %>%
    arrange(desc(tbv)) %>%
    mutate(
      orig_index = row_number(),
      rank_pct = orig_index / n()
    )
  
  # Top X% pool
  top_pool_idx <- ranked %>%
    filter(rank_pct <= top_pct) %>%
    pull(orig_index)
  
  # Random sample from top pool
  sire_idx <- sample(top_pool_idx, n_sires)
  dam_idx <- sample(setdiff(top_pool_idx, sire_idx), n_dams)
  
  return(list(sire_idx = sire_idx, dam_idx = dam_idx))
}

create_balanced_crosses <- function(sire_idx, dam_idx, dams_per_sire) {
  n_sires <- length(sire_idx)
  n_crosses <- n_sires * dams_per_sire
  
  crossPlan <- matrix(0, nrow = n_crosses, ncol = 2)
  cross_idx <- 1
  
  for(i in 1:n_sires) {
    dam_start <- (i - 1) * dams_per_sire + 1
    dam_end <- i * dams_per_sire
    
    for(j in dam_start:dam_end) {
      crossPlan[cross_idx, 1] <- dam_idx[j]  # female
      crossPlan[cross_idx, 2] <- sire_idx[i]  # male
      cross_idx <- cross_idx + 1
    }
  }
  
  return(crossPlan)
}

# ============================================================================
# OUTER LOOP: 5 ITERATIONS
# ============================================================================

for(iter in 1:n_iterations) {
  
  cat(sprintf("\n========================================\n"))
  cat(sprintf("ITERATION %d\n", iter))
  cat(sprintf("========================================\n"))
  
  # Set unique seed per iteration
  set.seed(1000 + iter * 100)
  
  # ------------------------------
  # SAMPLE QTL EFFECTS (PER ITERATION)
  # ------------------------------
  
  # Normal QTL: all from N(0,1)
  beta_normal <- rnorm(total_qtl, mean = 0, sd = 1)
  beta_normal <- beta_normal / sd(beta_normal)
  
  # Mixture QTL: 80-15-5
  n_small <- round(total_qtl * 0.80)
  n_medium <- round(total_qtl * 0.15)
  n_large <- total_qtl - n_small - n_medium
  
  beta_mixture <- c(
    rnorm(n_small, 0, 0.1),
    rnorm(n_medium, 0, 0.5),
    rnorm(n_large, 0, 1.0)
  )
  beta_mixture <- sample(beta_mixture)
  beta_mixture <- beta_mixture / sd(beta_mixture)
  
  # Save QTL effects
  qtl_info <- tibble(
    qtl_id = 1:total_qtl,
    snp_index = qtl_idx,
    beta_normal = beta_normal,
    beta_mixture = beta_mixture,
    effect_class = c(rep("small", n_small), 
                     rep("medium", n_medium), 
                     rep("large", n_large))
  )
  
  iter_dir <- file.path(output_dir, sprintf("iteration_%d", iter))
  dir.create(iter_dir, showWarnings = FALSE, recursive = TRUE)
  write_csv(qtl_info, file.path(iter_dir, sprintf("%s_qtl_effects.csv", organism)))
  cat(sprintf("  ✓ QTL effects sampled and saved: %d QTLs\n", total_qtl))
  cat(sprintf("    - Normal: sd(beta) = %.3f\n", sd(beta_normal)))
  cat(sprintf("    - Mixture: sd(beta) = %.3f (80%%/15%%/5%% = %d/%d/%d)\n", 
              sd(beta_mixture), n_small, n_medium, n_large))
  
  # ------------------------------
  # INITIALIZE ACCUMULATORS
  # ------------------------------
  
  genotypes_list <- list(founders_geno)
  all_phenotypes_normal <- tibble()
  all_phenotypes_mixture <- tibble()
  all_pedigree <- founders_ped
  
  pop_objects <- list()
  pop_objects[["founders"]] <- founders
  
  # ------------------------------
  # INNER LOOP: 3 GENERATIONS
  # ------------------------------
  
  for(gen in 1:n_generations) {
    
    cat(sprintf("\n--- Generation %d ---\n", gen))

    tryCatch({
        # PARENT SELECTION
        if(gen == 1) {
        # Gen 1: random from founders
        parent_pop <- founders
        sire_idx <- sample(1:nInd(founders), n_sires)
        #sire_idx <- sample(1:nInd(founders), n_sires, replace = FALSE)
        #remaining_founders <- setdiff(1:nInd(founders), sire_idx)
        #dam_idx <- sample(remaining_founders, n_families, replace = TRUE)
        dam_idx <- sample(setdiff(1:nInd(founders), sire_idx), n_families)
        
        } else {
        # Gen 2+: selection from previous generation
        prev_gen_pheno <- all_phenotypes_normal %>%
            filter(generation == (gen - 1))
        
        parent_pop <- pop_objects[[sprintf("gen%d", gen - 1)]]
        selected <- select_parents(prev_gen_pheno, parent_pop, n_sires, n_families, top_pct_selection)
        sire_idx <- selected$sire_idx
        dam_idx <- selected$dam_idx
        cat(sprintf("  ✓ Selected parents: %d sires, %d dams from top %.0f%% (n=%d)\n",
                    n_sires, n_families, top_pct_selection*100, nrow(prev_gen_pheno)))
        }
        
        # CREATE CROSS PLAN
        crossPlan <- create_balanced_crosses(sire_idx, dam_idx, n_dams_per_sire)
        cat(sprintf("  ✓ Cross plan created: %d crosses (%d sires × %d dams/sire)\n",
                    nrow(crossPlan), n_sires, n_dams_per_sire))
        
        # GENERATE OFFSPRING
        offspring <- makeCross2(
        females = parent_pop,
        males = parent_pop,
        crossPlan = crossPlan,
        nProgeny = n_offspring_per_family,
        simParam = SP_base
        )
        
        pop_objects[[sprintf("gen%d", gen)]] <- offspring
        cat(sprintf("  ✓ Offspring generated: %d individuals from %d families\n",
                    nInd(offspring), n_families))
        
        # EXTRACT GENOTYPES
        offspring_geno <- pullSegSiteGeno(offspring, simParam = SP_base)
        cat(sprintf("  ✓ Genotypes extracted: %d ind × %d SNPs\n",
                    nrow(offspring_geno), ncol(offspring_geno)))
        
        # CALCULATE PHENOTYPES (USING SAME BETA FROM THIS ITERATION)
        X_qtl <- offspring_geno[, qtl_idx]
        
        # Normal scenario
        tbv_normal <- as.vector(X_qtl %*% beta_normal)
        var_g_normal <- var(tbv_normal)
        var_e_normal <- var_g_normal * (1/h2_target - 1)
        e_normal <- rnorm(nInd(offspring), 0, sqrt(var_e_normal))
        pheno_normal <- tbv_normal + e_normal
        
        # Mixture scenario
        tbv_mixture <- as.vector(X_qtl %*% beta_mixture)
        var_g_mixture <- var(tbv_mixture)
        var_e_mixture <- var_g_mixture * (1/h2_target - 1)
        e_mixture <- rnorm(nInd(offspring), 0, sqrt(var_e_mixture))
        pheno_mixture <- tbv_mixture + e_mixture

        cat(sprintf("  ✓ Phenotypes calculated:\n"))
        cat(sprintf("    - Normal:  h² = %.3f, mean(TBV) = %.3f, sd(TBV) = %.3f\n",
                    var_g_normal / var(pheno_normal), mean(tbv_normal), sd(tbv_normal)))
        cat(sprintf("    - Mixture: h² = %.3f, mean(TBV) = %.3f, sd(TBV) = %.3f\n",
                    var_g_mixture / var(pheno_mixture), mean(tbv_mixture), sd(tbv_mixture)))
        
        # CREATE IDs
        family_ids <- rep(1:n_families, each = n_offspring_per_family)
        ind_in_family <- rep(1:n_offspring_per_family, times = n_families)
        offspring_ids <- sprintf("iter%d_gen%d_fam%03d_ind%02d", 
                                iter, gen, family_ids, ind_in_family)
        
        # CREATE DATAFRAMES
        offspring_pheno_normal <- tibble(
            individual_id = offspring_ids,
            tbv = tbv_normal,
            phenotype = pheno_normal,
            generation = gen,
            population_type = "reference",
            iteration = iter,
            family_id = family_ids,
            h2_realized = var_g_normal / var(pheno_normal)
        )
        
        offspring_pheno_mixture <- tibble(
            individual_id = offspring_ids,
            tbv = tbv_mixture,
            phenotype = pheno_mixture,
            generation = gen,
            population_type = "reference",
            iteration = iter,
            family_id = family_ids,
            h2_realized = var_g_mixture / var(pheno_mixture)
        )
        
        # Get parent IDs based on generation using cross plan
        mother_ids <- rep(NA_character_, nInd(offspring))
        father_ids <- rep(NA_character_, nInd(offspring))

        for(i in 1:nrow(crossPlan)) {
          # Each cross produces n_offspring_per_family offspring
          offspring_idx <- ((i-1) * n_offspring_per_family + 1):(i * n_offspring_per_family)
          
          if(gen == 1) {
            mother_ids[offspring_idx] <- founders_ped$individual_id[crossPlan[i, 1]]
            father_ids[offspring_idx] <- founders_ped$individual_id[crossPlan[i, 2]]
          } else {
            prev_gen_ids <- all_pedigree %>%
              filter(generation == (gen - 1), iteration == iter) %>%
              pull(individual_id)
            
            mother_ids[offspring_idx] <- prev_gen_ids[crossPlan[i, 1]]
            father_ids[offspring_idx] <- prev_gen_ids[crossPlan[i, 2]]
          }
        }

        offspring_ped <- tibble(
          individual_id = offspring_ids,
          generation = gen,
          population_type = "reference",
          mother_id = mother_ids,
          father_id = father_ids,
          iteration = iter,
          family_id = family_ids
        )

        # QC
        cat("  QC: Checking pedigree integrity...\n")
        # Build reference pedigree (all previous generations)
        if(gen == 1) {
        ref_pedigree <- founders_ped
        } else {
        ref_pedigree <- bind_rows(
            founders_ped,
            all_pedigree %>% filter(generation < gen)
        )
        }
        
        # Check mother IDs
        ped_check_mother <- offspring_ped %>%
        filter(!is.na(mother_id)) %>%
        anti_join(ref_pedigree, by = c("mother_id" = "individual_id"))
        
        # Check father IDs
        ped_check_father <- offspring_ped %>%
        filter(!is.na(father_id)) %>%
        anti_join(ref_pedigree, by = c("father_id" = "individual_id"))
        
        if(nrow(ped_check_mother) > 0) {
        stop(sprintf("Gen %d: %d offspring with missing mother IDs!", 
                    gen, nrow(ped_check_mother)))
        }
        
        if(nrow(ped_check_father) > 0) {
        stop(sprintf("Gen %d: %d offspring with missing father IDs!", 
                    gen, nrow(ped_check_father)))
        }
        
        cat(sprintf("  QC: ✓ All %d offspring have valid parent IDs\n", nrow(offspring_ped)))
        cat(sprintf("  ✓ Pedigree dataframe created: %d rows\n", nrow(offspring_ped)))
        
        # ACCUMULATE DATA
        colnames(offspring_geno) <- NULL
        colnames(offspring_geno) <- colnames(genotypes_list[[1]])
        genotypes_list[[length(genotypes_list) + 1]] <- offspring_geno
        all_phenotypes_normal <- rbind(all_phenotypes_normal, offspring_pheno_normal)
        all_phenotypes_mixture <- rbind(all_phenotypes_mixture, offspring_pheno_mixture)
        all_pedigree <- rbind(all_pedigree, offspring_ped)

        cat(sprintf("  ✓ Data accumulated:\n"))
        total_geno_rows <- sum(sapply(genotypes_list, nrow))
        cat(sprintf("    - Genotypes: %d → %d individuals\n",
                    total_geno_rows - nrow(offspring_geno), total_geno_rows))
        cat(sprintf("    - Phenotypes: %d → %d individuals\n",
                    nrow(all_phenotypes_normal) - nrow(offspring_pheno_normal), 
                    nrow(all_phenotypes_normal)))
        
        # Build ID vector safely
        if(FALSE){
        all_genotypes <- do.call(rbind, genotypes_list)
        all_ids <- c(
        founders_ped$individual_id,
        all_pedigree %>% 
            filter(generation >= 1) %>%
            arrange(iteration, generation, family_id) %>%
            pull(individual_id)
        )

        # Create SNP names (AX-00000001 format)
        snp_names <- sprintf("AX-%08d", 1:ncol(all_genotypes))

        geno_df <- as_tibble(all_genotypes)
        colnames(geno_df) <- snp_names

        rm(all_genotypes)
        gc()

        # Add metadata columns
        geno_df <- geno_df %>%
        mutate(
            ID = all_ids,
            generation = c(
            rep(0, nInd(founders)),
            all_pedigree %>% filter(generation >= 1) %>% pull(generation)
            ),
            population_type = ifelse(generation == 0, "founder", "reference"),
            .before = 1
        )

        rm(geno_df)
        gc()

        write_csv(geno_df, file.path(iter_dir, sprintf("%s_gen%d_genotypes.csv", organism, gen)))
        }

        # Combine genotypes only when writing
        all_genotypes <- do.call(rbind, genotypes_list)

        # Build ID vector safely
        all_ids <- c(
          founders_ped$individual_id,
          all_pedigree %>% 
            filter(generation >= 1) %>%
            arrange(iteration, generation, family_id) %>%
            pull(individual_id)
        )

        # Build generation vector
        gen_vector <- c(
          rep(0, nInd(founders)),
          all_pedigree %>% filter(generation >= 1) %>% pull(generation)
        )

        # Build population_type vector
        pop_type_vector <- ifelse(gen_vector == 0, "founder", "reference")

        # Create SNP names (AX-00000001 format)
        snp_names <- sprintf("AX-%08d", 1:ncol(all_genotypes))

        # Write using fwrite (MUCH faster and less memory)
        geno_dt <- data.table(
          ID = all_ids,
          generation = gen_vector,
          population_type = pop_type_vector,
          all_genotypes
        )
        setnames(geno_dt, 4:ncol(geno_dt), snp_names)

        fwrite(geno_dt, 
              file.path(iter_dir, sprintf("%s_gen%d_genotypes.csv", organism, gen)),
              sep = ",")
        n_geno_rows <- nrow(geno_dt)
        rm(all_genotypes, geno_dt)
        gc()

        write_csv(all_phenotypes_normal, file.path(iter_dir, sprintf("%s_gen%d_phenotypes_normal.csv", organism, gen)))
        write_csv(all_phenotypes_mixture, file.path(iter_dir, sprintf("%s_gen%d_phenotypes_mixture.csv", organism, gen)))
        write_csv(all_pedigree, file.path(iter_dir, sprintf("%s_pedigree_all.csv", organism)))
        
        cat(sprintf("  ✓ Files saved to: %s\n", iter_dir))
        cat(sprintf("    - gen%d_genotypes.csv (%d rows)\n", gen, n_geno_rows))
        cat(sprintf("    - gen%d_phenotypes_normal.csv (%d rows)\n", gen, nrow(all_phenotypes_normal)))
        cat(sprintf("    - gen%d_phenotypes_mixture.csv (%d rows)\n", gen, nrow(all_phenotypes_mixture)))
        cat(sprintf("    - pedigree_all.csv (%d rows)\n", nrow(all_pedigree)))
        
        cat(sprintf("  ✓ Generation %d completed: %d new offspring, %d total individuals\n",
                          gen, nInd(offspring), n_geno_rows))

    }, error = function(e) {
      cat(sprintf("\n❌ ERROR in Iteration %d, Generation %d:\n", iter, gen))
      cat(sprintf("   Message: %s\n", e$message))
      cat(sprintf("   Call: %s\n", deparse(e$call)))
      stop(e)  # Re-throw error after logging
    })
  }

  cat("\n--- Validation Checks ---\n")
  
  # Check 1: TBV Variance by Generation (Normal scenario)
  tbv_var_by_gen_normal <- all_phenotypes_normal %>%
    group_by(generation) %>%
    summarize(
      mean_tbv = mean(tbv),
      sd_tbv = sd(tbv),
      var_tbv = var(tbv),
      .groups = "drop"
    )
  
  cat("\nTBV Statistics (Normal QTL) by Generation:\n")
  print(tbv_var_by_gen_normal)
  
  # Check 2: TBV Variance by Generation (Mixture scenario)
  tbv_var_by_gen_mixture <- all_phenotypes_mixture %>%
    group_by(generation) %>%
    summarize(
      mean_tbv = mean(tbv),
      sd_tbv = sd(tbv),
      var_tbv = var(tbv),
      .groups = "drop"
    )
  
  cat("\nTBV Statistics (Mixture QTL) by Generation:\n")
  print(tbv_var_by_gen_mixture)
  
  # Check 3: Heritability Realized by Generation
  h2_by_gen_normal <- all_phenotypes_normal %>%
    group_by(generation) %>%
    summarize(
      var_g = var(tbv),
      var_p = var(phenotype),
      h2_realized = var(tbv) / var(phenotype),
      .groups = "drop"
    )
  
  h2_by_gen_mixture <- all_phenotypes_mixture %>%
    group_by(generation) %>%
    summarize(
      var_g = var(tbv),
      var_p = var(phenotype),
      h2_realized = var(tbv) / var(phenotype),
      .groups = "drop"
    )
  
  cat("\nHeritability (Normal QTL) by Generation:\n")
  print(h2_by_gen_normal)
  
  cat("\nHeritability (Mixture QTL) by Generation:\n")
  print(h2_by_gen_mixture)
  
  # Check 4: Warning if h² deviates too much from target
  h2_normal_range <- range(h2_by_gen_normal$h2_realized)
  h2_mixture_range <- range(h2_by_gen_mixture$h2_realized)
  
  if(any(abs(h2_by_gen_normal$h2_realized - h2_target) > 0.05)) {
    warning(sprintf("Iteration %d: h² (normal) deviates from target (%.3f-%.3f vs %.2f target)",
                    iter, h2_normal_range[1], h2_normal_range[2], h2_target))
  }
  
  if(any(abs(h2_by_gen_mixture$h2_realized - h2_target) > 0.05)) {
    warning(sprintf("Iteration %d: h² (mixture) deviates from target (%.3f-%.3f vs %.2f target)",
                    iter, h2_mixture_range[1], h2_mixture_range[2], h2_target))
  }
  
  # Check 5: Genetic trend
  if(n_generations >= 2) {
    tbv_trend_normal <- diff(tbv_var_by_gen_normal$mean_tbv)
    tbv_trend_mixture <- diff(tbv_var_by_gen_mixture$mean_tbv)
    
    cat(sprintf("\nGenetic Trend (Mean TBV change per generation):\n"))
    cat(sprintf("  Normal:  %s\n", paste(sprintf("%.3f", tbv_trend_normal), collapse = ", ")))
    cat(sprintf("  Mixture: %s\n", paste(sprintf("%.3f", tbv_trend_mixture), collapse = ", ")))
  }
  
  # Save validation summary
  validation_summary <- tibble(
    iteration = iter,
    scenario = rep(c("normal", "mixture"), each = n_generations),
    generation = rep(1:n_generations, 2),
    mean_tbv = c(tbv_var_by_gen_normal$mean_tbv, tbv_var_by_gen_mixture$mean_tbv),
    sd_tbv = c(tbv_var_by_gen_normal$sd_tbv, tbv_var_by_gen_mixture$sd_tbv),
    var_tbv = c(tbv_var_by_gen_normal$var_tbv, tbv_var_by_gen_mixture$var_tbv),
    h2_realized = c(h2_by_gen_normal$h2_realized, h2_by_gen_mixture$h2_realized)
  )
  
  write_csv(validation_summary, file.path(iter_dir, sprintf("%s_validation_summary.csv", organism)))
  cat("\n✓ Validation summary saved\n")
  
  cat("\n==================================================\n")
  cat(sprintf("✓ ITERATION %d COMPLETED\n", iter))
  cat(sprintf("  Total individuals: %d (%d founders + %d offspring)\n",
              n_geno_rows, nInd(founders), nrow(all_phenotypes_normal)))
  cat(sprintf("  Files saved: %s\n", iter_dir))
  cat(sprintf("  Time: %s\n", Sys.time()))
  cat("\n==================================================\n")
 
}

cat("\n========================================\n")
cat("ALL SIMULATIONS COMPLETED\n")
cat(sprintf("Total iterations: %d\n", n_iterations))
cat(sprintf("Total generations per iteration: %d\n", n_generations))
cat("========================================\n")

# ============================================================================
# PRINT EXAMPLE DATA
# ============================================================================

cat("\n========================================\n")
cat("EXAMPLE DATA PREVIEW\n")
cat("========================================\n")

# Use iteration 1, generation 3 as example
example_iter <- 1
example_gen <- 3
example_dir <- file.path(output_dir, sprintf("iteration_%d", example_iter))

# Read example files
cat("\n--- GENOTYPE DATA ---\n")
example_geno <- read_csv(
  file.path(example_dir, sprintf("%s_gen%d_genotypes.csv", organism, example_gen)),
  show_col_types = FALSE
)
cat(sprintf("Dimensions: %d rows × %d columns\n", nrow(example_geno), ncol(example_geno)))
print(example_geno[1:5, 1:8])

cat("\n--- PHENOTYPE DATA (NORMAL) ---\n")
example_pheno_normal <- read_csv(
  file.path(example_dir, sprintf("%s_gen%d_phenotypes_normal.csv", organism, example_gen)),
  show_col_types = FALSE
)
cat(sprintf("Dimensions: %d rows × %d columns\n", nrow(example_pheno_normal), ncol(example_pheno_normal)))
print(head(example_pheno_normal, 5))

cat("\n--- PHENOTYPE DATA (MIXTURE) ---\n")
example_pheno_mixture <- read_csv(
  file.path(example_dir, sprintf("%s_gen%d_phenotypes_mixture.csv", organism, example_gen)),
  show_col_types = FALSE
)
cat(sprintf("Dimensions: %d rows × %d columns\n", nrow(example_pheno_mixture), ncol(example_pheno_mixture)))
print(head(example_pheno_mixture, 5))

cat("\n--- PEDIGREE DATA ---\n")
example_ped <- read_csv(
  file.path(example_dir, sprintf("%s_pedigree_all.csv", organism)),
  show_col_types = FALSE
)
cat(sprintf("Dimensions: %d rows × %d columns\n", nrow(example_ped), ncol(example_ped)))
cat("First 5 rows (founders):\n")
print(head(example_ped %>% filter(generation == 0), 5))
cat("First 5 rows (generation 1):\n")
print(head(example_ped %>% filter(generation == 1), 5))

cat("\n--- SNP MAP DATA ---\n")
example_map <- read_tsv(
  file.path(output_dir, sprintf("%s_snp_map.txt", organism)),
  show_col_types = FALSE
)
cat(sprintf("Dimensions: %d rows × %d columns\n", nrow(example_map), ncol(example_map)))
print(head(example_map, 10))

cat("\n--- QTL EFFECTS DATA ---\n")
example_qtl <- read_csv(
  file.path(example_dir, sprintf("%s_qtl_effects.csv", organism)),
  show_col_types = FALSE
)
cat(sprintf("Dimensions: %d rows × %d columns\n", nrow(example_qtl), ncol(example_qtl)))
print(head(example_qtl, 10))

cat("\n--- SUMMARY STATISTICS ---\n")
cat(sprintf("Total founders: %d\n", sum(example_ped$generation == 0)))
cat(sprintf("Total offspring (all generations): %d\n", sum(example_ped$generation > 0)))
cat(sprintf("Total SNPs: %d\n", nrow(example_map)))
cat(sprintf("Total QTLs: %d\n", nrow(example_qtl)))
cat(sprintf("\nPhenotype correlation (normal vs mixture): %.3f\n",
            cor(example_pheno_normal$tbv, example_pheno_mixture$tbv)))

cat("\n========================================\n")
cat("END OF SIMULATION\n")
cat("========================================\n")