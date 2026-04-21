# data_prep.R
# Data preprocessing and feature engineering functions

suppressPackageStartupMessages({
  library(tidyverse)
})

load_genomic_data <- function(config) {
  genotype_ref <- read_csv(file.path(config$data_dir, config$genotype), show_col_types = FALSE)
  
  phenotype_ref <- read_csv(file.path(config$data_dir, config$phenotype), show_col_types = FALSE) %>%
    filter(generation == config$target_generation)
  
  pedigree <- read_csv(file.path(config$data_dir, config$pedigree), show_col_types = FALSE) %>%
    filter(
      population_type == config$population_type,
      generation == config$target_generation
    )
  
  # RENAME ID to individual_id in genotype for consistency
  if("ID" %in% names(genotype_ref) && !"individual_id" %in% names(genotype_ref)) {
    genotype_ref <- genotype_ref %>% rename(individual_id = ID)
    cat("Renamed genotype ID column to individual_id\n")
  }
  
  cat("\n=== Data Loading Diagnostics ===\n")
  cat("Genotype rows:", nrow(genotype_ref), "\n")
  cat("Genotype first 3 IDs:", paste(head(genotype_ref$individual_id, 3), collapse=", "), "\n")
  
  cat("\nPhenotype rows:", nrow(phenotype_ref), "\n")
  cat("Phenotype first 3 IDs:", paste(head(phenotype_ref$individual_id, 3), collapse=", "), "\n")
  
  cat("\nPedigree rows (after filter):", nrow(pedigree), "\n")
  
  # Join pedigree info
  phenotype_ref <- phenotype_ref %>%
    left_join(pedigree %>% select(individual_id, father_id, mother_id), 
              by = "individual_id")
  
  # FILTER and REORDER genotype to match phenotype
  cat("\n=== Matching Genotype to Phenotype ===\n")
  cat("Phenotype IDs in genotype:", 
      sum(phenotype_ref$individual_id %in% genotype_ref$individual_id), "/", nrow(phenotype_ref), "\n")
  
  # Keep only individuals in phenotype and reorder
  genotype_ref <- genotype_ref %>%
    filter(individual_id %in% phenotype_ref$individual_id) %>%
    slice(match(phenotype_ref$individual_id, individual_id))
  
  cat("Genotype rows after filtering:", nrow(genotype_ref), "\n")
  cat("Order matches:", identical(phenotype_ref$individual_id, genotype_ref$individual_id), "\n")
  cat("================================\n\n")
  
  list(
    genotype = genotype_ref,
    phenotype = phenotype_ref,
    ids = genotype_ref$individual_id
  )
}

quality_control <- function(geno_matrix, config) {
  if(!is.numeric(geno_matrix)) {
    stop("Genotype matrix must be numeric. Check if ID column was removed.")
  }
  
  allele_freq <- colMeans(geno_matrix, na.rm = TRUE) / 2
  maf <- pmin(allele_freq, 1 - allele_freq)
  callrate <- colMeans(!is.na(geno_matrix))
  
  keep_snps <- (maf >= config$min_allele_freq) & 
               (callrate >= config$callrate_threshold)
  
  cat("=== Quality Control ===\n")
  cat("SNPs retained:", sum(keep_snps), "out of", ncol(geno_matrix), "\n")
  
  list(
    filtered = geno_matrix[, keep_snps],
    keep_snps = keep_snps,
    n_retained = sum(keep_snps)
  )
}

impute_and_standardize <- function(geno_filtered) {
  geno_imputed <- apply(geno_filtered, 2, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    x
  })

  snp_means <- colMeans(geno_imputed)
  snp_sds <- apply(geno_imputed, 2, sd)
  zero_var <- snp_sds < 1e-10

  if(sum(zero_var) > 0) {
    cat("Removing", sum(zero_var), "zero-variance SNPs\n")
    geno_imputed <- geno_imputed[, !zero_var]
    snp_means <- snp_means[!zero_var]
    snp_sds <- snp_sds[!zero_var]
  }
  
  geno_std <- scale(geno_imputed, center = snp_means, scale = snp_sds)

  if(any(is.na(geno_std))) {
    stop("Standardization produced NA values - check for zero variance SNPs")
  }
  
  cat("Clean matrix:", ncol(geno_std), "SNPs\n\n")
  
  list(
    standardized = geno_std,
    imputed = geno_imputed
  )
}

create_grm <- function(geno_imputed) {
  allele_freq <- colMeans(geno_imputed, na.rm = TRUE) / 2
  Z <- geno_imputed - matrix(2 * allele_freq, 
                             nrow = nrow(geno_imputed), 
                             ncol = ncol(geno_imputed), byrow = TRUE)
  k <- 2 * sum(allele_freq * (1 - allele_freq))
  GRM <- tcrossprod(Z) / k

  high_rel <- which(GRM > 0.99 & lower.tri(GRM), arr.ind = TRUE)
  if(nrow(high_rel) > 0) {
    cat("WARNING: Found", nrow(high_rel), "pairs with GRM > 0.99 (possible duplicates)\n")
    cat("First 5 pairs:\n")
    print(head(high_rel, 5))
  }

  bad_diag <- which(abs(diag(GRM) - 1) > 0.1)
  if(length(bad_diag) > 0) {
    cat("WARNING:", length(bad_diag), "individuals with diagonal != 1 ± 0.1\n")
    cat("Range:", range(diag(GRM)[bad_diag]), "\n")
  }

  # GRM evaluation
  # Diagonal (self-relationship)
  diag_mean <- mean(diag(GRM))
  diag_range <- range(diag(GRM))

  # Off-diagonal (relatedness)
  offdiag <- GRM[lower.tri(GRM)]
  offdiag_mean <- mean(offdiag)
  offdiag_range <- range(offdiag)
  
  cat("=== GRM Evaluation ===\n")
  cat("Diagonal mean:", round(diag_mean, 3), "| range: [", 
      round(diag_range[1], 3), ",", round(diag_range[2], 3), "]\n")
  cat("Off-diagonal mean:", round(offdiag_mean, 3), "| range: [", 
      round(offdiag_range[1], 3), ",", round(offdiag_range[2], 3), "]\n\n")

  GRM
}

setup_cv_folds <- function(phenotype_ref, config) {
  set.seed(config$seed)
  
  family_check <- phenotype_ref %>%
    count(father_id) %>%
    summarise(
      n_sires = n(),
      median_offspring = median(n),
      sires_with_5plus = sum(n >= 5)
    )
  
  min_off   <- if(!is.null(config$min_offspring_per_sire)) config$min_offspring_per_sire else 5
  cv_method <- config$cv_method
  feasible  <- family_check$n_sires > 1 && family_check$median_offspring >= min_off

  if(cv_method %in% c("within_family", "across_family") && !feasible) {
    cv_method <- "random"
    reason    <- if(family_check$n_sires <= 1) "only_one_sire"
                 else if(family_check$median_offspring < min_off) "insufficient_offspring"
                 else "too_few_adequate_sires"
  } else {
    reason <- NULL
  }

  if(cv_method == "within_family") {
    folds          <- caret::createFolds(y=phenotype_ref$father_id, k=config$k_folds, list=FALSE, returnTrain=FALSE)
    sires_per_fold <- phenotype_ref %>% mutate(fold=folds) %>%
                      group_by(fold) %>% summarise(n_sires=n_distinct(father_id)) %>% pull(n_sires)
    avg_folds_per_sire <- mean(table(folds))
    method <- "within_family"

  } else if(cv_method == "across_family") {
    sire_ids           <- unique(phenotype_ref$father_id)
    sire_folds         <- setNames(sample(rep(1:config$k_folds, length.out=length(sire_ids))), sire_ids)
    folds              <- sire_folds[phenotype_ref$father_id]
    sires_per_fold     <- phenotype_ref %>% mutate(fold=folds) %>%
                          group_by(fold) %>% summarise(n_sires=n_distinct(father_id)) %>% pull(n_sires)
    avg_folds_per_sire <- mean(table(folds))
    method <- "across_family"

  } else {
    folds              <- sample(rep(1:config$k_folds, length.out=nrow(phenotype_ref)))
    sires_per_fold     <- NULL
    avg_folds_per_sire <- NULL
    method             <- "random"
    reason             <- if(is.null(reason)) "cv_method_random" else reason
  }
  
  list(
    folds = folds,
    method = method,
    reason = reason,
    sires_per_fold = sires_per_fold,
    avg_folds_per_sire = avg_folds_per_sire
  )
}

remove_duplicates <- function(geno_imputed, phenotype_ref, config) {
  threshold <- config$duplicate_threshold
  allele_freq <- colMeans(geno_imputed) / 2
  Z <- geno_imputed - matrix(2 * allele_freq, 
                             nrow = nrow(geno_imputed), 
                             ncol = ncol(geno_imputed), byrow = TRUE)
  k <- 2 * sum(allele_freq * (1 - allele_freq))
  GRM_temp <- tcrossprod(Z) / k
  
  high_rel <- which(GRM_temp > threshold & lower.tri(GRM_temp), arr.ind = TRUE)
  
  if(nrow(high_rel) == 0) {
    return(list(
      keep_idx = 1:nrow(geno_imputed),
      n_removed = 0
    ))
  }
  
  to_remove <- integer(0)
  for(i in 1:nrow(high_rel)) {
    idx1 <- high_rel[i, 1]
    idx2 <- high_rel[i, 2]
    if(!(idx1 %in% to_remove) && !(idx2 %in% to_remove)) {
      to_remove <- c(to_remove, idx2)
    }
  }
  
  to_remove <- unique(to_remove)
  cat("Removing", length(to_remove), "duplicates (GRM >", threshold, "):",
  paste(phenotype_ref$individual_id[to_remove], collapse = ", "), "\n")
  
  list(
    keep_idx = setdiff(1:nrow(geno_imputed), to_remove),
    n_removed = length(to_remove)
  )
}