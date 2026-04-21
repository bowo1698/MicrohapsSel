suppressPackageStartupMessages({
  library(tidyverse)
  library(sommer)
  library(coda)
  library(masbayes)
  library(RhpcBLASctl)
})

select <- dplyr::select
filter <- dplyr::filter
blas_set_num_threads(1)

config_file <- Sys.getenv("GENOMIC_CONFIG", "config.R")
source(config_file)
setwd(config$base_dir)

source("R/data_preparation.R")
source("R/matrix_construction.R")
source("R/evaluation.R")
source("R/models/gblup.R")
source("R/models/bayesR.R")
source("R/models/bayesA.R")

pipeline_start <- Sys.time()
cat(sprintf("=== Pipeline started: %s ===\n", Sys.time()))

# ── 1. Load genotype ONCE (shared semua skenario) ─────────────────────────────
hap_geno_list <- lapply(1:config$n_chromosomes, function(chr) {
  read_table(
    file.path(config$data_dir, config$block_dir, paste0("hap_geno_", chr)),
    show_col_types = FALSE
  )
})
hap_geno <- hap_geno_list[[1]]
for (chr in 2:config$n_chromosomes)
  hap_geno <- bind_cols(hap_geno, hap_geno_list[[chr]] %>% select(-ID))
hap_geno <- hap_geno %>% rename(individual_id = ID)
cat(sprintf("Genotype loaded: %d individuals, %d chromosomes\n", nrow(hap_geno), config$n_chromosomes))

pedigree <- read_csv(
  file.path(config$data_dir, config$pedigree_file),
  show_col_types = FALSE
) %>% filter(population_type == config$population_type,
             generation      == config$target_generation)

# ── 2. Construct W matrix ONCE ────────────────────────────────────────────────
hap_cols <- hap_geno %>% select(-individual_id)
allele_freq <- compute_allele_frequencies(hap_cols, config)

W_result <- construct_W_matrix(hap_cols, allele_freq, drop_baseline = config$drop_baseline)

# Filter zero-variance + scale
col_var      <- apply(W_result$W_ah, 2, var)
keep         <- col_var > 0
W_result$W_ah <- W_result$W_ah[, keep]
if (!is.null(W_result$allele_info))
  W_result$allele_info <- lapply(W_result$allele_info, function(x) x[keep])

scale_factor  <- sqrt(mean(apply(W_result$W_ah, 2, var)))
W_result$W_ah <- W_result$W_ah / scale_factor

W <- W_result$W_ah
n <- nrow(W)

# GRM
G     <- tcrossprod(W)
k_ah  <- sum(diag(G)) / n
A_mat <- G / k_ah
rownames(A_mat) <- colnames(A_mat) <- as.character(hap_geno$individual_id)
WtW_diag <- colSums(W^2)
cat(sprintf("W matrix ready: %d ind x %d alleles | k_ah=%.4f | scale_factor=%.4f\n",
            nrow(W), ncol(W), k_ah, scale_factor))

# ── 3. Loop skenario ──────────────────────────────────────────────────────────
for (sc in names(config$scenarios)) {

  pheno_file <- config$scenarios[[sc]]
  out_dir    <- file.path(config$output_dir, sc)
  sc_start   <- Sys.time()

  pheno <- read_csv(file.path(config$data_dir, pheno_file), show_col_types = FALSE) %>%
    filter(generation == config$target_generation) %>%
    left_join(pedigree %>% select(individual_id, father_id, mother_id),
              by = "individual_id")

  common_ids <- intersect(pheno$individual_id, hap_geno$individual_id)
  pheno <- pheno %>%
    filter(individual_id %in% common_ids) %>%
    arrange(match(individual_id, common_ids))

  idx   <- match(pheno$individual_id, hap_geno$individual_id)
  A_sc <- A_mat[idx, idx]
  W_sc  <- W[idx, ]

  Wty <- crossprod(W_sc, pheno$phenotype)

  matrices <- list(
    W_train         = W_sc,
    W_test          = W_sc[0, ],   # empty — tidak dipakai saat fit
    allele_info     = W_result$allele_info,
    dropped_alleles = W_result$dropped_alleles,
    A_train         = A_sc,
    A_combined      = A_sc,
    k_ah            = k_ah,
    WtW_diag        = WtW_diag,
    Wty             = Wty
  )

  split <- list(
    train = list(phenotypes = pheno,
                 haplotypes = hap_geno %>%
                   filter(individual_id %in% pheno$individual_id) %>%
                   arrange(match(individual_id, pheno$individual_id))),
    test  = list(phenotypes = pheno[0, ],
                 haplotypes = hap_geno[0, ])
  )

  cat(sprintf("\n[%s] Scenario: %s\n", sc, Sys.time()))
  cat(sprintf("  N individuals: %d | N alleles: %d\n", length(idx), ncol(W_sc)))

  # ── GBLUP ──
  cat("  Running GBLUP...\n")
  dir.create(file.path(out_dir, "gblup"), recursive = TRUE, showWarnings = FALSE)

  gblup_res <- run_gblup(matrices, split, config)

  gblup_model <- list(
    mu_hat          = mean(pheno$phenotype),
    varcomp         = gblup_res$varcomp,
    k_ah            = k_ah,
    A_mat           = A_sc,
    A_ref_inv       = solve(A_sc + diag(gblup_res$varcomp$lambda, length(idx))),
    W_train         = W_sc,
    GEBV_train      = gblup_res$predictions$train,
    allele_freq     = allele_freq,
    W_ref_structure = list(
      allele_info     = W_result$allele_info,
      dropped_alleles = W_result$dropped_alleles
    ),
    scale_factor    = scale_factor,
    keep_cols       = keep,
    ind_ids         = as.character(pheno$individual_id)
  )
  saveRDS(gblup_model, file.path(out_dir, "gblup", "model.rds"))

  gebv_gblup <- tibble(
    individual_id = pheno$individual_id,
    tbv           = pheno$tbv,
    GEBV          = gblup_res$predictions$train
  ) %>% arrange(desc(GEBV)) %>% mutate(rank = row_number())
  saveRDS(gebv_gblup, file.path(out_dir, "gblup", "gebv_ranked.rds"))

  cat(sprintf("  GBLUP done | h2=%.3f | cor_train_tbv=%.3f\n",
              gblup_res$varcomp$h2, gblup_res$accuracy$train_tbv))

  # ── BayesR ──
  cat("  Running BayesR...\n")
  dir.create(file.path(out_dir, "bayesR"), recursive = TRUE, showWarnings = FALSE)
  bayesR_res <- run_bayesR(matrices, split, gblup_res$varcomp, config, fold = 0)

  bayesR_model <- list(
    beta_hat        = bayesR_res$coef$beta_hat,
    mu_hat          = bayesR_res$coef$mu_hat,
    varcomp         = bayesR_res$varcomp,
    markers         = bayesR_res$markers,
    samples         = bayesR_res$samples,
    allele_freq     = allele_freq,
    W_ref_structure = list(
      allele_info     = W_result$allele_info,
      dropped_alleles = W_result$dropped_alleles
    ),
    scale_factor    = scale_factor,
    keep_cols       = keep,
    ind_ids         = as.character(pheno$individual_id)
  )
  saveRDS(bayesR_model, file.path(out_dir, "bayesR", "model.rds"))

  gebv_bayesR <- tibble(
    individual_id = pheno$individual_id,
    tbv           = pheno$tbv,
    GEBV          = bayesR_res$predictions$train
  ) %>% arrange(desc(GEBV)) %>% mutate(rank = row_number())
  saveRDS(gebv_bayesR, file.path(out_dir, "bayesR", "gebv_ranked.rds"))

  cat(sprintf("  BayesR done | h2=%.3f | cor_train_tbv=%.3f\n",
              bayesR_res$varcomp$h2, bayesR_res$accuracy$train_tbv))

  # ── BayesA ──
  cat("  Running BayesA...\n")
  dir.create(file.path(out_dir, "bayesA"), recursive = TRUE, showWarnings = FALSE)
  bayesA_res <- run_bayesA(matrices, split, gblup_res$varcomp, config, fold = 0)

  bayesA_model <- list(
    beta_hat        = bayesA_res$coef$beta_hat,
    mu_hat          = bayesA_res$coef$mu_hat,
    varcomp         = bayesA_res$varcomp,
    markers         = bayesA_res$markers,
    samples         = bayesA_res$samples,
    allele_freq     = allele_freq,
    W_ref_structure = list(
      allele_info     = W_result$allele_info,
      dropped_alleles = W_result$dropped_alleles
    ),
    scale_factor    = scale_factor,
    keep_cols       = keep,
    ind_ids         = as.character(pheno$individual_id)
  )
  saveRDS(bayesA_model, file.path(out_dir, "bayesA", "model.rds"))

  gebv_bayesA <- tibble(
    individual_id = pheno$individual_id,
    tbv           = pheno$tbv,
    GEBV          = bayesA_res$predictions$train
  ) %>% arrange(desc(GEBV)) %>% mutate(rank = row_number())
  saveRDS(gebv_bayesA, file.path(out_dir, "bayesA", "gebv_ranked.rds"))

  cat(sprintf("  BayesA done | h2=%.3f | cor_train_tbv=%.3f\n",
              bayesA_res$varcomp$h2, bayesA_res$accuracy$train_tbv))

  cat(sprintf("Scenario %s done in %.1f mins\n", sc,
              as.numeric(difftime(Sys.time(), sc_start, units = "mins"))))
}

cat(sprintf("\nAll scenarios complete | Total: %.1f mins\n",
            as.numeric(difftime(Sys.time(), pipeline_start, units = "mins"))))