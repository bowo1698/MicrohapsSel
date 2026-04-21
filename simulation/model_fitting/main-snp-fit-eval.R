suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(sommer)
  library(masbayes)
  library(RhpcBLASctl)
})

select <- dplyr::select
filter <- dplyr::filter
blas_set_num_threads(1)

config_file <- Sys.getenv("GENOMIC_CONFIG", "config.R")
source(config_file)
setwd(config$base_dir)

source("R-snp/data_prep.R")
source("R-snp/models.R")
source("R-snp/evaluation.R")

pipeline_start <- Sys.time()
cat(sprintf("=== Pipeline started: %s ===\n", Sys.time()))

# ── 1. Load genotype ONCE ─────────────────────────────────────────────────────
cat("Loading genotype...\n")
genotype_raw <- as_tibble(fread(file.path(config$data_dir, config$genotype)))
if ("ID" %in% names(genotype_raw) && !"individual_id" %in% names(genotype_raw))
  genotype_raw <- genotype_raw %>% rename(individual_id = ID)
cat(sprintf("Genotype loaded: %d individuals x %d columns\n", nrow(genotype_raw), ncol(genotype_raw)))

pedigree <- read_csv(file.path(config$data_dir, config$pedigree), show_col_types = FALSE) %>%
  filter(population_type == config$population_type,
         generation      == config$target_generation)

# ── 2. QC + standardize ONCE ─────────────────────────────────────────────────
geno_matrix <- genotype_raw %>% select(-individual_id, -any_of(c("generation", "population_type"))) %>% as.matrix()
ind_ids_raw <- genotype_raw$individual_id

qc_result   <- quality_control(geno_matrix, config)
geno_proc   <- impute_and_standardize(qc_result$filtered)

# Remove duplicates
dedup <- remove_duplicates(geno_proc$imputed, 
                           data.frame(individual_id = ind_ids_raw), 
                           config)
keep_idx     <- dedup$keep_idx
geno_std     <- geno_proc$standardized[keep_idx, ]
geno_imputed <- geno_proc$imputed[keep_idx, ]
ind_ids      <- ind_ids_raw[keep_idx]

# Standardization params — needed for predict
snp_means <- colMeans(geno_imputed)
snp_sds   <- apply(geno_imputed, 2, sd)
keep_snps <- qc_result$keep_snps

# GRM full (semua individu setelah QC)
GRM_full <- create_grm(geno_imputed)
rownames(GRM_full) <- colnames(GRM_full) <- as.character(ind_ids)

# W matrix (VanRaden centering)
allele_freq_full <- colMeans(geno_imputed) / 2
W_full <- sweep(geno_imputed, 2, 2 * allele_freq_full, "-")
storage.mode(W_full) <- "double"

cat(sprintf("Genotype ready: %d individuals | %d SNPs retained | %d duplicates removed\n",
            length(ind_ids), qc_result$n_retained, dedup$n_removed))

# ── 3. Loop skenario ──────────────────────────────────────────────────────────
for (sc in names(config$scenarios)) {

  pheno_file <- config$scenarios[[sc]]
  out_dir    <- file.path(config$output_dir, sc)
  sc_start   <- Sys.time()

  pheno <- read_csv(file.path(config$data_dir, pheno_file), show_col_types = FALSE) %>%
    filter(generation == config$target_generation) %>%
    left_join(pedigree %>% select(individual_id, father_id, mother_id), by = "individual_id")

  common_ids <- intersect(pheno$individual_id, ind_ids)
  pheno <- pheno %>%
    filter(individual_id %in% common_ids) %>%
    arrange(match(individual_id, common_ids))

  idx      <- match(pheno$individual_id, ind_ids)
  geno_sc  <- geno_std[idx, ]
  GRM_sc   <- GRM_full[idx, idx]
  y        <- pheno$phenotype
  train_ids <- as.character(pheno$individual_id)

  cat(sprintf("\n[%s] Scenario: %s\n", sc, Sys.time()))
  cat(sprintf("  N individuals: %d | N SNPs: %d\n", length(idx), ncol(geno_sc)))

  # ── GBLUP ──
  cat("  Running GBLUP...\n")
  dir.create(file.path(out_dir, "gblup"), recursive = TRUE, showWarnings = FALSE)
  gblup_fit <- train_gblup(y, GRM_sc, n_threads = 1)

  var_g  <- gblup_fit$var_g
  var_e  <- gblup_fit$var_e
  lambda <- var_e / var_g
  GRM_inv <- solve(GRM_sc + diag(lambda, nrow(GRM_sc)))
  pred_train_gblup <- predict_gblup_train(gblup_fit, GRM_sc, y)$pred

  # Compute u_snp for cross-generation prediction
  allele_freq_sc <- snp_means[keep_snps] / 2
  k_sc           <- 2 * sum(allele_freq_sc * (1 - allele_freq_sc))
  Z_train        <- sweep(geno_imputed[idx, ], 2, 2 * allele_freq_sc, "-")
  u_snp          <- as.vector(crossprod(Z_train, GRM_inv %*% gblup_fit$u)) / k_sc

  gblup_model <- list(
    fit         = gblup_fit$fit,
    mu          = as.numeric(gblup_fit$fit$b),
    var_g       = var_g,
    var_e       = var_e,
    lambda      = lambda,
    GRM_train   = GRM_sc,
    GRM_inv     = GRM_inv,
    u_train     = gblup_fit$u,
    u_snp       = u_snp,
    allele_freq = allele_freq_sc,
    k           = k_sc,
    snp_means   = snp_means,
    snp_sds     = snp_sds,
    keep_snps   = keep_snps,
    ind_ids     = train_ids
  )
  saveRDS(gblup_model, file.path(out_dir, "gblup", "model.rds"))

  gebv_gblup <- tibble(
    individual_id = pheno$individual_id,
    tbv           = pheno$tbv,
    GEBV          = pred_train_gblup
  ) %>% arrange(desc(GEBV)) %>% mutate(rank = row_number())
  saveRDS(gebv_gblup, file.path(out_dir, "gblup", "gebv_ranked.rds"))

  cat(sprintf("  GBLUP done | h2=%.3f | cor_train_tbv=%.3f\n",
              gblup_fit$h2, cor(pred_train_gblup, pheno$tbv, use="complete.obs")))
  
  # ── BayesR ──
  cat("  Running BayesR...\n")
  dir.create(file.path(out_dir, "bayesR"), recursive = TRUE, showWarnings = FALSE)
  bayesR_fit <- train_bayesR(W_sc, y, var_g, var_e, config)

  bayesR_model <- list(
    beta_hat    = bayesR_fit$beta_hat,
    mu_hat      = bayesR_fit$mu_hat,
    markers     = bayesR_fit$markers,
    samples     = bayesR_fit$samples,
    varcomp     = list(sigma2_e = bayesR_fit$sigma2_e_hat,
                      sigma2_g = bayesR_fit$sigma2_g_hat,
                      h2       = bayesR_fit$h2),
    allele_freq = allele_freq_full,
    keep_snps   = keep_snps,
    ind_ids     = train_ids
  )
  saveRDS(bayesR_model, file.path(out_dir, "bayesR", "model.rds"))

  pred_train_bayesR <- bayesR_fit$pred_train
  gebv_bayesR <- tibble(
    individual_id = pheno$individual_id,
    tbv           = pheno$tbv,
    GEBV          = pred_train_bayesR
  ) %>% arrange(desc(GEBV)) %>% mutate(rank = row_number())
  saveRDS(gebv_bayesR, file.path(out_dir, "bayesR", "gebv_ranked.rds"))

  cat(sprintf("  BayesR done | h2=%.3f | cor_train_tbv=%.3f\n",
              bayesR_fit$h2, cor(pred_train_bayesR, pheno$tbv, use="complete.obs")))

  # ── BayesA ──
  cat("  Running BayesA...\n")
  dir.create(file.path(out_dir, "bayesA"), recursive = TRUE, showWarnings = FALSE)
  W_sc <- W_full[idx, ]
  bayesA_fit <- train_bayesA(W_sc, y, var_g, var_e, config)

  bayesA_model <- list(
    beta_hat    = bayesA_fit$beta_hat,
    mu_hat      = bayesA_fit$mu_hat,
    markers     = bayesA_fit$markers,
    samples     = bayesA_fit$samples,
    varcomp     = list(sigma2_e = bayesA_fit$sigma2_e_hat,
                      sigma2_g = bayesA_fit$sigma2_g_hat,
                      h2       = bayesA_fit$h2),
    allele_freq = allele_freq_full,
    keep_snps   = keep_snps,
    ind_ids     = train_ids
  )
  saveRDS(bayesA_model, file.path(out_dir, "bayesA", "model.rds"))

  pred_train_bayesA <- bayesA_fit$pred_train
  gebv_bayesA <- tibble(
    individual_id = pheno$individual_id,
    tbv           = pheno$tbv,
    GEBV          = pred_train_bayesA
  ) %>% arrange(desc(GEBV)) %>% mutate(rank = row_number())
  saveRDS(gebv_bayesA, file.path(out_dir, "bayesA", "gebv_ranked.rds"))

  cat(sprintf("  BayesA done | h2=%.3f | cor_train_tbv=%.3f\n",
              bayesA_fit$h2, cor(pred_train_bayesA, pheno$tbv, use="complete.obs")))

  # ── Metrics ──
  metrics_gblup <- data.frame(
    model    = "GBLUP",
    scenario = sc,
    cor_tbv  = cor(pred_train_gblup, pheno$tbv, use = "complete.obs"),
    bias_tbv = coef(lm(pheno$tbv ~ pred_train_gblup))[2],
    rmse_tbv = sqrt(mean((pred_train_gblup - pheno$tbv)^2, na.rm = TRUE)),
    h2       = gblup_fit$h2
  )

  metrics_bayesR <- data.frame(
    model    = "BayesR",
    scenario = sc,
    cor_tbv  = cor(pred_train_bayesR, pheno$tbv, use = "complete.obs"),
    bias_tbv = coef(lm(pheno$tbv ~ pred_train_bayesR))[2],
    rmse_tbv = sqrt(mean((pred_train_bayesR - pheno$tbv)^2, na.rm = TRUE)),
    h2       = bayesR_fit$h2
  )

  metrics_bayesA <- data.frame(
    model    = "BayesA",
    scenario = sc,
    cor_tbv  = cor(pred_train_bayesA, pheno$tbv, use = "complete.obs"),
    bias_tbv = coef(lm(pheno$tbv ~ pred_train_bayesA))[2],
    rmse_tbv = sqrt(mean((pred_train_bayesA - pheno$tbv)^2, na.rm = TRUE)),
    h2       = bayesA_fit$h2
  )

  metrics_all <- bind_rows(metrics_gblup, metrics_bayesA, metrics_bayesR)
  write.csv(metrics_all, file.path(out_dir, "train_metrics.csv"), row.names = FALSE)

  cat(sprintf("Scenario %s done in %.1f mins\n", sc,
              as.numeric(difftime(Sys.time(), sc_start, units = "mins"))))
}

cat(sprintf("\nAll scenarios complete | Total: %.1f mins\n",
            as.numeric(difftime(Sys.time(), pipeline_start, units = "mins"))))