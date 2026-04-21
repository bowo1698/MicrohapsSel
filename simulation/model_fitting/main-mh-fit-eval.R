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
cat(sprintf("=== CV Pipeline started: %s ===\n", Sys.time()))

# ── 1. Load genotype
load_hap_geno <- function(data_dir, n_chr, block_dir) {
  geno_list <- lapply(1:n_chr, function(chr) {
    read_table(
      file.path(data_dir, block_dir, paste0("hap_geno_", chr)),
      show_col_types = FALSE
    )
  })
  geno <- geno_list[[1]]
  for (chr in 2:n_chr)
    geno <- bind_cols(geno, geno_list[[chr]] %>% select(-ID))
  geno %>% rename(individual_id = ID)
}

hap_geno_current <- load_hap_geno(config$data_dir, config$n_chromosomes, config$block_dir)

cat(sprintf("Genotype loaded: %d individuals\n", nrow(hap_geno_current)))

# ── 2. Split per family
pedigree_current <- read_csv(
  file.path(config$data_dir, config$pedigree_file),
  show_col_types = FALSE
) %>% filter(population_type == config$population_type)

set.seed(config$seed_split)
split_ids <- pedigree_current %>%
  group_by(family_id) %>%
  summarise(
    train_ids = list(sample(individual_id, round(n() * config$train_prop))),
    test_ids  = list(setdiff(individual_id, train_ids[[1]])),
    .groups   = "drop"
  )

current_train_ids <- unlist(split_ids$train_ids)
test_ids          <- unlist(split_ids$test_ids)

train_ids <- current_train_ids

cat(sprintf("Split: train=%d | test=%d\n", length(train_ids), length(test_ids)))

# ── 3. Construct W matrix 
hap_train <- hap_geno_current %>%
  filter(individual_id %in% train_ids) %>%
  arrange(match(individual_id, train_ids))
hap_test <- hap_geno_current %>%
  filter(individual_id %in% test_ids) %>%
  arrange(match(individual_id, test_ids))

allele_freq <- compute_allele_frequencies(hap_train %>% select(-individual_id), config)

W_result_train <- construct_W_matrix(hap_train %>% select(-individual_id),
                                     allele_freq, drop_baseline = config$drop_baseline)
W_result_test  <- construct_W_matrix(hap_test %>% select(-individual_id),
                                     allele_freq_filtered = NULL,
                                     reference_structure  = W_result_train,
                                     drop_baseline        = config$drop_baseline)

col_var <- apply(W_result_train$W_ah, 2, var)
keep    <- col_var > 0
W_result_train$W_ah <- W_result_train$W_ah[, keep]
W_result_test$W_ah  <- W_result_test$W_ah[, keep]

scale_factor        <- sqrt(mean(apply(W_result_train$W_ah, 2, var)))
W_train <- W_result_train$W_ah / scale_factor
W_test  <- W_result_test$W_ah  / scale_factor
rownames(W_train) <- train_ids
rownames(W_test)  <- test_ids
n_train <- nrow(W_train)

G_train  <- tcrossprod(W_train)
k_ah     <- sum(diag(G_train)) / n_train
A_train  <- G_train / k_ah

W_combined <- rbind(W_train, W_test)
G_combined <- tcrossprod(W_combined)
A_combined <- G_combined / k_ah
all_ids    <- c(train_ids, test_ids)
rownames(A_combined) <- colnames(A_combined) <- all_ids

allele_info_keep <- lapply(W_result_train$allele_info, function(x) x[keep])
cat(sprintf("W matrix ready: %d train x %d alleles | k_ah=%.4f | scale_factor=%.4f\n",
            n_train, ncol(W_train), k_ah, scale_factor))

# ── 4. Loop scenario 
for (sc in names(config$scenarios)) {

  out_dir    <- file.path(config$output_dir, sc)
  sc_start   <- Sys.time()

  # Load pheno
  pheno_all <- read_csv(
    file.path(config$data_dir, config$scenarios[[sc]]),
    show_col_types = FALSE
  )

  pheno_train <- pheno_all %>%
    filter(individual_id %in% train_ids) %>%
    arrange(match(individual_id, train_ids))

  pheno_test <- pheno_all %>%
    filter(individual_id %in% test_ids) %>%
    arrange(match(individual_id, test_ids))

  cat(sprintf("\n[%s] Scenario: %s | train=%d | test=%d\n",
              sc, Sys.time(), nrow(pheno_train), nrow(pheno_test)))

  WtW_diag <- colSums(W_train^2)
  Wty      <- crossprod(W_train, pheno_train$phenotype)

  matrices <- list(
    W_train         = W_train,
    W_test          = W_test,
    allele_info     = allele_info_keep,
    dropped_alleles = W_result_train$dropped_alleles,
    A_train         = A_train,
    A_combined      = A_combined,
    k_ah            = k_ah,
    WtW_diag        = WtW_diag,
    Wty             = Wty
  )

  split <- list(
    train = list(
      phenotypes = pheno_train,
      haplotypes = hap_train
    ),
    test = list(
      phenotypes = pheno_test,
      haplotypes = hap_test
    )
  )

  # Run GBLUP once for all
  cat("  Running GBLUP...\n")
  gblup_res <- run_gblup(matrices, split, config)
  gblup_vc  <- gblup_res$varcomp

  for (model_name in c("gblup", "bayesR", "bayesA")) {

    dir.create(file.path(out_dir, model_name), recursive = TRUE, showWarnings = FALSE)

    if (model_name == "gblup") {
      GEBV_test <- gblup_res$predictions$test
      h2        <- gblup_vc$h2

    } else if (model_name == "bayesR") {
      cat("  Running BayesR...\n")
      res       <- run_bayesR(matrices, split, gblup_vc, config, fold = 0)
      GEBV_test <- res$predictions$test
      h2        <- res$varcomp$h2

    } else {
      cat("  Running BayesA...\n")
      res       <- run_bayesA(matrices, split, gblup_vc, config, fold = 0)
      GEBV_test <- res$predictions$test
      h2        <- res$varcomp$h2
    }

    acc  <- cor(GEBV_test, pheno_test$tbv, use = "complete.obs")
    bias <- coef(lm(pheno_test$tbv ~ GEBV_test))[2]
    rmse <- sqrt(mean((GEBV_test - pheno_test$tbv)^2, na.rm = TRUE))

    predict_cv <- tibble(
      individual_id = pheno_test$individual_id,
      tbv           = pheno_test$tbv,
      GEBV          = GEBV_test,
      scenario      = sc,
      model         = model_name,
      marker        = "mh",
      split         = "test",
      accuracy      = acc,
      bias          = bias,
      rmse          = rmse,
      h2            = h2,
      n_train       = n_train,
      n_test        = nrow(pheno_test)
    )

    # Training accuracy
    GEBV_train <- switch(model_name,
      gblup  = gblup_res$predictions$train,
      bayesR = res$predictions$train,
      bayesA = res$predictions$train
    )
    acc_train  <- cor(GEBV_train, pheno_train$tbv, use = "complete.obs")
    bias_train <- coef(lm(pheno_train$tbv ~ GEBV_train))[2]
    rmse_train <- sqrt(mean((GEBV_train - pheno_train$tbv)^2, na.rm = TRUE))

    saveRDS(predict_cv, file.path(out_dir, model_name, "predict_cv.rds"))
    cat(sprintf("  %s | train: acc=%.4f bias=%.4f rmse=%.4f | test: acc=%.4f bias=%.4f rmse=%.4f | h2=%.3f\n",
                model_name, acc_train, bias_train, rmse_train, acc, bias, rmse, h2))
  }

  cat(sprintf("Scenario %s done in %.1f mins\n", sc,
              as.numeric(difftime(Sys.time(), sc_start, units = "mins"))))
}

cat(sprintf("\nAll scenarios complete | Total: %.1f mins\n",
            as.numeric(difftime(Sys.time(), pipeline_start, units = "mins"))))