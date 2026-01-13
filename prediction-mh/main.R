# Main Pipeline for Genomic Prediction with Multiple Models

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(sommer)
  library(xgboost)
  library(coda)
  library(gridExtra)
  library(caret)
  library(foreach)
  library(doParallel)
  library(masbayes)
})

select <- dplyr::select
filter <- dplyr::filter

# Load configuration
config_file <- Sys.getenv("GENOMIC_CONFIG", "config.R")
source(config_file)

# Confirm configuration loaded
cat("\n=== Configuration Loaded ===\n")
cat("Config file:", config_file, "\n\n")

cat("Data Sources:\n")
cat("  Data directory:", config$data_dir, "\n")
cat("  Block directory:", config$block_dir, "\n")
cat("  Phenotype file:", file.path(config$data_dir, config$phenotype_file), "\n")
cat("  Pedigree file:", file.path(config$data_dir, config$pedigree_file), "\n")
cat("  N chromosomes:", config$n_chromosomes, "\n\n")

cat("Output directory:", config$output_dir, "\n")
cat("CV folds:", config$k_folds, "\n")
cat("MCMC iterations:", config$mcmc$n_iter, "\n")
cat("============================\n\n")

# Record start time
pipeline_start_time <- Sys.time()
cat("Pipeline started at:", as.character(pipeline_start_time), "\n\n")

setwd(config$base_dir)

# Load functions
source("R/data_preparation.R")
source("R/matrix_construction.R")
source("R/cv_setup.R")
source("R/evaluation.R")
source("R/visualization.R")

# Load models
source("R/models/gblup.R")
source("R/models/bayesR.R")
source("R/models/bayesA.R")
source("R/models/xgboost_model.R")

# Create output directory
dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(config$output_dir, "logs"), recursive = TRUE, showWarnings = FALSE)

# Load data
cat("Loading genomic data...\n")
data <- load_genomic_data(config)

cat("Total individuals:", nrow(data$phenotypes), "\n")
cat("Total haplotype markers:", ncol(data$haplotypes) - 1, "\n\n")

# Setup cross-validation
cat("Setting up cross-validation...\n")
cv_setup <- setup_stratified_cv(data, config)

cat("CV method:", cv_setup$method, "\n")
cat("Number of folds:", config$k_folds, "\n")

if(cv_setup$method == "stratified") {
  cat("Sires per fold:", cv_setup$sires_per_fold, "\n")
  cat("Avg folds per sire:", round(cv_setup$avg_folds_per_sire, 2), "\n")
} else {
  cat("Reason for random CV:", cv_setup$reason, "\n")
}

cat("Fold distribution:\n")
print(table(cv_setup$folds))
cat("\n")

# Setup parallel processing
n_cores <- min(config$k_folds, parallel::detectCores() - 1)
cat("Setting up parallel processing with", n_cores, "cores...\n\n")

cl <- makeCluster(n_cores)
registerDoParallel(cl)

clusterEvalQ(cl, {
  library(tidyverse)
  library(sommer)
  library(xgboost)
  library(coda)
  library(gridExtra)
  library(caret)
  library(masbayes)
})

# Run cross-validation in parallel
all_results <- foreach(fold = 1:config$k_folds,
                       .packages = c("tidyverse", "sommer", "xgboost", "coda", "masbayes"),
                       .errorhandling = "pass",
                       .verbose = FALSE) %dopar% {
  
  # Source all required functions inside worker
  source("R/data_preparation.R")
  source("R/matrix_construction.R")
  source("R/cv_setup.R")
  source("R/evaluation.R")
  source("R/models/gblup.R")
  source("R/models/bayesR.R")
  source("R/models/bayesA.R")
  source("R/models/xgboost_model.R")
  
  select <- dplyr::select
  filter <- dplyr::filter
  
  tryCatch({

    fold_start <- Sys.time()
    
    # Setup fold-specific log file
    log_file <- file.path(config$output_dir, "logs", paste0("fold_", fold, ".log"))
    log_con <- file(log_file, open = "wt")
    sink(log_con, type = "output")
    sink(log_con, type = "message")
    
    cat(sprintf("\n========== FOLD %d/%d - START ==========\n", fold, config$k_folds))
    cat("Started at:", as.character(Sys.time()), "\n\n")
    
    # Get fold info
    fold_info <- get_fold_split_info(data, cv_setup, fold)
    cat(sprintf("Train N=%d (%.1f%%), Test N=%d (%.1f%%)\n",
                fold_info$n_train, fold_info$prop_train * 100,
                fold_info$n_test, fold_info$prop_test * 100))
    
    # Split data
    split <- split_fold_data(data, cv_setup$folds, fold)
    
    # Compute allele frequencies ONCE per fold
    hap_train <- split$train$haplotypes %>% select(-individual_id)
    allele_freq <- compute_allele_frequencies(hap_train, config)

    # Construct matrices
    cat("\nConstructing W matrices and GRM...\n")
    matrices <- construct_matrices(split$train, split$test, allele_freq, config)

    cat("\n=== Matrix Construction Summary ===\n")
    cat("Genotype preview (first 5 ind, first 8 blocks):\n")
    print(head(split$train$haplotypes[, 1:min(9, ncol(split$train$haplotypes))], 5))
    cat("\nW_ah preview (first 5 ind, first 10 alleles):\n")
    print(head(matrices$W_train[, 1:min(10, ncol(matrices$W_train))], 5))
    cat("\n")
    
    cat("W_train dimension:", dim(matrices$W_train), "\n")
    cat("Columns with zero variance:", sum(apply(matrices$W_train, 2, var) < 1e-10), "\n")
    cat("NA value numbers:", sum(is.na(matrices$W_train)), "\n")
    cat("W_test dimension:", dim(matrices$W_test), "\n")
    cat("Columns with zero variance:", sum(apply(matrices$W_test, 2, var) < 1e-10), "\n")
    cat("NA value numbers:", sum(is.na(matrices$W_test)), "\n")
    cat("Total alleles:", ncol(matrices$W_train), "\n")
    cat("k_αh =", round(matrices$k_ah, 4), "\n\n")
    
    # Run GBLUP
    cat("Running GBLUP...\n")
    gblup_start <- Sys.time()
    gblup_res <- run_gblup(matrices, split, config)
    gblup_time <- as.numeric(difftime(Sys.time(), gblup_start, units = "hours"))
    gblup_res$runtime <- gblup_time
    cat("GBLUP Training Accuracy (Phenotype):", round(gblup_res$accuracy$train_pheno, 4), "\n")
    cat("GBLUP Training Accuracy (TBV):", round(gblup_res$accuracy$train_tbv, 4), "\n")
    cat("GBLUP Test Accuracy (Phenotype):", round(gblup_res$accuracy$test_pheno, 4), "\n")
    cat("GBLUP Test Accuracy (TBV):", round(gblup_res$accuracy$test_tbv, 4), "\n")
    cat("h² =", round(gblup_res$varcomp$h2, 4), "\n\n")
    cat("Runtime:", sprintf("%.1f hours", gblup_time), "\n\n")
    
    # Run BayesR
    cat("Running BayesR...\n")
    bayesR_start <- Sys.time()
    bayesR_res <- run_bayesR(matrices, split, gblup_res$varcomp, config, fold = fold)
    bayesR_time <- as.numeric(difftime(Sys.time(), bayesR_start, units = "hours"))
    bayesR_res$runtime <- bayesR_time
    cat("BayesR Training Accuracy (Phenotype):", round(bayesR_res$accuracy$train_pheno, 4), "\n")
    cat("BayesR Training Accuracy (TBV):", round(bayesR_res$accuracy$train_tbv, 4), "\n")
    cat("BayesR Test Accuracy (Phenotype):", round(bayesR_res$accuracy$test_pheno, 4), "\n")
    cat("BayesR Test Accuracy (TBV):", round(bayesR_res$accuracy$test_tbv, 4), "\n")
    cat("h² =", round(bayesR_res$varcomp$h2, 4), "\n")
    cat("ESS:", round(bayesR_res$diagnostics$ess), 
        "| Geweke Z:", round(bayesR_res$diagnostics$geweke_z, 3), "\n\n")
    cat("Runtime:", sprintf("%.1f hours", bayesR_time), "\n\n")
    
    # Run BayesA
    cat("Running BayesA...\n")
    bayesA_start <- Sys.time()
    bayesA_res <- run_bayesA(matrices, split, gblup_res$varcomp, config, fold = fold)
    bayesA_time <- as.numeric(difftime(Sys.time(), bayesA_start, units = "hours"))
    bayesA_res$runtime <- bayesA_time
    cat("BayesA Training Accuracy (Phenotype):", round(bayesA_res$accuracy$train_pheno, 4), "\n")
    cat("BayesA Training Accuracy (TBV):", round(bayesA_res$accuracy$train_tbv, 4), "\n")
    cat("BayesA Test Accuracy (Phenotype):", round(bayesA_res$accuracy$test_pheno, 4), "\n")
    cat("BayesA Test Accuracy (TBV):", round(bayesA_res$accuracy$test_tbv, 4), "\n")
    cat("h² =", round(bayesA_res$varcomp$h2, 4), "\n")
    cat("ESS:", round(bayesA_res$diagnostics$ess),
        "| Geweke Z:", round(bayesA_res$diagnostics$geweke_z, 3), "\n\n")
    cat("Runtime:", sprintf("%.1f hours", bayesA_time), "\n\n")
    
    # Run XGBoost
    cat("Running XGBoost...\n")
    xgb_start <- Sys.time()
    xgb_res <- run_xgboost(matrices, split, config)
    xgb_time <- as.numeric(difftime(Sys.time(), xgb_start, units = "hours"))
    xgb_res$runtime <- xgb_time
    cat("XGBoost Training Accuracy (Phenotype):", round(xgb_res$accuracy$train_pheno, 4), "\n")
    cat("XGBoost Training Accuracy (TBV):", round(xgb_res$accuracy$train_tbv, 4), "\n")
    cat("XGBoost Test Accuracy (Phenotype):", round(xgb_res$accuracy$test_pheno, 4), "\n")
    cat("XGBoost Test Accuracy (TBV):", round(xgb_res$accuracy$test_tbv, 4), "\n\n")
    cat("Runtime:", sprintf("%.1f hours", xgb_time), "\n\n")
    
    # Save fold results
    cat("Saving fold results...\n")
    save_fold_results(fold, gblup_res, bayesR_res, bayesA_res,
                     xgb_res, matrices, config$output_dir)

    save_gebv_with_regression(fold, gblup_res, bayesR_res, bayesA_res,
                          xgb_res, split, config$output_dir)
    
    # Calculate fold runtime
    fold_end <- Sys.time()
    fold_runtime <- as.numeric(difftime(fold_end, fold_start, units = "hours"))

    # Model runtime summary
    cat("\n========== MODEL RUNTIME SUMMARY ==========\n")
    cat("GBLUP:  ", sprintf("%6.1f hours", gblup_time), "\n")
    cat("BayesR: ", sprintf("%6.1f hours", bayesR_time), "\n")
    cat("BayesA: ", sprintf("%6.1f hours", bayesA_time), "\n")
    cat("XGBoost:", sprintf("%6.1f hours", xgb_time), "\n")
    cat("===========================================\n\n")
    
    cat(sprintf("========== FOLD %d - DONE ==========\n", fold))
    cat("Fold runtime:", sprintf("%.1f hours", fold_runtime), "\n")
    cat("Completed at:", as.character(fold_end), "\n")
    
    # Close log file
    sink(type = "message")
    sink(type = "output")
    close(log_con)
    
    gc()
    
    # Return aggregated results
    result <- aggregate_fold_results(gblup_res, bayesR_res, bayesA_res,
                                    xgb_res, fold)
    
    result$runtime <- list(
      fold_total = fold_runtime,
      gblup = gblup_time,
      bayesR = bayesR_time,
      bayesA = bayesA_time,
      xgboost = xgb_time
    )

    result
    
  }, error = function(e) {
    # Close log on error
    sink(type = "message")
    sink(type = "output")
    if(exists("log_con")) close(log_con)
    
    list(fold = fold, error = e$message)
  })
}

# Stop cluster
stopCluster(cl)

# Calculate parallel processing time
parallel_end_time <- Sys.time()
parallel_runtime_hours <- as.numeric(difftime(parallel_end_time, pipeline_start_time, units = "hours"))
parallel_runtime_mins <- as.numeric(difftime(parallel_end_time, pipeline_start_time, units = "mins"))

cat("\n=== Parallel Processing Complete ===\n")
cat("Total folds processed:", length(all_results), "/", config$k_folds, "\n")
cat("Parallel runtime:", sprintf("%.2f hours (%.1f mins)", parallel_runtime_hours, parallel_runtime_mins), "\n")

# Check for failures
failed_folds <- sapply(all_results, function(x) !is.null(x$error))
if(any(failed_folds)) {
  cat("\n⚠ Failed folds:", which(failed_folds), "\n")
  for(i in which(failed_folds)) {
    cat("  Fold", i, "error:", all_results[[i]]$error, "\n")
  }
  all_results <- all_results[!failed_folds]
}

# Calculate average runtime per model across folds
if(length(all_results) > 0) {
  runtime_summary <- data.frame(
    fold = sapply(all_results, function(x) x$fold),
    gblup = sapply(all_results, function(x) x$runtime$gblup),
    bayesR = sapply(all_results, function(x) x$runtime$bayesR),
    bayesA = sapply(all_results, function(x) x$runtime$bayesA),
    xgboost = sapply(all_results, function(x) x$runtime$xgboost),
    fold_total = sapply(all_results, function(x) x$runtime$fold_total)
  )

  cat("\n=== AVERAGE RUNTIME PER MODEL ===\n")
  cat("GBLUP:  ", sprintf("%.1f ± %.1f hours", mean(runtime_summary$gblup), sd(runtime_summary$gblup)), "\n")
  cat("BayesR: ", sprintf("%.1f ± %.1f hours", mean(runtime_summary$bayesR), sd(runtime_summary$bayesR)), "\n")
  cat("BayesA: ", sprintf("%.1f ± %.1f hours", mean(runtime_summary$bayesA), sd(runtime_summary$bayesA)), "\n")
  cat("XGBoost:", sprintf("%.1f ± %.1f hours", mean(runtime_summary$xgboost), sd(runtime_summary$xgboost)), "\n")
  cat("\nAverage fold runtime:", sprintf("%.1f ± %.1f hours", mean(runtime_summary$fold_total), sd(runtime_summary$fold_total)), "\n")
  cat("=================================\n")
  
  # Save runtime summary
  write.csv(runtime_summary, 
            file.path(config$output_dir, "runtime_summary.csv"),
            row.names = FALSE)
}

cat("\nFold-specific logs saved in:", file.path(config$output_dir, "logs"), "\n")
cat("====================================\n\n")

cat("Completed", length(all_results), "/", config$k_folds, "folds\n\n")

# Summarize results
cat("Summarizing cross-validation results...\n")
cv_summary <- save_cv_summary(all_results, config$output_dir)

# Print summary
print_cv_summary(cv_summary$summary)

# Create visualizations
cat("\nCreating visualizations...\n")
plot_accuracy_comparison(cv_summary$results, config$output_dir)
plot_heritability_comparison(cv_summary$results, config$output_dir)
plot_fold_variability(cv_summary$results, config$output_dir)
plot_variance_components(cv_summary$results, config$output_dir)
plot_train_vs_test_accuracy(cv_summary$results, config$output_dir)

cat("\nAll outputs saved to:", config$output_dir, "\n")
# Calculate total runtime
pipeline_end_time <- Sys.time()
total_runtime_hours <- as.numeric(difftime(pipeline_end_time, pipeline_start_time, units = "hours"))
total_runtime_mins <- as.numeric(difftime(pipeline_end_time, pipeline_start_time, units = "mins"))

cat("\n=== Pipeline Complete ===\n")
cat("Total runtime:", sprintf("%.2f hours (%.1f mins)", total_runtime_hours, total_runtime_mins), "\n")
cat("Started:", as.character(pipeline_start_time), "\n")
cat("Ended:", as.character(pipeline_end_time), "\n")