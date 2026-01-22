# main.R
# Main orchestration script with parallel processing

suppressPackageStartupMessages({
  library(tidyverse)
  library(parallel)
  library(doParallel)
  library(caret)
})

select <- dplyr::select
filter <- dplyr::filter

config_file <- Sys.getenv("GENOMIC_CONFIG", "config.R")
source(config_file)
setwd(config$base_dir)
source("R-snp/data_prep.R")
source("R-snp/models.R")
source("R-snp/evaluation.R")
source("R-snp/visualization.R")

dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(config$output_dir, "logs"), recursive = TRUE, showWarnings = FALSE)

# Load data
cat("Loading genomic data")
data <- load_genomic_data(config)

cat("Total individuals:", nrow(data$phenotype), "\n")
cat("Total SNP markers:", ncol(data$genotype) - 1, "\n")

# Setup parallel processing
n_cores <- min(config$k_folds, parallel::detectCores() - 1)
cat("Setting up parallel processing with", n_cores, "cores...\n\n")

cl <- makeCluster(n_cores)
registerDoParallel(cl)

clusterEvalQ(cl, {
  library(tidyverse)
  library(sommer)
  library(coda)
  library(gridExtra)
  library(caret)
})

cat("Preparing genomic data...\n")
prep_start <- Sys.time()
prep_result <- local({
  geno_matrix <- data$genotype %>%
    select(-individual_id, -generation, -population_type) %>%
    as.matrix()
  
  cat("Genotype matrix:", nrow(geno_matrix), "individuals ×", ncol(geno_matrix), "SNPs\n")
  
  qc_result <- quality_control(geno_matrix, config)
  geno_processed <- impute_and_standardize(qc_result$filtered)
  dedup <- remove_duplicates(geno_processed$imputed, data$phenotype, config)

  if(dedup$n_removed > 0) {
    geno_processed$imputed <- geno_processed$imputed[dedup$keep_idx, ]
    geno_processed$standardized <- geno_processed$standardized[dedup$keep_idx, ]
    data$phenotype <- data$phenotype[dedup$keep_idx, ]
  }

  GRM_full <- create_grm(geno_processed$imputed)
  
  cat("\nSetting up cross-validation...\n")
  cv_setup <- setup_cv_folds(data$phenotype, config)
  cat("CV method:", cv_setup$method, "| Folds:", config$k_folds, "\n\n")

  list(
    phenotype = data$phenotype,
    geno_std = geno_processed$standardized,
    GRM_full = GRM_full,
    cv_folds = cv_setup$folds,
    qc_n_retained = qc_result$n_retained
  )
})

saveRDS(prep_result, file.path(config$output_dir, "shared_data.rds"))
prep_time <- as.numeric(difftime(Sys.time(), prep_start, units = "hours"))
cat("Data preparation complete | N final =", nrow(prep_result$phenotype), 
    "| SNPs retained =", prep_result$qc_n_retained, 
    sprintf("| Time: %.2f min\n\n", prep_time * 60))

total_cpus <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", parallel::detectCores()))
n_threads_per_fold <- max(1, floor(total_cpus / (n_cores + 2)))
cat("Parallel setup: Threads per fold =", n_threads_per_fold, "\n\n")

clusterExport(cl, c("config", "n_threads_per_fold"))

clusterExport(cl, c("train_gblup", "predict_gblup", "predict_gblup_train",
                    "train_bayesA", "predict_bayesA", "predict_bayesA_train",
                    "train_bayesR", "predict_bayesR", "predict_bayesR_train",
                    "evaluate_model", "evaluate_model_train"))

# Run cross-validation in parallel
all_results <- foreach(fold = 1:config$k_folds,
                       .packages = c("tidyverse", "coda", "sommer", "hibayes"),
                       .errorhandling = "pass", #pass
                       .verbose = TRUE) %dopar% {
  
  source("R-snp/data_prep.R")
  source("R-snp/models.R")
  source("R-snp/evaluation.R")
  source("R-snp/visualization.R")

  select <- dplyr::select
  filter <- dplyr::filter

  tryCatch({
    # Setup logging for this fold
    log_file <- file.path(config$output_dir, "logs", sprintf("fold_%d.log", fold))
    log_con <- file(log_file, open = "wt")
    sink(log_con, type = "output")
    sink(log_con, type = "message")
    
    cat(sprintf("=== FOLD %d START ===\n", fold))
    cat(sprintf("Time: %s\n\n", Sys.time()))

    cat("Loading shared data...\n")
    shared_data <- readRDS(file.path(config$output_dir, "shared_data.rds"))
    cat("Data loaded.\n\n")
    
    cat("Splitting data...\n")
    val_idx <- which(shared_data$cv_folds == fold)
    train_idx <- which(shared_data$cv_folds != fold)
    cat("  Train N:", length(train_idx), "| Test N:", length(val_idx), "\n")
    
    train_ids <- shared_data$phenotype$individual_id[train_idx]
    val_ids <- shared_data$phenotype$individual_id[val_idx]

    cat("Extracting phenotypes...\n")
    y_train <- shared_data$phenotype$phenotype[train_idx]
    y_val <- shared_data$phenotype$phenotype[val_idx]
    tbv_train <- shared_data$phenotype$tbv[train_idx]
    tbv_val <- shared_data$phenotype$tbv[val_idx]
    
    cat("Subsetting genotype matrices...\n")
    geno_train <- shared_data$geno_std[train_idx, ]
    cat("  geno_train dim:", dim(geno_train), "\n")
    geno_val <- shared_data$geno_std[val_idx, ]
    cat("  geno_val dim:", dim(geno_val), "\n")
    geno_train_standardized <- shared_data$geno_std[train_idx, ]
    cat("  geno_train_standardized dim:", dim(geno_train_standardized), "\n")
    geno_val_standardized <- shared_data$geno_std[val_idx, ]
    cat("  geno_val_standardized dim:", dim(geno_val_standardized), "\n")
    
    cat("Subsetting GRM...\n")
    GRM_train <- shared_data$GRM_full[train_idx, train_idx]
    cat("  GRM_train dim:", dim(GRM_train), "\n")
    GRM_val_train <- shared_data$GRM_full[val_idx, train_idx]
    cat("  GRM_val_train dim:", dim(GRM_val_train), "\n")

    rm(shared_data)
    invisible(gc(full = TRUE))
    
    cat("Setup complete, starting models...\n\n")
    
    results_fold <- list()
    predictions_fold <- data.frame(individual_id = val_ids, True_TBV = tbv_val, True_Pheno = y_val)
    
    # GBLUP
    cat("Starting GBLUP...\n")
    gblup_fit <- train_gblup(y_train, GRM_train, n_threads_per_fold)
    cat("  GBLUP fit success:", gblup_fit$success, "\n")
    if(gblup_fit$success) {
      gblup_pred_train <- predict_gblup_train(gblup_fit, GRM_train, y_train)
      gblup_pred <- predict_gblup(gblup_fit, GRM_val_train, y_train, GRM_train)
      
      if(gblup_pred$success) {
        eval_train <- evaluate_model_train(gblup_pred_train$pred, tbv_train, y_train, "GBLUP")
        eval_val <- evaluate_model(gblup_pred$pred, tbv_val, y_val, "GBLUP")
        results_fold$GBLUP <- data.frame(
          Fold = fold,
          eval_val,
          eval_train,
          Train_Time = gblup_fit$train_time,
          Pred_Time = gblup_pred$pred_time,
          h2 = gblup_fit$h2,
          stringsAsFactors = FALSE
        )
        predictions_fold$GBLUP <- gblup_pred$pred

        cat(sprintf("Fold %d - GBLUP | TBV: Train=%.3f Test=%.3f | Pheno: Train=%.3f Test=%.3f | h2=%.3f | Time: %.2f min\n",
              fold, eval_train$Cor_TBV_Train, eval_val$Cor_TBV,
              eval_train$Cor_Pheno_Train, eval_val$Cor_Pheno, gblup_fit$h2,
              gblup_fit$train_time/60))
      }
      rm(gblup_pred_train, gblup_pred, eval_train, eval_val)
    }
    vg_init <- if(gblup_fit$success) gblup_fit$fit$Vu else var(y_train) * 0.5
    ve_init <- if(gblup_fit$success) gblup_fit$fit$Ve else var(y_train) * 0.5
    rm(gblup_fit); invisible(gc(full = TRUE))
    
    # BayesA
    cat("Starting BayesA...\n")
    pheno_data <- data.frame(id = train_ids, phenotype = y_train)
    bayesA_fit <- train_bayesA(pheno_data, geno_train_standardized, train_ids, config, n_threads_per_fold)
    cat("  BayesA fit success:", bayesA_fit$success, "\n")
    if(bayesA_fit$success) {
      bayesA_pred_train <- predict_bayesA_train(bayesA_fit$fit, geno_train_standardized)
      bayesA_pred <- predict_bayesA(bayesA_fit$fit, geno_val_standardized)
      
      if(bayesA_pred$success) {
        eval_train <- evaluate_model_train(bayesA_pred_train$pred, tbv_train, y_train, "BayesA")
        eval_val <- evaluate_model(bayesA_pred$pred, tbv_val, y_val, "BayesA")
        results_fold$BayesA <- data.frame(
          Fold = fold,
          eval_val,
          eval_train,
          Train_Time = bayesA_fit$train_time,
          Pred_Time = bayesA_pred$pred_time,
          h2 = bayesA_fit$h2,
          stringsAsFactors = FALSE
        )
        predictions_fold$BayesA <- bayesA_pred$pred

        cat(sprintf("Fold %d - BayesA | TBV: Train=%.3f Test=%.3f | Pheno: Train=%.3f Test=%.3f | h2=%.3f | Time: %.2f min\n",
              fold, eval_train$Cor_TBV_Train, eval_val$Cor_TBV,
              eval_train$Cor_Pheno_Train, eval_val$Cor_Pheno, bayesA_fit$h2,
              bayesA_fit$train_time/60))
      }
      rm(bayesA_pred_train, bayesA_pred, eval_train, eval_val)
    }
    rm(bayesA_fit); invisible(gc(full = TRUE))
    
    # BayesR
    cat("Starting BayesR...\n")
    bayesR_fit <- train_bayesR(pheno_data, geno_train_standardized, train_ids, vg_init, ve_init, config, n_threads_per_fold)
    cat("  BayesR fit success:", bayesR_fit$success, "\n")
    if(bayesR_fit$success) {
      bayesR_pred_train <- predict_bayesR_train(bayesR_fit$fit, geno_train_standardized)
      bayesR_pred <- predict_bayesR(bayesR_fit$fit, geno_val_standardized)
      
      if(bayesR_pred$success) {
        eval_train <- evaluate_model_train(bayesR_pred_train$pred, tbv_train, y_train, "BayesR")
        eval_val <- evaluate_model(bayesR_pred$pred, tbv_val, y_val, "BayesR")
        results_fold$BayesR <- data.frame(
          Fold = fold,
          eval_val,
          eval_train,
          Train_Time = bayesR_fit$train_time,
          Pred_Time = bayesR_pred$pred_time,
          h2 = bayesR_fit$h2,
          stringsAsFactors = FALSE
        )
        predictions_fold$BayesR <- bayesR_pred$pred

        cat(sprintf("Fold %d - BayesR | TBV: Train=%.3f Test=%.3f | Pheno: Train=%.3f Test=%.3f | h2=%.3f | Time: %.2f min\n",
              fold, eval_train$Cor_TBV_Train, eval_val$Cor_TBV,
              eval_train$Cor_Pheno_Train, eval_val$Cor_Pheno, bayesR_fit$h2,
              bayesR_fit$train_time/60))
      }
      rm(bayesR_pred_train, bayesR_pred, eval_train, eval_val)
    }
    rm(bayesR_fit); invisible(gc(full = TRUE))

    # Cleanup fold-level objects
    rm(geno_train, geno_val, geno_train_standardized, geno_val_standardized,
       GRM_train, GRM_val_train, pheno_data,
       y_train, y_val, tbv_train, tbv_val, train_ids, val_ids)
    invisible(gc(full = TRUE))

    cat(sprintf("\n=== FOLD %d COMPLETE ===\n", fold))
    cat(sprintf("Time: %s\n", Sys.time()))
    sink(type = "message")
    sink(type = "output")
    close(log_con)
    
    result_list <- list(
      results = do.call(rbind, results_fold), 
      predictions = predictions_fold
    )
    rm(results_fold, predictions_fold)
    invisible(gc(full = TRUE))
    
    result_list
  }, error = function(e) {

    # Cleanup on error
      sink(type = "message")
      sink(type = "output")
      if(exists("log_con")) close(log_con)
      
      # Return error info
      list(fold = fold, error = e$message, results = NULL, predictions = NULL)
    })
}

# Stop cluster
stopCluster(cl)

# Check for errors
failed_folds <- sapply(all_results, function(x) !is.null(x$error))
if(any(failed_folds)) {
  cat("\n⚠ WARNING: Some folds failed:\n")
  for(i in which(failed_folds)) {
    cat(sprintf("  Fold %d: %s\n", all_results[[i]]$fold, all_results[[i]]$error))
  }
  cat("\n")
}

# Filter out failed results
valid_results <- all_results[!failed_folds]

if(length(valid_results) == 0) {
  stop("ERROR: All folds failed. Check log files in ", config$output_dir, "/logs/")
}

cat(sprintf("Successfully completed %d out of %d folds\n", 
            length(valid_results), config$k_folds))

# Calculate actual wall-clock CV time (parallel execution)
cv_total_time <- max(sapply(valid_results, function(x) {
  sum(x$results$Train_Time, na.rm = TRUE)
}))

# Also calculate cumulative time (if all were run serially)
cv_cumulative_time <- sum(sapply(valid_results, function(x) {
  sum(x$results$Train_Time, na.rm = TRUE)
}))

cat(sprintf("\n=== Timing Summary ===\n"))
cat(sprintf("Data preparation: %.2f hours\n", prep_time))
cat(sprintf("Actual CV time (parallel): %.2f hours\n", cv_total_time/3600))
cat(sprintf("Cumulative time (serial equivalent): %.2f hours\n", cv_cumulative_time/3600))
cat(sprintf("Parallel speedup: %.1fx\n", cv_cumulative_time/cv_total_time))
cat(sprintf("Total pipeline: %.2f hours\n\n", prep_time + cv_total_time/3600))

# Aggregate results (only from valid folds)
all_cv_results <- do.call(rbind, lapply(valid_results, function(x) x$results))
all_predictions_combined <- do.call(rbind, lapply(valid_results, function(x) x$predictions))

rm(all_results)
gc()

# Evaluate and visualize
aggregated_results <- aggregate_cv_results(all_cv_results)
save_results(all_cv_results, all_predictions_combined, aggregated_results, config)
save_all_plots(all_predictions_combined, all_cv_results, aggregated_results, config)
generate_summary_report(aggregated_results, all_cv_results, 
                        nrow(prep_result$phenotype), 
                        prep_result$qc_n_retained, config)

cat("\n=== Pipeline completed ===\n")