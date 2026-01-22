# Model Evaluation and Results Aggregation Functions

calculate_regression_slope <- function(predicted, observed) {
  fit <- lm(observed ~ predicted)
  coef(fit)[2]  # Return slope (b)
}

aggregate_fold_results <- function(gblup_res, bayesR_res, bayesA_res, fold_id) {
  
  list(
    fold = fold_id,
    gblup_train_pheno = gblup_res$accuracy$train_pheno,
    gblup_train_tbv = gblup_res$accuracy$train_tbv,
    gblup_test_pheno = gblup_res$accuracy$test_pheno,
    gblup_test_tbv = gblup_res$accuracy$test_tbv,
    gblup_h2 = gblup_res$varcomp$h2,
    gblup_mean_diag_train = gblup_res$varcomp$mean_diag_train,
    gblup_mean_diag_combined = gblup_res$varcomp$mean_diag_combined,
    bayesR_train_pheno = bayesR_res$accuracy$train_pheno,
    bayesR_train_tbv = bayesR_res$accuracy$train_tbv,
    bayesR_test_pheno = bayesR_res$accuracy$test_pheno,
    bayesR_test_tbv = bayesR_res$accuracy$test_tbv,
    bayesR_h2 = bayesR_res$varcomp$h2,
    bayesR_mean_ess = bayesR_res$diagnostics$mean_ess,
    bayesR_min_ess = bayesR_res$diagnostics$min_ess,
    bayesR_ess_sigma2_e = if(is.list(bayesR_res$diagnostics$ess)) bayesR_res$diagnostics$ess$sigma2_e else bayesR_res$diagnostics$ess,
    bayesR_geweke_sigma2_e = if(is.list(bayesR_res$diagnostics$geweke_z)) bayesR_res$diagnostics$geweke_z$sigma2_e else bayesR_res$diagnostics$geweke_z,
    bayesA_train_pheno = bayesA_res$accuracy$train_pheno,
    bayesA_train_tbv = bayesA_res$accuracy$train_tbv,
    bayesA_test_pheno = bayesA_res$accuracy$test_pheno,
    bayesA_test_tbv = bayesA_res$accuracy$test_tbv,
    bayesA_h2 = bayesA_res$varcomp$h2,
    gblup_b_train_pheno = gblup_res$regression_slopes$train_pheno,
    gblup_b_train_tbv = gblup_res$regression_slopes$train_tbv,
    gblup_b_test_pheno = gblup_res$regression_slopes$test_pheno,
    gblup_b_test_tbv = gblup_res$regression_slopes$test_tbv,
    bayesR_b_train_pheno = bayesR_res$regression_slopes$train_pheno,
    bayesR_b_train_tbv = bayesR_res$regression_slopes$train_tbv,
    bayesR_b_test_pheno = bayesR_res$regression_slopes$test_pheno,
    bayesR_b_test_tbv = bayesR_res$regression_slopes$test_tbv,
    bayesA_b_train_pheno = bayesA_res$regression_slopes$train_pheno,
    bayesA_b_train_tbv = bayesA_res$regression_slopes$train_tbv,
    bayesA_b_test_pheno = bayesA_res$regression_slopes$test_pheno,
    bayesA_b_test_tbv = bayesA_res$regression_slopes$test_tbv,
    bayesA_mean_ess = bayesA_res$diagnostics$mean_ess,
    bayesA_min_ess = bayesA_res$diagnostics$min_ess,
    bayesA_ess_sigma2_e = if(is.list(bayesA_res$diagnostics$ess)) bayesA_res$diagnostics$ess$sigma2_e else bayesA_res$diagnostics$ess,
    bayesA_geweke_sigma2_e = if(is.list(bayesA_res$diagnostics$geweke_z)) bayesA_res$diagnostics$geweke_z$sigma2_e else bayesA_res$diagnostics$geweke_z,
    gblup_rmse_train_pheno = gblup_res$rmse$train_pheno,
    gblup_rmse_train_tbv = gblup_res$rmse$train_tbv,
    gblup_rmse_test_pheno = gblup_res$rmse$test_pheno,
    gblup_rmse_test_tbv = gblup_res$rmse$test_tbv,
    bayesR_rmse_train_pheno = bayesR_res$rmse$train_pheno,
    bayesR_rmse_train_tbv = bayesR_res$rmse$train_tbv,
    bayesR_rmse_test_pheno = bayesR_res$rmse$test_pheno,
    bayesR_rmse_test_tbv = bayesR_res$rmse$test_tbv,
    bayesA_rmse_train_pheno = bayesA_res$rmse$train_pheno,
    bayesA_rmse_train_tbv = bayesA_res$rmse$train_tbv,
    bayesA_rmse_test_pheno = bayesA_res$rmse$test_pheno,
    bayesA_rmse_test_tbv = bayesA_res$rmse$test_tbv,
    gblup_r2_train_pheno = gblup_res$r2$train_pheno,
    gblup_r2_train_tbv = gblup_res$r2$train_tbv,
    gblup_r2_test_pheno = gblup_res$r2$test_pheno,
    gblup_r2_test_tbv = gblup_res$r2$test_tbv,
    bayesR_r2_train_pheno = bayesR_res$r2$train_pheno,
    bayesR_r2_train_tbv = bayesR_res$r2$train_tbv,
    bayesR_r2_test_pheno = bayesR_res$r2$test_pheno,
    bayesR_r2_test_tbv = bayesR_res$r2$test_tbv,
    bayesA_r2_train_pheno = bayesA_res$r2$train_pheno,
    bayesA_r2_train_tbv = bayesA_res$r2$train_tbv,
    bayesA_r2_test_pheno = bayesA_res$r2$test_pheno,
    bayesA_r2_test_tbv = bayesA_res$r2$test_tbv,
    gblup_n_cores = gblup_res$n_cores,
    bayesR_n_cores = bayesR_res$n_cores,
    bayesA_n_cores = bayesA_res$n_cores,
    gblup_core_hours = gblup_res$runtime * gblup_res$n_cores,
    bayesR_core_hours = bayesR_res$runtime * bayesR_res$n_cores,
    bayesA_core_hours = bayesA_res$runtime * bayesA_res$n_cores,
    gblup_efficiency = (gblup_res$runtime * gblup_res$n_cores) / gblup_res$accuracy$test_tbv,
    bayesR_efficiency = (bayesR_res$runtime * bayesR_res$n_cores) / bayesR_res$accuracy$test_tbv,
    bayesA_efficiency = (bayesA_res$runtime * bayesA_res$n_cores) / bayesA_res$accuracy$test_tbv
  )
}

summarize_cv_results <- function(results_list) {
  
  results_df <- do.call(rbind, lapply(results_list, as.data.frame))
  
  summary_stats <- results_df %>%
    summarise(
      GBLUP_train_pheno_mean = mean(gblup_train_pheno),
      GBLUP_train_pheno_sd = sd(gblup_train_pheno),
      GBLUP_train_tbv_mean = mean(gblup_train_tbv),
      GBLUP_train_tbv_sd = sd(gblup_train_tbv),
      GBLUP_test_pheno_mean = mean(gblup_test_pheno),
      GBLUP_test_pheno_sd = sd(gblup_test_pheno),
      GBLUP_test_tbv_mean = mean(gblup_test_tbv),
      GBLUP_test_tbv_sd = sd(gblup_test_tbv),
      GBLUP_h2_mean = mean(gblup_h2),
      GBLUP_h2_sd = sd(gblup_h2),
      GBLUP_mean_diag_train_mean = mean(gblup_mean_diag_train),
      GBLUP_mean_diag_train_sd = sd(gblup_mean_diag_train),
      GBLUP_mean_diag_combined_mean = mean(gblup_mean_diag_combined),
      GBLUP_mean_diag_combined_sd = sd(gblup_mean_diag_combined),
      BayesR_train_pheno_mean = mean(bayesR_train_pheno),
      BayesR_train_pheno_sd = sd(bayesR_train_pheno),
      BayesR_train_tbv_mean = mean(bayesR_train_tbv),
      BayesR_train_tbv_sd = sd(bayesR_train_tbv),
      BayesR_test_pheno_mean = mean(bayesR_test_pheno),
      BayesR_test_pheno_sd = sd(bayesR_test_pheno),
      BayesR_test_tbv_mean = mean(bayesR_test_tbv),
      BayesR_test_tbv_sd = sd(bayesR_test_tbv),
      BayesR_h2_mean = mean(bayesR_h2),
      BayesR_h2_sd = sd(bayesR_h2),
      BayesR_mean_ess_mean = mean(bayesR_mean_ess, na.rm = TRUE),
      BayesR_mean_ess_sd = sd(bayesR_mean_ess, na.rm = TRUE),
      BayesR_min_ess_mean = mean(bayesR_min_ess, na.rm = TRUE),
      BayesR_min_ess_sd = sd(bayesR_min_ess, na.rm = TRUE),
      BayesR_ess_sigma2_e_mean = mean(bayesR_ess_sigma2_e, na.rm = TRUE),
      BayesR_geweke_sigma2_e_mean = mean(abs(bayesR_geweke_sigma2_e), na.rm = TRUE),
      BayesA_train_pheno_mean = mean(bayesA_train_pheno),
      BayesA_train_pheno_sd = sd(bayesA_train_pheno),
      BayesA_train_tbv_mean = mean(bayesA_train_tbv),
      BayesA_train_tbv_sd = sd(bayesA_train_tbv),
      BayesA_test_pheno_mean = mean(bayesA_test_pheno),
      BayesA_test_pheno_sd = sd(bayesA_test_pheno),
      BayesA_test_tbv_mean = mean(bayesA_test_tbv),
      BayesA_test_tbv_sd = sd(bayesA_test_tbv),
      BayesA_h2_mean = mean(bayesA_h2),
      BayesA_h2_sd = sd(bayesA_h2),
      BayesA_mean_ess_mean = mean(bayesA_mean_ess, na.rm = TRUE),
      BayesA_mean_ess_sd = sd(bayesA_mean_ess, na.rm = TRUE),
      BayesA_min_ess_mean = mean(bayesA_min_ess, na.rm = TRUE),
      BayesA_min_ess_sd = sd(bayesA_min_ess, na.rm = TRUE),
      BayesA_ess_sigma2_e_mean = mean(bayesA_ess_sigma2_e, na.rm = TRUE),
      BayesA_geweke_sigma2_e_mean = mean(abs(bayesA_geweke_sigma2_e), na.rm = TRUE),
      GBLUP_b_train_pheno_mean = mean(gblup_b_train_pheno),
      GBLUP_b_train_pheno_sd = sd(gblup_b_train_pheno),
      GBLUP_b_train_tbv_mean = mean(gblup_b_train_tbv),
      GBLUP_b_train_tbv_sd = sd(gblup_b_train_tbv),
      GBLUP_b_test_pheno_mean = mean(gblup_b_test_pheno),
      GBLUP_b_test_pheno_sd = sd(gblup_b_test_pheno),
      GBLUP_b_test_tbv_mean = mean(gblup_b_test_tbv),
      GBLUP_b_test_tbv_sd = sd(gblup_b_test_tbv),
      BayesR_b_train_pheno_mean = mean(bayesR_b_train_pheno),
      BayesR_b_train_pheno_sd = sd(bayesR_b_train_pheno),
      BayesR_b_train_tbv_mean = mean(bayesR_b_train_tbv),
      BayesR_b_train_tbv_sd = sd(bayesR_b_train_tbv),
      BayesR_b_test_pheno_mean = mean(bayesR_b_test_pheno),
      BayesR_b_test_pheno_sd = sd(bayesR_b_test_pheno),
      BayesR_b_test_tbv_mean = mean(bayesR_b_test_tbv),
      BayesR_b_test_tbv_sd = sd(bayesR_b_test_tbv),
      BayesA_b_train_pheno_mean = mean(bayesA_b_train_pheno),
      BayesA_b_train_pheno_sd = sd(bayesA_b_train_pheno),
      BayesA_b_train_tbv_mean = mean(bayesA_b_train_tbv),
      BayesA_b_train_tbv_sd = sd(bayesA_b_train_tbv),
      BayesA_b_test_pheno_mean = mean(bayesA_b_test_pheno),
      BayesA_b_test_pheno_sd = sd(bayesA_b_test_pheno),
      BayesA_b_test_tbv_mean = mean(bayesA_b_test_tbv),
      BayesA_b_test_tbv_sd = sd(bayesA_b_test_tbv),
      GBLUP_rmse_train_pheno_mean = mean(gblup_rmse_train_pheno),
      GBLUP_rmse_train_pheno_sd = sd(gblup_rmse_train_pheno),
      GBLUP_rmse_test_pheno_mean = mean(gblup_rmse_test_pheno),
      GBLUP_rmse_test_pheno_sd = sd(gblup_rmse_test_pheno),
      GBLUP_rmse_train_tbv_mean = mean(gblup_rmse_train_tbv),
      GBLUP_rmse_train_tbv_sd = sd(gblup_rmse_train_tbv),
      GBLUP_rmse_test_tbv_mean = mean(gblup_rmse_test_tbv),
      GBLUP_rmse_test_tbv_sd = sd(gblup_rmse_test_tbv),
      BayesR_rmse_train_pheno_mean = mean(bayesR_rmse_train_pheno),
      BayesR_rmse_train_pheno_sd = sd(bayesR_rmse_train_pheno),
      BayesR_rmse_test_pheno_mean = mean(bayesR_rmse_test_pheno),
      BayesR_rmse_test_pheno_sd = sd(bayesR_rmse_test_pheno),
      BayesR_rmse_train_tbv_mean = mean(bayesR_rmse_train_tbv),
      BayesR_rmse_train_tbv_sd = sd(bayesR_rmse_train_tbv),
      BayesR_rmse_test_tbv_mean = mean(bayesR_rmse_test_tbv),
      BayesR_rmse_test_tbv_sd = sd(bayesR_rmse_test_tbv),
      BayesA_rmse_train_pheno_mean = mean(bayesA_rmse_train_pheno),
      BayesA_rmse_train_pheno_sd = sd(bayesA_rmse_train_pheno),
      BayesA_rmse_test_pheno_mean = mean(bayesA_rmse_test_pheno),
      BayesA_rmse_test_pheno_sd = sd(bayesA_rmse_test_pheno),
      BayesA_rmse_train_tbv_mean = mean(bayesA_rmse_train_tbv),
      BayesA_rmse_train_tbv_sd = sd(bayesA_rmse_train_tbv),
      BayesA_rmse_test_tbv_mean = mean(bayesA_rmse_test_tbv),
      BayesA_rmse_test_tbv_sd = sd(bayesA_rmse_test_tbv),
      GBLUP_r2_train_pheno_mean = mean(gblup_r2_train_pheno),
      GBLUP_r2_train_pheno_sd = sd(gblup_r2_train_pheno),
      GBLUP_r2_test_pheno_mean = mean(gblup_r2_test_pheno),
      GBLUP_r2_test_pheno_sd = sd(gblup_r2_test_pheno),
      GBLUP_r2_train_tbv_mean = mean(gblup_r2_train_tbv),
      GBLUP_r2_train_tbv_sd = sd(gblup_r2_train_tbv),
      GBLUP_r2_test_tbv_mean = mean(gblup_r2_test_tbv),
      GBLUP_r2_test_tbv_sd = sd(gblup_r2_test_tbv),
      BayesR_r2_train_pheno_mean = mean(bayesR_r2_train_pheno),
      BayesR_r2_train_pheno_sd = sd(bayesR_r2_train_pheno),
      BayesR_r2_test_pheno_mean = mean(bayesR_r2_test_pheno),
      BayesR_r2_test_pheno_sd = sd(bayesR_r2_test_pheno),
      BayesR_r2_train_tbv_mean = mean(bayesR_r2_train_tbv),
      BayesR_r2_train_tbv_sd = sd(bayesR_r2_train_tbv),
      BayesR_r2_test_tbv_mean = mean(bayesR_r2_test_tbv),
      BayesR_r2_test_tbv_sd = sd(bayesR_r2_test_tbv),
      BayesA_r2_train_pheno_mean = mean(bayesA_r2_train_pheno),
      BayesA_r2_train_pheno_sd = sd(bayesA_r2_train_pheno),
      BayesA_r2_test_pheno_mean = mean(bayesA_r2_test_pheno),
      BayesA_r2_test_pheno_sd = sd(bayesA_r2_test_pheno),
      BayesA_r2_train_tbv_mean = mean(bayesA_r2_train_tbv),
      BayesA_r2_train_tbv_sd = sd(bayesA_r2_train_tbv),
      BayesA_r2_test_tbv_mean = mean(bayesA_r2_test_tbv),
      BayesA_r2_test_tbv_sd = sd(bayesA_r2_test_tbv),
      GBLUP_core_hours_mean = mean(gblup_core_hours),
      GBLUP_core_hours_sd = sd(gblup_core_hours),
      BayesR_core_hours_mean = mean(bayesR_core_hours),
      BayesR_core_hours_sd = sd(bayesR_core_hours),
      BayesA_core_hours_mean = mean(bayesA_core_hours),
      BayesA_core_hours_sd = sd(bayesA_core_hours),
      GBLUP_efficiency_mean = mean(gblup_efficiency),
      GBLUP_efficiency_sd = sd(gblup_efficiency),
      BayesR_efficiency_mean = mean(bayesR_efficiency),
      BayesR_efficiency_sd = sd(bayesR_efficiency),
      BayesA_efficiency_mean = mean(bayesA_efficiency),
      BayesA_efficiency_sd = sd(bayesA_efficiency)
    )
  
  list(
    results = results_df,
    summary = summary_stats
  )
}

save_gebv_with_regression <- function(fold_id, gblup_res, bayesR_res, 
                                      bayesA_res, split, output_dir) {
  
  fold_dir <- file.path(output_dir, paste0("fold_", fold_id))
  
  # Training set
  train_gebv <- split$train$phenotypes %>%
    select(individual_id, tbv, phenotype) %>%
    mutate(
      GEBV_GBLUP = gblup_res$predictions$train,
      GEBV_BayesR = bayesR_res$predictions$train,
      GEBV_BayesA = bayesA_res$predictions$train,
      b_GBLUP_tbv = gblup_res$regression_slopes$train_tbv,
      b_BayesR_tbv = bayesR_res$regression_slopes$train_tbv,
      b_BayesA_tbv = bayesA_res$regression_slopes$train_tbv
    )
  
  # Test set
  test_gebv <- split$test$phenotypes %>%
    select(individual_id, tbv, phenotype) %>%
    mutate(
      GEBV_GBLUP = gblup_res$predictions$test,
      GEBV_BayesR = bayesR_res$predictions$test,
      GEBV_BayesA = bayesA_res$predictions$test,
      b_GBLUP_tbv = gblup_res$regression_slopes$test_tbv,
      b_BayesR_tbv = bayesR_res$regression_slopes$test_tbv,
      b_BayesA_tbv = bayesA_res$regression_slopes$test_tbv
    )
  
  write.csv(train_gebv, file.path(fold_dir, "train_gebv.csv"), row.names = FALSE)
  write.csv(test_gebv, file.path(fold_dir, "test_gebv.csv"), row.names = FALSE)
}

save_fold_results <- function(fold_id, gblup_res, bayesR_res, bayesA_res, matrices, output_dir) {
  
  fold_dir <- file.path(output_dir, paste0("fold_", fold_id))
  dir.create(fold_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save aggregated results
  fold_results <- aggregate_fold_results(
    gblup_res, bayesR_res, bayesA_res, fold_id
  )
  saveRDS(fold_results, file.path(fold_dir, "results.rds"))
  
  # Save BayesR markers
  if(!is.null(bayesR_res$markers)) {
    write.csv(
      bayesR_res$markers,
      file.path(fold_dir, "bayesR_all_markers.csv"),
      row.names = FALSE
    )
    
    significant <- bayesR_res$markers %>% filter(significant == TRUE)
    if(nrow(significant) > 0) {
      write.csv(
        significant,
        file.path(fold_dir, "bayesR_significant_markers.csv"),
        row.names = FALSE
      )
    }
    
    high_pip <- bayesR_res$markers %>% filter(PIP > 0.5)
    if(nrow(high_pip) > 0) {
      write.csv(
        high_pip,
        file.path(fold_dir, "bayesR_high_pip_markers.csv"),
        row.names = FALSE
      )
    }
  }
  
  # Save BayesA markers
  if(!is.null(bayesA_res$markers)) {
    write.csv(
      bayesA_res$markers,
      file.path(fold_dir, "bayesA_all_markers.csv"),
      row.names = FALSE
    )
    
    significant_A <- bayesA_res$markers %>% filter(significant == TRUE)
    if(nrow(significant_A) > 0) {
      write.csv(
        significant_A,
        file.path(fold_dir, "bayesA_significant_markers.csv"),
        row.names = FALSE
      )
    }
  }
  
  # Save baseline alleles
  if(!is.null(matrices$dropped_alleles) && nrow(matrices$dropped_alleles) > 0) {
    write.csv(
      matrices$dropped_alleles,
      file.path(fold_dir, "baseline_alleles.csv"),
      row.names = FALSE
    )
  }
  
  # Save convergence diagnostics
  diagnostics_summary <- data.frame(
    model = character(0),
    h2 = numeric(0),
    ess = numeric(0),
    geweke_z = numeric(0),
    stringsAsFactors = FALSE
  )

  # GBLUP
  diagnostics_summary <- rbind(diagnostics_summary, data.frame(
    model = "GBLUP",
    h2 = gblup_res$varcomp$h2,
    ess = NA_real_,
    geweke_z = NA_real_
  ))

  # BayesR
  diagnostics_summary <- rbind(diagnostics_summary, data.frame(
    model = "BayesR",
    h2 = bayesR_res$varcomp$h2,
    ess = bayesR_res$diagnostics$ess,
    geweke_z = bayesR_res$diagnostics$geweke_z
  ))

  # BayesA
  diagnostics_summary <- rbind(diagnostics_summary, data.frame(
    model = "BayesA",
    h2 = bayesA_res$varcomp$h2,
    ess = bayesA_res$diagnostics$ess,
    geweke_z = bayesA_res$diagnostics$geweke_z
  ))

  write.csv(
    diagnostics_summary,
    file.path(fold_dir, "diagnostics_summary.csv"),
    row.names = FALSE
  )
  
  invisible(fold_results)
}

save_cv_summary <- function(cv_results, output_dir) {
  
  summary <- summarize_cv_results(cv_results)
  
  write.csv(
    summary$results,
    file.path(output_dir, "all_folds_results.csv"),
    row.names = FALSE
  )
  
  write.csv(
    summary$summary,
    file.path(output_dir, "summary_statistics.csv"),
    row.names = FALSE
  )
  
  summary
}

print_cv_summary <- function(summary_stats) {
  
  cat("\n=== Cross-Validation Summary ===\n")
  
  cat("\nTEST Accuracy ± SD (TBV):\n")
  cat("GBLUP:  ", sprintf("%.4f ± %.4f", 
                         summary_stats$GBLUP_test_tbv_mean, 
                         summary_stats$GBLUP_test_tbv_sd), "\n")
  cat("BayesR: ", sprintf("%.4f ± %.4f", 
                         summary_stats$BayesR_test_tbv_mean, 
                         summary_stats$BayesR_test_tbv_sd), "\n")
  cat("BayesA: ", sprintf("%.4f ± %.4f", 
                         summary_stats$BayesA_test_tbv_mean, 
                         summary_stats$BayesA_test_tbv_sd), "\n")
  
  cat("\nMean h² ± SD:\n")
  cat("GBLUP:  ", sprintf("%.4f ± %.4f", 
                         summary_stats$GBLUP_h2_mean, 
                         summary_stats$GBLUP_h2_sd), "\n")
  cat("BayesR: ", sprintf("%.4f ± %.4f", 
                         summary_stats$BayesR_h2_mean, 
                         summary_stats$BayesR_h2_sd), "\n")
  cat("BayesA: ", sprintf("%.4f ± %.4f", 
                         summary_stats$BayesA_h2_mean, 
                         summary_stats$BayesA_h2_sd), "\n")
  
  cat("\nTRAINING RMSE ± SD (TBV):\n")
  cat("GBLUP:  ", sprintf("%.4f ± %.4f", 
                         summary_stats$GBLUP_rmse_train_tbv_mean, 
                         summary_stats$GBLUP_rmse_train_tbv_sd), "\n")
  cat("BayesR: ", sprintf("%.4f ± %.4f", 
                         summary_stats$BayesR_rmse_train_tbv_mean, 
                         summary_stats$BayesR_rmse_train_tbv_sd), "\n")
  cat("BayesA: ", sprintf("%.4f ± %.4f", 
                         summary_stats$BayesA_rmse_train_tbv_mean, 
                         summary_stats$BayesA_rmse_train_tbv_sd), "\n")
  
  cat("\nTRAINING R² ± SD (TBV):\n")
  cat("GBLUP:  ", sprintf("%.4f ± %.4f", 
                         summary_stats$GBLUP_r2_train_tbv_mean, 
                         summary_stats$GBLUP_r2_train_tbv_sd), "\n")
  cat("BayesR: ", sprintf("%.4f ± %.4f", 
                         summary_stats$BayesR_r2_train_tbv_mean, 
                         summary_stats$BayesR_r2_train_tbv_sd), "\n")
  cat("BayesA: ", sprintf("%.4f ± %.4f", 
                         summary_stats$BayesA_r2_train_tbv_mean, 
                         summary_stats$BayesA_r2_train_tbv_sd), "\n")
  
  cat("\nTRAINING Regression Slope (b) ± SD (TBV):\n")
  cat("GBLUP:  ", sprintf("%.4f ± %.4f", 
                        summary_stats$GBLUP_b_train_tbv_mean, 
                        summary_stats$GBLUP_b_train_tbv_sd), "\n")
  cat("BayesR: ", sprintf("%.4f ± %.4f", 
                        summary_stats$BayesR_b_train_tbv_mean, 
                        summary_stats$BayesR_b_train_tbv_sd), "\n")
  cat("BayesA: ", sprintf("%.4f ± %.4f", 
                        summary_stats$BayesA_b_train_tbv_mean, 
                        summary_stats$BayesA_b_train_tbv_sd), "\n")
  
  cat("\nEfficiency (Core-Hours / Test Accuracy) ± SD:\n")
  cat("GBLUP:  ", sprintf("%.2f ± %.2f", 
                         summary_stats$GBLUP_efficiency_mean, 
                         summary_stats$GBLUP_efficiency_sd), 
      " (Lower is better)\n")
  cat("BayesR: ", sprintf("%.2f ± %.2f", 
                         summary_stats$BayesR_efficiency_mean, 
                         summary_stats$BayesR_efficiency_sd), 
      " (Lower is better)\n")
  cat("BayesA: ", sprintf("%.2f ± %.2f", 
                         summary_stats$BayesA_efficiency_mean, 
                         summary_stats$BayesA_efficiency_sd), 
      " (Lower is better)\n")
}