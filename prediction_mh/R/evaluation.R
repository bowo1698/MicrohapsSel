# Model Evaluation and Results Aggregation Functions

calculate_regression_slope <- function(predicted, observed) {
  fit <- lm(observed ~ predicted)
  coef(fit)[2]  # Return slope (b)
}

aggregate_fold_results <- function(gblup_res, bayesR_res, bayesA_res,
                                   xgb_res, fold_id) {
  
  list(
    fold = fold_id,
    gblup_train_pheno = gblup_res$accuracy$train_pheno,
    gblup_train_tbv = gblup_res$accuracy$train_tbv,
    gblup_test_pheno = gblup_res$accuracy$test_pheno,
    gblup_test_tbv = gblup_res$accuracy$test_tbv,
    gblup_h2 = gblup_res$varcomp$h2,
    bayesR_train_pheno = bayesR_res$accuracy$train_pheno,
    bayesR_train_tbv = bayesR_res$accuracy$train_tbv,
    bayesR_test_pheno = bayesR_res$accuracy$test_pheno,
    bayesR_test_tbv = bayesR_res$accuracy$test_tbv,
    bayesR_h2 = bayesR_res$varcomp$h2,
    bayesA_train_pheno = bayesA_res$accuracy$train_pheno,
    bayesA_train_tbv = bayesA_res$accuracy$train_tbv,
    bayesA_test_pheno = bayesA_res$accuracy$test_pheno,
    bayesA_test_tbv = bayesA_res$accuracy$test_tbv,
    bayesA_h2 = bayesA_res$varcomp$h2,
    xgb_train_pheno = xgb_res$accuracy$train_pheno,
    xgb_train_tbv = xgb_res$accuracy$train_tbv,
    xgb_test_pheno = xgb_res$accuracy$test_pheno,
    xgb_test_tbv = xgb_res$accuracy$test_tbv,
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
    xgb_b_train_pheno = xgb_res$regression_slopes$train_pheno,
    xgb_b_train_tbv = xgb_res$regression_slopes$train_tbv,
    xgb_b_test_pheno = xgb_res$regression_slopes$test_pheno,
    xgb_b_test_tbv = xgb_res$regression_slopes$test_tbv,
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
    xgb_rmse_train_pheno = xgb_res$rmse$train_pheno,
    xgb_rmse_train_tbv = xgb_res$rmse$train_tbv,
    xgb_rmse_test_pheno = xgb_res$rmse$test_pheno,
    xgb_rmse_test_tbv = xgb_res$rmse$test_tbv,
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
    xgb_r2_train_pheno = xgb_res$r2$train_pheno,
    xgb_r2_train_tbv = xgb_res$r2$train_tbv,
    xgb_r2_test_pheno = xgb_res$r2$test_pheno,
    xgb_r2_test_tbv = xgb_res$r2$test_tbv
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
      XGBoost_train_pheno_mean = mean(xgb_train_pheno),
      XGBoost_train_pheno_sd = sd(xgb_train_pheno),
      XGBoost_train_tbv_mean = mean(xgb_train_tbv),
      XGBoost_train_tbv_sd = sd(xgb_train_tbv),
      XGBoost_test_pheno_mean = mean(xgb_test_pheno),
      XGBoost_test_pheno_sd = sd(xgb_test_pheno),
      XGBoost_test_tbv_mean = mean(xgb_test_tbv),
      XGBoost_test_tbv_sd = sd(xgb_test_tbv),
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
      XGBoost_b_train_pheno_mean = mean(xgb_b_train_pheno),
      XGBoost_b_train_pheno_sd = sd(xgb_b_train_pheno),
      XGBoost_b_train_tbv_mean = mean(xgb_b_train_tbv),
      XGBoost_b_train_tbv_sd = sd(xgb_b_train_tbv),
      XGBoost_b_test_pheno_mean = mean(xgb_b_test_pheno),
      XGBoost_b_test_pheno_sd = sd(xgb_b_test_pheno),
      XGBoost_b_test_tbv_mean = mean(xgb_b_test_tbv),
      XGBoost_b_test_tbv_sd = sd(xgb_b_test_tbv),
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
      XGBoost_rmse_train_pheno_mean = mean(xgb_rmse_train_pheno),
      XGBoost_rmse_train_pheno_sd = sd(xgb_rmse_train_pheno),
      XGBoost_rmse_test_pheno_mean = mean(xgb_rmse_test_pheno),
      XGBoost_rmse_test_pheno_sd = sd(xgb_rmse_test_pheno),
      XGBoost_rmse_train_tbv_mean = mean(xgb_rmse_train_tbv),
      XGBoost_rmse_train_tbv_sd = sd(xgb_rmse_train_tbv),
      XGBoost_rmse_test_tbv_mean = mean(xgb_rmse_test_tbv),
      XGBoost_rmse_test_tbv_sd = sd(xgb_rmse_test_tbv),
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
      XGBoost_r2_train_pheno_mean = mean(xgb_r2_train_pheno),
      XGBoost_r2_train_pheno_sd = sd(xgb_r2_train_pheno),
      XGBoost_r2_test_pheno_mean = mean(xgb_r2_test_pheno),
      XGBoost_r2_test_pheno_sd = sd(xgb_r2_test_pheno),
      XGBoost_r2_train_tbv_mean = mean(xgb_r2_train_tbv),
      XGBoost_r2_train_tbv_sd = sd(xgb_r2_train_tbv),
      XGBoost_r2_test_tbv_mean = mean(xgb_r2_test_tbv),
      XGBoost_r2_test_tbv_sd = sd(xgb_r2_test_tbv)
    )
  
  list(
    results = results_df,
    summary = summary_stats
  )
}

save_gebv_with_regression <- function(fold_id, gblup_res, bayesR_res, 
                                      bayesA_res, xgb_res, split, output_dir) {
  
  fold_dir <- file.path(output_dir, paste0("fold_", fold_id))
  
  # Training set
  train_gebv <- data.frame(
    individual_id = split$train$phenotypes$individual_id,
    tbv = split$train$phenotypes$tbv,
    phenotype = split$train$phenotypes$phenotype,
    GEBV_GBLUP = gblup_res$predictions$train,
    GEBV_BayesR = bayesR_res$predictions$train,
    GEBV_BayesA = bayesA_res$predictions$train,
    GEBV_XGBoost = xgb_res$predictions$train,
    b_GBLUP_tbv = gblup_res$regression_slopes$train_tbv,
    b_BayesR_tbv = bayesR_res$regression_slopes$train_tbv,
    b_BayesA_tbv = bayesA_res$regression_slopes$train_tbv,
    b_XGBoost_tbv = xgb_res$regression_slopes$train_tbv
  )
  
  # Test set
  test_gebv <- data.frame(
    individual_id = split$test$phenotypes$individual_id,
    tbv = split$test$phenotypes$tbv,
    phenotype = split$test$phenotypes$phenotype,
    GEBV_GBLUP = gblup_res$predictions$test,
    GEBV_BayesR = bayesR_res$predictions$test,
    GEBV_BayesA = bayesA_res$predictions$test,
    GEBV_XGBoost = xgb_res$predictions$test,
    b_GBLUP_tbv = gblup_res$regression_slopes$test_tbv,
    b_BayesR_tbv = bayesR_res$regression_slopes$test_tbv,
    b_BayesA_tbv = bayesA_res$regression_slopes$test_tbv,
    b_XGBoost_tbv = xgb_res$regression_slopes$test_tbv
  )
  
  write.csv(train_gebv, file.path(fold_dir, "train_gebv.csv"), row.names = FALSE)
  write.csv(test_gebv, file.path(fold_dir, "test_gebv.csv"), row.names = FALSE)
}

save_fold_results <- function(fold_id, gblup_res, bayesR_res, bayesA_res,
                              xgb_res, matrices, output_dir) {
  
  fold_dir <- file.path(output_dir, paste0("fold_", fold_id))
  dir.create(fold_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save aggregated results
  fold_results <- aggregate_fold_results(
    gblup_res, bayesR_res, bayesA_res, xgb_res, fold_id
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
  
  # Save XGBoost importance
  if(!is.null(xgb_res$importance)) {
    write.csv(
      xgb_res$importance,
      file.path(fold_dir, "xgboost_feature_importance.csv"),
      row.names = FALSE
    )
    
    if(!is.null(xgb_res$high_importance) && nrow(xgb_res$high_importance) > 0) {
      write.csv(
        xgb_res$high_importance,
        file.path(fold_dir, "xgboost_high_importance_markers.csv"),
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
  cat("XGBoost:", sprintf("%.4f ± %.4f", 
                         summary_stats$XGBoost_test_tbv_mean, 
                         summary_stats$XGBoost_test_tbv_sd), "\n")
  
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
  cat("XGBoost:", sprintf("%.4f ± %.4f", 
                         summary_stats$XGBoost_rmse_train_tbv_mean, 
                         summary_stats$XGBoost_rmse_train_tbv_sd), "\n")
  
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
  cat("XGBoost:", sprintf("%.4f ± %.4f", 
                         summary_stats$XGBoost_r2_train_tbv_mean, 
                         summary_stats$XGBoost_r2_train_tbv_sd), "\n")
  
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
  cat("XGBoost:", sprintf("%.4f ± %.4f", 
                        summary_stats$XGBoost_b_train_tbv_mean, 
                        summary_stats$XGBoost_b_train_tbv_sd), "\n")
}