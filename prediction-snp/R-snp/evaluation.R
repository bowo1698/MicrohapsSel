# evaluation.R
# Model evaluation and metrics calculation functions

suppressPackageStartupMessages({
  library(tidyverse)
})

evaluate_model <- function(predicted, true_tbv, true_pheno, model_name) {
  r_tbv <- cor(predicted, true_tbv, use = "complete.obs")
  r2_tbv <- r_tbv^2
  rmse_tbv <- sqrt(mean((predicted - true_tbv)^2, na.rm = TRUE))
  mae_tbv <- mean(abs(predicted - true_tbv), na.rm = TRUE)
  bias_lm_tbv <- lm(predicted ~ true_tbv)
  bias_slope_tbv <- coef(bias_lm_tbv)[2]
  rank_cor_tbv <- cor(predicted, true_tbv, method = "spearman", use = "complete.obs")
  
  r_pheno <- cor(predicted, true_pheno, use = "complete.obs")
  r2_pheno <- r_pheno^2
  rmse_pheno <- sqrt(mean((predicted - true_pheno)^2, na.rm = TRUE))
  mae_pheno <- mean(abs(predicted - true_pheno), na.rm = TRUE)
  bias_lm_pheno <- lm(predicted ~ true_pheno)
  bias_slope_pheno <- coef(bias_lm_pheno)[2]
  rank_cor_pheno <- cor(predicted, true_pheno, method = "spearman", use = "complete.obs")
  
  data.frame(
    Model = model_name,
    Cor_TBV = r_tbv,
    R2_TBV = r2_tbv,
    RMSE_TBV = rmse_tbv,
    MAE_TBV = mae_tbv,
    Bias_Slope_TBV = bias_slope_tbv,
    Rank_Cor_TBV = rank_cor_tbv,
    Cor_Pheno = r_pheno,
    R2_Pheno = r2_pheno,
    RMSE_Pheno = rmse_pheno,
    MAE_Pheno = mae_pheno,
    Bias_Slope_Pheno = bias_slope_pheno,
    Rank_Cor_Pheno = rank_cor_pheno
  )
}

evaluate_model_train <- function(predicted, true_tbv, true_pheno, model_name) {
  r_tbv <- cor(predicted, true_tbv, use = "complete.obs")
  rmse_tbv <- sqrt(mean((predicted - true_tbv)^2, na.rm = TRUE))
  
  r_pheno <- cor(predicted, true_pheno, use = "complete.obs")
  rmse_pheno <- sqrt(mean((predicted - true_pheno)^2, na.rm = TRUE))
  
  data.frame(
    Cor_TBV_Train = r_tbv,
    RMSE_TBV_Train = rmse_tbv,
    Cor_Pheno_Train = r_pheno,
    RMSE_Pheno_Train = rmse_pheno
  )
}

aggregate_cv_results <- function(all_cv_results) {
  results_tbv <- do.call(rbind, lapply(split(all_cv_results, all_cv_results$Model), function(x) {
    data.frame(
      Model = unique(x$Model),
      Correlation = mean(x$Cor_TBV),
      Correlation_SD = sd(x$Cor_TBV),
      R2 = mean(x$R2_TBV),
      RMSE = mean(x$RMSE_TBV),
      RMSE_SD = sd(x$RMSE_TBV),
      MAE = mean(x$MAE_TBV),
      Bias_Slope = mean(x$Bias_Slope_TBV),
      Rank_Correlation = mean(x$Rank_Cor_TBV),
      Cor_Train = mean(x$Cor_TBV_Train),
      RMSE_Train = mean(x$RMSE_TBV_Train)
    )
  }))
  rownames(results_tbv) <- NULL
  
  results_pheno <- do.call(rbind, lapply(split(all_cv_results, all_cv_results$Model), function(x) {
    data.frame(
      Model = unique(x$Model),
      Correlation = mean(x$Cor_Pheno),
      Correlation_SD = sd(x$Cor_Pheno),
      R2 = mean(x$R2_Pheno),
      RMSE = mean(x$RMSE_Pheno),
      RMSE_SD = sd(x$RMSE_Pheno),
      MAE = mean(x$MAE_Pheno),
      Bias_Slope = mean(x$Bias_Slope_Pheno),
      Rank_Correlation = mean(x$Rank_Cor_Pheno),
      Cor_Train = mean(x$Cor_Pheno_Train),
      RMSE_Train = mean(x$RMSE_Pheno_Train)
    )
  }))
  rownames(results_pheno) <- NULL
  
  train_time_summary <- do.call(rbind, lapply(split(all_cv_results, all_cv_results$Model), function(x) {
    data.frame(
      Model = unique(x$Model),
      Train_Time = mean(x$Train_Time),
      Pred_Time = mean(x$Pred_Time)
    )
  }))
  rownames(train_time_summary) <- NULL
  
  results_tbv <- merge(results_tbv, train_time_summary, by = "Model")
  results_pheno <- merge(results_pheno, train_time_summary, by = "Model")
  
  results_ranked_tbv <- results_tbv[order(-results_tbv$Correlation), ]
  results_ranked_pheno <- results_pheno[order(-results_pheno$Correlation), ]
  
  list(
    results_tbv = results_tbv,
    results_pheno = results_pheno,
    results_ranked_tbv = results_ranked_tbv,
    results_ranked_pheno = results_ranked_pheno
  )
}

save_results <- function(all_cv_results, all_predictions_combined, aggregated_results, config) {
  write.csv(aggregated_results$results_pheno, 
            file.path(config$output_dir, "model_performance_vs_Phenotype.csv"),
            row.names = FALSE)
  
  write.csv(aggregated_results$results_tbv, 
            file.path(config$output_dir, "model_performance_vs_TBV.csv"),
            row.names = FALSE)
  
  write.csv(aggregated_results$results_ranked_pheno, 
            file.path(config$output_dir, "models_ranked_by_Phenotype.csv"),
            row.names = FALSE)
  
  write.csv(aggregated_results$results_ranked_tbv, 
            file.path(config$output_dir, "models_ranked_by_TBV.csv"),
            row.names = FALSE)
  
  write.csv(all_predictions_combined, 
            file.path(config$output_dir, "cv_predictions_all_folds.csv"), 
            row.names = FALSE)
  
  write.csv(all_cv_results, 
            file.path(config$output_dir, "cv_results_by_fold.csv"), 
            row.names = FALSE)
}

generate_summary_report <- function(aggregated_results, all_cv_results, n_individuals, n_snps, config) {
  report_output <- capture.output({
    cat(strrep("=", 60), "\n")
    cat("         GENOMIC PREDICTION MODEL COMPARISON REPORT\n")
    cat(strrep("=", 60), "\n\n")
    
    cat("Dataset Information:\n")
    cat("  Reference Population:", n_individuals, "individuals\n")
    cat("  CV Configuration:", config$k_folds, "-fold cross-validation\n")
    cat("  SNPs after QC:", n_snps, "\n")
    cat("  Average fold size:", round(n_individuals / config$k_folds), "individuals\n")
    
    cat("\n=== vs TBV ===\n")
    cat("  Best Model (Correlation):", aggregated_results$results_ranked_tbv$Model[1], 
        "(r =", round(aggregated_results$results_ranked_tbv$Correlation[1], 3), ")\n")
    cat("  Best Model (RMSE):", aggregated_results$results_tbv[which.min(aggregated_results$results_tbv$RMSE), ]$Model,
        "(RMSE =", round(min(aggregated_results$results_tbv$RMSE), 3), ")\n")
    
    cat("\n=== vs Phenotype ===\n")
    cat("  Best Model (Correlation):", aggregated_results$results_ranked_pheno$Model[1], 
        "(r =", round(aggregated_results$results_ranked_pheno$Correlation[1], 3), ")\n")
    cat("  Best Model (RMSE):", aggregated_results$results_pheno[which.min(aggregated_results$results_pheno$RMSE), ]$Model,
        "(RMSE =", round(min(aggregated_results$results_pheno$RMSE), 3), ")\n")
    
    train_time_summary <- all_cv_results %>%
      group_by(Model) %>%
      summarise(Train_Time = mean(Train_Time), Pred_Time = mean(Pred_Time))
    
    cat("\n=== Computational Efficiency ===\n")
    cat("  Fastest Training:", train_time_summary[which.min(train_time_summary$Train_Time), ]$Model, 
        "(", round(min(train_time_summary$Train_Time), 2), "sec)\n")
    cat("  Fastest Prediction:", train_time_summary[which.min(train_time_summary$Pred_Time), ]$Model, 
        "(", round(min(train_time_summary$Pred_Time), 4), "sec)\n\n")
    
    cat("\nDetailed Performance vs TBV:\n")
    print(aggregated_results$results_tbv)
    
    cat("\n\nDetailed Performance vs Phenotype:\n")
    print(aggregated_results$results_pheno)
    
    cat(strrep("=", 60), "\n")
  })
  
  writeLines(report_output, file.path(config$output_dir, "summary_report.txt"))
}