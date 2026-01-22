# visualization.R
# Visualization functions for model comparison

suppressPackageStartupMessages({
  library(tidyverse)
  library(gridExtra)
})

plot_scatter <- function(all_predictions_combined, config) {
  available_models <- c()
  pred_cols <- c()
  
  if(!all(is.na(all_predictions_combined$GBLUP))) {
    available_models <- c(available_models, "GBLUP")
    pred_cols <- c(pred_cols, all_predictions_combined$GBLUP)
  }
  if(!all(is.na(all_predictions_combined$BayesA))) {
    available_models <- c(available_models, "BayesA")
    pred_cols <- c(pred_cols, all_predictions_combined$BayesA)
  }
  if(!all(is.na(all_predictions_combined$BayesR))) {
    available_models <- c(available_models, "BayesR")
    pred_cols <- c(pred_cols, all_predictions_combined$BayesR)
  }
  
  n_models <- length(available_models)
  if(n_models == 0) {
    cat("WARNING: No valid predictions found for plotting\n")
    return(NULL)
  }
  
  scatter_data <- data.frame(
    True_TBV = rep(all_predictions_combined$True_TBV, n_models),
    Predicted = pred_cols,
    Model = rep(available_models, each = nrow(all_predictions_combined))
  )
  
  p_scatter <- ggplot(scatter_data, aes(x = True_TBV, y = Predicted)) +
    geom_point(alpha = 0.5, size = 1.5, color = "#2E86AB") +
    geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, color = "#E63946") +
    facet_wrap(~ Model, scales = "free", ncol = 3) +
    theme_minimal() +
    labs(title = "Predicted vs true breeding values",
         x = "True breeding value (TBV)",
         y = "Predicted breeding value") +
    theme(strip.text = element_text(face = "bold"))
  
  ggsave(file.path(config$output_dir, "cv_predicted_vs_true_scatter.png"), 
         plot = p_scatter,
         width = 18, height = 8, dpi = 300, bg = "white")
}

plot_performance_bars <- function(results_pheno, config) {
  p_corr <- ggplot(results_pheno, aes(x = Model, y = Correlation, fill = Model)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = Correlation - Correlation_SD,
                      ymax = Correlation + Correlation_SD),
                  width = 0.2) +
    theme_minimal() +
    labs(title = "Model prediction accuracy",
         x = "Model", y = "Correlation with Phenotype") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    geom_text(aes(label = sprintf("%.3f±%.3f", Correlation, Correlation_SD)), 
              vjust = -0.5, size = 3)
  
  p_rmse <- ggplot(results_pheno, aes(x = Model, y = RMSE, fill = Model)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = RMSE - RMSE_SD,
                      ymax = RMSE + RMSE_SD),
                  width = 0.2) +
    theme_minimal() +
    labs(title = "Model prediction error",
         x = "Model",
         y = "RMSE") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    geom_text(aes(label = sprintf("%.3f±%.3f", RMSE, RMSE_SD)),
              vjust = -0.5, size = 3)
  
  p_combined <- grid.arrange(p_corr, p_rmse, ncol = 2)
  
  ggsave(
    filename = file.path(config$output_dir, "cv_performance_comparison_barplots.png"),
    plot = p_combined,
    width = 18,
    height = 6,
    dpi = 300,
    bg = "white"
  )
}

plot_fold_performance <- function(all_cv_results, config) {
  fold_summary <- all_cv_results %>%
    select(Model, Fold, Cor_Pheno, RMSE_Pheno) %>%
    rename(Correlation = Cor_Pheno, RMSE = RMSE_Pheno)
  
  p_fold_corr <- ggplot(fold_summary, aes(x = as.factor(Fold), y = Correlation, 
                                           color = Model, group = Model)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    theme_minimal() +
    labs(title = "Correlation across folds",
         x = "Fold",
         y = "Correlation") +
    theme(legend.position = "bottom")
  
  p_fold_rmse <- ggplot(fold_summary, aes(x = as.factor(Fold), y = RMSE, 
                                          color = Model, group = Model)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    theme_minimal() +
    labs(title = "RMSE across folds",
         x = "Fold",
         y = "RMSE") +
    theme(legend.position = "bottom")
  
  p_fold_combined <- grid.arrange(p_fold_corr, p_fold_rmse, ncol = 2)
  
  ggsave(
    filename = file.path(config$output_dir, "cv_performance_by_fold.png"),
    plot = p_fold_combined,
    width = 18,
    height = 6,
    dpi = 300,
    bg = "white"
  )
}

plot_training_time <- function(all_cv_results, config) {
  train_time_data <- all_cv_results %>%
    group_by(Model) %>%
    summarise(Train_Time = mean(Train_Time))
  
  p_time <- ggplot(train_time_data, aes(x = reorder(Model, Train_Time), y = Train_Time, fill = Model)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    theme_minimal() +
    labs(title = "Training time comparison",
         x = "Model",
         y = "Training time (seconds)") +
    theme(legend.position = "none") +
    geom_text(aes(label = round(Train_Time, 2)), hjust = -0.1)
  
  ggsave(file.path(config$output_dir, "cv_training_time_comparison.png"), 
         plot = p_time,
         width = 15, height = 6, dpi = 300, bg = "white")
}

save_all_plots <- function(all_predictions_combined, all_cv_results, aggregated_results, config) {
  plot_scatter(all_predictions_combined, config)
  plot_performance_bars(aggregated_results$results_pheno, config)
  plot_fold_performance(all_cv_results, config)
  plot_training_time(all_cv_results, config)
}