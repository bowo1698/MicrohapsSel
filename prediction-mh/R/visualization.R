# Visualization Functions for Genomic Prediction Results

plot_accuracy_comparison <- function(results_df, output_dir = NULL) {
  
  results_long <- results_df %>%
    select(fold, gblup_test_tbv, bayesR_test_tbv, bayesA_test_tbv, xgb_test_tbv) %>%
    pivot_longer(-fold, names_to = "Method", values_to = "Accuracy") %>%
    mutate(Method = recode(Method,
                          gblup_test_tbv = "GBLUP",
                          bayesR_test_tbv = "BayesR",
                          bayesA_test_tbv = "BayesA",
                          xgb_test_tbv = "XGBoost"))
  
  p <- ggplot(results_long, aes(x = Method, y = Accuracy, fill = Method)) +
    stat_summary(fun = mean, geom = "bar", alpha = 0.7) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
    labs(
      title = "Prediction accuracy across folds",
      y = "Correlation with TBV",
      x = ""
    ) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 0, hjust = 0.5))
  
  if(!is.null(output_dir)) {
    ggsave(
      file.path(output_dir, "accuracy_comparison_barplot.png"),
      p, width = 10, height = 6
    )
  }
}

plot_heritability_comparison <- function(results_df, output_dir = NULL) {
  
  h2_long <- results_df %>%
    select(fold, gblup_h2, bayesR_h2, bayesA_h2) %>%
    pivot_longer(-fold, names_to = "Method", values_to = "h2") %>%
    mutate(Method = recode(Method,
                          gblup_h2 = "GBLUP",
                          bayesR_h2 = "BayesR",
                          bayesA_h2 = "BayesA"))
  
  p <- ggplot(h2_long, aes(x = Method, y = h2, fill = Method)) +
    stat_summary(fun = mean, geom = "bar", alpha = 0.7) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
    labs(
      title = "Heritability estimates across folds",
      y = "h²",
      x = ""
    ) +
    theme_minimal() +
    theme(legend.position = "none")
  
  if(!is.null(output_dir)) {
    ggsave(
      file.path(output_dir, "heritability_comparison.png"),
      p, width = 10, height = 6
    )
  }
}

plot_fold_variability <- function(results_df, output_dir = NULL) {
  
  results_long <- results_df %>%
    select(fold, gblup_test_tbv, bayesR_test_tbv, bayesA_test_tbv, xgb_test_tbv) %>%
    pivot_longer(-fold, names_to = "Method", values_to = "Accuracy") %>%
    mutate(Method = recode(Method,
                          gblup_test_tbv = "GBLUP",
                          bayesR_test_tbv = "BayesR",
                          bayesA_test_tbv = "BayesA",
                          xgb_test_tbv = "XGBoost"))
  
  p <- ggplot(results_long, aes(x = factor(fold), y = Accuracy, 
                                 color = Method, group = Method)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    labs(
      title = "Prediction accuracy by fold",
      x = "Fold",
      y = "Correlation with TBV"
    ) +
    theme_minimal() +
    theme(legend.position = "right")
  
  if(!is.null(output_dir)) {
    ggsave(
      file.path(output_dir, "fold_variability.png"),
      p, width = 10, height = 6
    )
  }
}

plot_variance_components <- function(results_df, output_dir = NULL) {
  
  vc_long <- results_df %>%
    select(fold, gblup_h2, bayesR_h2, bayesA_h2) %>%
    pivot_longer(-fold, names_to = "Method", values_to = "h2") %>%
    mutate(Method = recode(Method,
                          gblup_h2 = "GBLUP",
                          bayesR_h2 = "BayesR",
                          bayesA_h2 = "BayesA"))
  
  p <- ggplot(vc_long, aes(x = Method, y = h2)) +
    geom_boxplot(aes(fill = Method), alpha = 0.7) +
    geom_jitter(width = 0.1, alpha = 0.5) +
    labs(
      title = "Heritability aistribution across folds",
      x = "",
      y = "h²"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
  
  if(!is.null(output_dir)) {
    ggsave(
      file.path(output_dir, "heritability_boxplot.png"),
      p, width = 10, height = 6
    )
  }
}

plot_train_vs_test_accuracy <- function(results_df, output_dir = NULL) {
  
  train_test <- results_df %>%
    select(fold, 
           gblup_train_tbv, gblup_test_tbv,
           bayesR_train_tbv, bayesR_test_tbv,
           bayesA_train_tbv, bayesA_test_tbv,
           xgb_train_tbv, xgb_test_tbv) %>%
    pivot_longer(-fold, names_to = "metric", values_to = "accuracy") %>%
    separate(metric, into = c("method", "set", "type"), sep = "_") %>%
    mutate(
      method = recode(method,
                     gblup = "GBLUP",
                     bayesR = "BayesR",
                     bayesA = "BayesA",
                     xgb = "XGBoost"),
      set = factor(set, levels = c("train", "test"))
    )
  
  p <- ggplot(train_test, aes(x = method, y = accuracy, fill = set)) +
    geom_boxplot(position = position_dodge(0.8)) +
    labs(
      title = "Training vs test accuracy",
      x = "",
      y = "Correlation with TBV",
      fill = "Dataset"
    ) +
    theme_minimal() +
    scale_fill_manual(values = c("train" = "lightblue", "test" = "coral"))
  
  if(!is.null(output_dir)) {
    ggsave(
      file.path(output_dir, "train_vs_test_accuracy.png"),
      p, width = 10, height = 6
    )
  }
}