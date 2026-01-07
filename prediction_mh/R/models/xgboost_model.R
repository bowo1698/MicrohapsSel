# XGBoost Model for Genomic Prediction

run_xgboost <- function(matrices, split, config) {
  
  W_train <- matrices$W_train
  W_test <- matrices$W_test
  
  pheno_train <- split$train$phenotypes
  pheno_test <- split$test$phenotypes
  
  set.seed(config$mcmc$seed)
  
  # Prepare DMatrix
  dtrain <- xgb.DMatrix(data = W_train, label = pheno_train$phenotype)
  dtest <- xgb.DMatrix(data = W_test, label = pheno_test$phenotype)
  
  # Set parameters
  xgb_params <- list(
    objective = "reg:squarederror",
    max_depth = config$xgboost$max_depth,
    eta = config$xgboost$eta,
    lambda = config$xgboost$lambda,
    colsample_bytree = config$xgboost$colsample_bytree,
    subsample = config$xgboost$subsample,
    nthread = 1
  )
  
  # Train model
  model_xgb <- xgb.train(
    params = xgb_params,
    data = dtrain,
    nrounds = config$xgboost$nrounds,
    verbose = 0
  )
  
  # Predictions
  GEBV_train <- predict(model_xgb, dtrain)
  GEBV_test <- predict(model_xgb, dtest)
  
  # Accuracy
  accuracy_train_pheno <- cor(GEBV_train, pheno_train$phenotype, use="complete.obs")
  accuracy_train_tbv <- cor(GEBV_train, pheno_train$tbv, use="complete.obs")
  accuracy_test_pheno <- cor(GEBV_test, pheno_test$phenotype, use="complete.obs")
  accuracy_test_tbv <- cor(GEBV_test, pheno_test$tbv, use="complete.obs")

  # Regression slopes
  b_train_pheno <- calculate_regression_slope(GEBV_train, pheno_train$phenotype)
  b_train_tbv <- calculate_regression_slope(GEBV_train, pheno_train$tbv)
  b_test_pheno <- calculate_regression_slope(GEBV_test, pheno_test$phenotype)
  b_test_tbv <- calculate_regression_slope(GEBV_test, pheno_test$tbv)

  # RMSE calculations
  rmse_train_pheno <- sqrt(mean((GEBV_train - pheno_train$phenotype)^2))
  rmse_train_tbv <- sqrt(mean((GEBV_train - pheno_train$tbv)^2))
  rmse_test_pheno <- sqrt(mean((GEBV_test - pheno_test$phenotype)^2))
  rmse_test_tbv <- sqrt(mean((GEBV_test - pheno_test$tbv)^2))

  # R² calculations
  r2_train_pheno <- cor(GEBV_train, pheno_train$phenotype, use="complete.obs")^2
  r2_train_tbv <- cor(GEBV_train, pheno_train$tbv, use="complete.obs")^2
  r2_test_pheno <- cor(GEBV_test, pheno_test$phenotype, use="complete.obs")^2
  r2_test_tbv <- cor(GEBV_test, pheno_test$tbv, use="complete.obs")^2
  
  # Feature importance
  importance_matrix <- xgb.importance(model = model_xgb)
  
  # Map feature indices to allele IDs
  if(is.null(colnames(W_train))) {
    colnames(W_train) <- paste0("allele_", 1:ncol(W_train))
  }

  feature_idx <- suppressWarnings(
    as.integer(gsub("f", "", importance_matrix$Feature))
  )
  feature_idx <- ifelse(is.na(feature_idx), 0, feature_idx) + 1
  feature_idx <- pmin(pmax(feature_idx, 1), ncol(W_train))

  importance_matrix$allele_id <- colnames(W_train)[feature_idx]
  
  # Merge with allele info
  importance_full <- importance_matrix %>%
    left_join(matrices$allele_info, by = "allele_id") %>%
    arrange(desc(Gain))
  
  # High importance markers
  gain_threshold <- quantile(importance_full$Gain, 0.95)
  high_importance <- importance_full %>% filter(Gain > gain_threshold)
  
  list(
    predictions = list(
      train = GEBV_train,
      test = GEBV_test
    ),
    accuracy = list(
      train_pheno = accuracy_train_pheno,
      train_tbv = accuracy_train_tbv,
      test_pheno = accuracy_test_pheno,
      test_tbv = accuracy_test_tbv
    ),
    regression_slopes = list(
      train_pheno = b_train_pheno,
      train_tbv = b_train_tbv,
      test_pheno = b_test_pheno,
      test_tbv = b_test_tbv
    ),
    rmse = list(
      train_pheno = rmse_train_pheno,
      train_tbv = rmse_train_tbv,
      test_pheno = rmse_test_pheno,
      test_tbv = rmse_test_tbv
    ),
    r2 = list(
      train_pheno = r2_train_pheno,
      train_tbv = r2_train_tbv,
      test_pheno = r2_test_pheno,
      test_tbv = r2_test_tbv
    ),
    importance = importance_full,
    high_importance = high_importance,
    model = model_xgb,
    n_cores = 1,
    runtime = NA
  )
}