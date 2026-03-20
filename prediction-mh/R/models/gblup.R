# GBLUP Model with GREML Variance Component Estimation

run_gblup <- function(matrices, split, config) {
  
  # Prepare phenotype data for reference and test
  pheno_train <- split$train$phenotypes
  pheno_test <- split$test$phenotypes

  pheno_train <- pheno_train %>%
    dplyr::select(individual_id, phenotype, tbv) %>%
    dplyr::rename(y = phenotype)

  pheno_train$individual_id <- as.factor(pheno_train$individual_id)

  pheno_test <- pheno_test %>%
    dplyr::select(individual_id, phenotype, tbv) %>%
    dplyr::rename(y = phenotype)
  
  pheno_test$individual_id <- as.factor(pheno_test$individual_id)
  
  # Add rownames/colnames to A_train
  A_train <- matrices$A_train
  A_test <- matrices$A_combined[
    as.character(pheno_test$individual_id),
    as.character(pheno_test$individual_id)
  ]
  rownames(A_train) <- pheno_train$individual_id
  colnames(A_train) <- pheno_train$individual_id
  rownames(A_test) <- levels(pheno_test$individual_id)
  colnames(A_test) <- levels(pheno_test$individual_id)
  
  # GREML variance component estimation
  model_greml <- mmes(
    fixed = y ~ 1,
    random = ~ vsm(ism(individual_id), Gu = A_train),
    rcov = ~ vsm(ism(units)),
    data = pheno_train,
    verbose = FALSE
  )

  # Extract variance components
  summary(model_greml)

  h2_test_gblup <- tryCatch({
  model_test <- mmes(
    fixed = y ~ 1,
    random = ~ vsm(ism(individual_id), Gu = A_test),
    rcov = ~ vsm(ism(units)),
    data = pheno_test,
    verbose = FALSE
  )
  vc <- model_test$theta
  sg <- vc[[1]][1,1]; se <- vc[[2]][1,1]
  h2_raw <- sg / (sg + se)
  if (is.na(h2_raw) || h2_raw < 0 || h2_raw > 1) NA else h2_raw
}, error = function(e) NA)
  
  # Calculate heritability
  var_comp <- model_greml$theta
  sigma2_g <- var_comp[[1]][1,1]
  sigma2_e <- var_comp[[2]][1,1]
  h2 <- sigma2_g / (sigma2_g + sigma2_e)
  
  # Extract BLUP
  if(is.list(model_greml$u)) {
    blup <- model_greml$u[[1]]
  } else {
    blup <- model_greml$u
  }
  if(is.null(dim(blup))) blup <- matrix(blup, ncol = 1)
  
  # Extract IDs
  if(!is.null(rownames(blup))) {
    blup_ids <- rownames(blup)
  } else if(!is.null(dimnames(blup)[[1]])) {
    blup_ids <- dimnames(blup)[[1]]
  } else {
    blup_ids <- as.character(pheno_train$individual_id)
  }
  
  # Match and reorder to phenotype ID order
  matched_idx <- match(as.character(pheno_train$individual_id), blup_ids)
  mu_hat <- mean(pheno_train$y)
  pred_train <- as.vector(blup[matched_idx]) + mu_hat
  
  # Extract relationship blocks
  A_ref_ref <- matrices$A_combined[
    as.character(pheno_train$individual_id), 
    as.character(pheno_train$individual_id)
  ]
  A_test_ref <- matrices$A_combined[
    as.character(pheno_test$individual_id), 
    as.character(pheno_train$individual_id)
  ]
  
  lambda <- sigma2_e / sigma2_g
  A_ref_inv <- solve(A_ref_ref + diag(lambda, nrow(A_ref_ref)))
  pred_test <- as.vector(A_test_ref %*% A_ref_inv %*% (pred_train - mu_hat)) + mu_hat
  
  # Accuracy
  accuracy_train_pheno <- cor(pred_train, pheno_train$y)
  accuracy_train_tbv <- cor(pred_train, pheno_train$tbv)
  accuracy_test_pheno <- cor(pred_test, pheno_test$y) #phenotype
  accuracy_test_tbv <- cor(pred_test, pheno_test$tbv)
  accuracy_test_gebv <- if (!is.na(h2_test_gblup)) accuracy_test_pheno / sqrt(h2_test_gblup) else NA

  b_train_pheno <- calculate_regression_slope(pred_train, pheno_train$y)
  b_train_tbv <- calculate_regression_slope(pred_train, pheno_train$tbv)
  b_test_pheno <- calculate_regression_slope(pred_test, pheno_test$y)
  b_test_tbv <- calculate_regression_slope(pred_test, pheno_test$tbv)

  # RMSE calculations
  rmse_train_pheno <- sqrt(mean((pred_train - pheno_train$y)^2, na.rm=TRUE))
  rmse_train_tbv <- sqrt(mean((pred_train - pheno_train$tbv)^2, na.rm=TRUE))
  rmse_test_pheno <- sqrt(mean((pred_test - pheno_test$y)^2, na.rm=TRUE))
  rmse_test_tbv <- sqrt(mean((pred_test - pheno_test$tbv)^2, na.rm=TRUE))

  # R² calculations
  r2_train_pheno <- cor(pred_train, pheno_train$y)^2
  r2_train_tbv <- cor(pred_train, pheno_train$tbv)^2
  r2_test_pheno <- cor(pred_test, pheno_test$y)^2
  r2_test_tbv <- cor(pred_test, pheno_test$tbv)^2
  
  list(
    predictions = list(
      train = pred_train,
      test = pred_test
    ),
    accuracy = list(
      train_pheno = accuracy_train_pheno,
      train_tbv = accuracy_train_tbv,
      test_pheno = accuracy_test_pheno,
      test_tbv = accuracy_test_tbv,
      test_gebv = accuracy_test_gebv
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
    varcomp = list(
      sigma2_g = sigma2_g,
      sigma2_e = sigma2_e,
      h2 = h2,
      h2_test = h2_test_gblup,
      lambda = lambda,
      mean_diag_train = mean(diag(A_train)),
      mean_diag_combined = mean(diag(matrices$A_combined))
    ),
    model = model_greml,
    n_cores = 1,
    runtime = NA
  )
}