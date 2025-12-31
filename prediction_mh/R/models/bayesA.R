# BayesA Model with Marker-Specific Variance

run_bayesA <- function(matrices, split, gblup_varcomp, config, fold = 0) {
  
  W_train <- matrices$W_train
  W_test <- matrices$W_test
  WtW_diag <- matrices$WtW_diag
  Wty <- matrices$Wty
  
  y <- split$train$phenotypes$phenotype
  n_alleles <- ncol(W_train)
  
  # Hyperparameters
  nu <- config$bayesA$nu
  sigma2_g <- gblup_varcomp$sigma2_g
  sigma2_e_init <- gblup_varcomp$sigma2_e
  S_squared <- sigma2_g / n_alleles
  
  a0_e <- config$bayesA$prior_df_residual
  b0_e <- sigma2_e_init * (a0_e - 1)

  prior_params <- list(
    a0_e = a0_e,
    b0_e = b0_e
  )
  
  mcmc_params <- list(
    n_iter = as.integer(config$mcmc$n_iter),
    n_burn = as.integer(config$mcmc$n_burn),
    n_thin = as.integer(config$mcmc$n_thin),
    seed = as.integer(config$mcmc$seed)
  )

  cat("BayesA Initialization Check:\n")
  cat("n_alleles:", n_alleles, "\n")
  cat("Storage dimensions will be:", (config$mcmc$n_iter - config$mcmc$n_burn) / config$mcmc$n_thin, "x", n_alleles, "\n\n")
  
  # === CALL RUST BACKEND ===
  cat("Running BayesA MCMC...\n")
  
  rust_result <- masbayes::run_bayesa_mcmc(
    w = W_train,
    y = as.numeric(y),
    wtw_diag = as.numeric(WtW_diag),
    wty = as.numeric(Wty),
    nu = as.numeric(nu),
    s_squared = as.numeric(S_squared),
    sigma2_e_init = as.numeric(sigma2_e_init),
    prior_params = prior_params,
    mcmc_params = mcmc_params,
    fold_id = as.integer(fold)
  )
  
  cat("\nMCMC completed!\n")

  beta_A_samples <- rust_result$beta_samples
  sigma2_j_samples <- rust_result$sigma2_j_samples
  sigma2_e_A_samples <- rust_result$sigma2_e_samples

  sigma2_e_A_mcmc <- mcmc(sigma2_e_A_samples)
  ess_A <- effectiveSize(sigma2_e_A_mcmc)
  geweke_A <- geweke.diag(sigma2_e_A_mcmc)$z
  
  # Posterior means
  beta_A_hat <- colMeans(beta_A_samples, na.rm = TRUE)
  sigma2_j_hat <- colMeans(sigma2_j_samples, na.rm = TRUE)
  sigma2_e_A_hat <- mean(sigma2_e_A_samples, na.rm = TRUE)
  
  # Direct marker effect prediction
  GEBV_train <- W_train %*% beta_A_hat
  GEBV_test <- W_test %*% beta_A_hat
  
  # Calculate heritability
  sigma2_g_bayesa <- var(GEBV_train)
  h2_bayesa <- sigma2_g_bayesa / (sigma2_g_bayesa + sigma2_e_A_hat)

  cat("\nBayesA Posterior Means:\n")
  cat("σ²_e:", round(sigma2_e_A_hat, 6), "\n")
  cat("Mean σ²_j:", round(mean(sigma2_j_hat), 8), "\n")
  cat("Range σ²_j:", round(range(sigma2_j_hat), 8), "\n\n")
  cat("σ²_g (genetic):", round(sigma2_g_bayesa, 6), "\n")
  cat("h² (heritability):", round(h2_bayesa, 4), "\n\n")
  
  # Accuracy
  pheno_train <- split$train$phenotypes
  pheno_test <- split$test$phenotypes
  
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
  rmse_train_pheno <- sqrt(mean((GEBV_train - pheno_train$phenotype)^2, na.rm=TRUE))
  rmse_train_tbv <- sqrt(mean((GEBV_train - pheno_train$tbv)^2, na.rm=TRUE))
  rmse_test_pheno <- sqrt(mean((GEBV_test - pheno_test$phenotype)^2, na.rm=TRUE))
  rmse_test_tbv <- sqrt(mean((GEBV_test - pheno_test$tbv)^2, na.rm=TRUE))

  # R² calculations
  r2_train_pheno <- cor(GEBV_train, pheno_train$phenotype, use="complete.obs")^2
  r2_train_tbv <- cor(GEBV_train, pheno_train$tbv, use="complete.obs")^2
  r2_test_pheno <- cor(GEBV_test, pheno_test$phenotype, use="complete.obs")^2
  r2_test_tbv <- cor(GEBV_test, pheno_test$tbv, use="complete.obs")^2

  # Marker statistics
  ## Compute Credible Intervals
  beta_ci_lower <- apply(beta_A_samples, 2, quantile, probs = 0.025, na.rm = TRUE)
  beta_ci_upper <- apply(beta_A_samples, 2, quantile, probs = 0.975, na.rm = TRUE)

  ## Identify Important Alleles
  significant <- (beta_ci_lower > 0 | beta_ci_upper < 0)

  marker_info <- tryCatch({
    m_names <- colnames(W_train)
    if(is.null(m_names)) m_names <- paste0("M", 1:ncol(W_train))

    data.frame(
      allele_id = m_names,
      beta_hat = as.vector(beta_A_hat),
      beta_abs = abs(as.vector(beta_A_hat)),
      ci_lower = as.vector(beta_ci_lower),
      ci_upper = as.vector(beta_ci_upper),
      sigma2_j = as.vector(sigma2_j_hat),
      significant = as.logical(significant),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    warning("Marker info for BayesA failed: ", e$message)
    return(NULL)
  })
  
  list(
    predictions = list(
      train = as.vector(GEBV_train),
      test = as.vector(GEBV_test)
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
    varcomp = list(
      sigma2_e = sigma2_e_A_hat,
      mean_sigma2_j = mean(sigma2_j_hat),
      range_sigma2_j = range(sigma2_j_hat),
      sigma2_g = sigma2_g_bayesa,
      h2 = h2_bayesa
    ),
    markers = marker_info,
    samples = if(config$save_mcmc_samples) {
      list(beta = beta_A_samples, sigma2_j = sigma2_j_samples)
    } else NULL,
    diagnostics = list(
      ess = ess_A,
      geweke_z = geweke_A
    )
  )
}