# BayesR Model with Mixture Prior

run_bayesR <- function(matrices, split, gblup_varcomp, config, fold = 0) {
  
  W_train <- matrices$W_train
  W_test <- matrices$W_test
  WtW_diag <- matrices$WtW_diag
  Wty <- matrices$Wty
  
  y <- split$train$phenotypes$phenotype
  n_alleles <- ncol(W_train)
  
  # Hyperparameters from config and GBLUP
  pi_vec <- config$bayesR$pi
  sigma2_ah <- gblup_varcomp$sigma2_g
  sigma2_e_init <- gblup_varcomp$sigma2_e
  
  # Variance components per mixture using variance from GBLUP
  sigma2_0 <- 0
  sigma2_small <- config$bayesR$variance_scaling[2] * sigma2_ah
  sigma2_medium <- config$bayesR$variance_scaling[3] * sigma2_ah
  sigma2_large <- config$bayesR$variance_scaling[4] * sigma2_ah
  sigma2_vec <- c(sigma2_0, sigma2_small, sigma2_medium, sigma2_large)
  
  # Prior hyperparameters for Inverse-Gamma
  a0_e <- config$bayesR$prior_df$residual
  b0_e <- sigma2_e_init * (a0_e - 1)
  
  a0_small <- config$bayesR$prior_df$small
  b0_small <- sigma2_vec[2] * (a0_small - 1)
  
  a0_medium <- config$bayesR$prior_df$medium
  b0_medium <- sigma2_vec[3] * (a0_medium - 1)
  
  a0_large <- config$bayesR$prior_df$large
  b0_large <- sigma2_vec[4] * (a0_large - 1)

  prior_params <- list(
    a0_e = a0_e,
    b0_e = b0_e,
    a0_small = a0_small,
    b0_small = b0_small,
    a0_medium = a0_medium,
    b0_medium = b0_medium,
    a0_large = a0_large,
    b0_large = b0_large
  )

  cat("BayesR Settings:\n")
  cat("Algorithm:", toupper(config$bayes_algo), "\n")
  cat("π =", pi_vec, "\n")
  cat("σ²_αh =", round(sigma2_ah, 6), "| σ²_e =", round(sigma2_e_init, 6), "\n")
  cat("σ² mix:", round(sigma2_vec, 8), "\n\n")

  if (config$bayes_algo == "em") {

    cat("Running BayesR with EM algorithm...\n\n")

    res <- masbayes::run_bayesr(
      w = W_train,
      y = y,
      wtw_diag = WtW_diag,
      wty = Wty,
      pi_vec = pi_vec,
      sigma2_vec = sigma2_vec,
      sigma2_e_init = sigma2_e_init,
      em_params = list(
        max_iter = as.integer(config$n_iter),
        tol = config$em_tol
      ),
      method = "em",
      fold_id = as.integer(fold)
    )

    beta_hat <- res$beta_samples[1, ]
    sigma2_e_hat <- res$sigma2_e_samples[1]
    sigma2_small_hat <- res$sigma2_small_samples[1]
    sigma2_medium_hat <- res$sigma2_medium_samples[1]
    sigma2_large_hat <- res$sigma2_large_samples[1]
    pi_hat <- res$pi_samples[1, ]
    gamma_mode <- res$gamma_samples[1, ]

    ess <- NA
    geweke <- NA

  } else if (config$bayes_algo == "mcmc") {
    
    cat("Running BayesR with MCMC algorithm...\n\n")

    res <- masbayes::run_bayesr(
      w = W_train,
      y = y,
      wtw_diag = WtW_diag,
      wty = Wty,
      pi_vec = pi_vec,
      sigma2_vec = sigma2_vec,
      sigma2_e_init = sigma2_e_init,
      sigma2_ah = sigma2_ah,
      prior_params = prior_params,
      mcmc_params = list(
        n_iter = as.integer(config$n_iter),
        n_burn = as.integer(config$n_burn),
        n_thin = as.integer(config$n_thin),
        seed = as.integer(config$seed)
      ),
      method = "mcmc",
      fold_id = as.integer(fold)
    )

    beta_hat <- colMeans(res$beta_samples, na.rm = TRUE)
    sigma2_e_hat <- mean(res$sigma2_e_samples, na.rm = TRUE)
    sigma2_small_hat <- mean(res$sigma2_small_samples, na.rm = TRUE)
    sigma2_medium_hat <- mean(res$sigma2_medium_samples, na.rm = TRUE)
    sigma2_large_hat <- mean(res$sigma2_large_samples, na.rm = TRUE)
    pi_hat <- colMeans(res$pi_samples, na.rm = TRUE)

    # Mode gamma
    gamma_mode <- apply(res$gamma_samples, 2, function(x) {
      tab <- table(x)
      if(length(tab) == 0) return(0)
      as.numeric(names(sort(tab, decreasing = TRUE)[1]))
    })

    # Diagnostics
    sigma2_e_mcmc <- mcmc(res$sigma2_e_samples)
    ess <- effectiveSize(sigma2_e_mcmc)
    geweke <- geweke.diag(sigma2_e_mcmc)$z
    
  } else {
    stop("config$bayes_algo must be 'mcmc' or 'em'")
  }
  
  # Direct marker effect prediction
  GEBV_train <- W_train %*% beta_hat
  GEBV_test <- W_test %*% beta_hat
  
  # Calculate heritability
  sigma2_g_bayesr <- pi_hat[2] * sigma2_small_hat +
                   pi_hat[3] * sigma2_medium_hat +
                   pi_hat[4] * sigma2_large_hat
  h2_bayesr <- sigma2_g_bayesr / (sigma2_g_bayesr + sigma2_e_hat)

  cat("Posterior Means:\n")
  cat("σ²_e:", round(sigma2_e_hat, 6), "\n")
  cat("σ²_small:", round(sigma2_small_hat, 8), "\n")
  cat("σ²_medium:", round(sigma2_medium_hat, 8), "\n")
  cat("σ²_large:", round(sigma2_large_hat, 8), "\n")
  cat("π:", round(pi_hat, 3), "\n\n")
  cat("σ²_g (genetic):", round(sigma2_g_bayesr, 6), "\n")
  cat("h² (heritability):", round(h2_bayesr, 4), "\n\n")
  
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

  ## Compute Credible Intervals
  # Marker statistics
  if (config$bayes_algo == "mcmc") {
    beta_ci_lower <- apply(res$beta_samples, 2, quantile, probs = 0.025, na.rm = TRUE)
    beta_ci_upper <- apply(res$beta_samples, 2, quantile, probs = 0.975, na.rm = TRUE)
    PIP <- colMeans(res$gamma_samples != 0)
  } else {
    # Untuk EM, tidak ada CI atau PIP
    beta_ci_lower <- rep(NA, n_alleles)
    beta_ci_upper <- rep(NA, n_alleles)
    PIP <- as.numeric(gamma_mode != 0)
  }

  ## Identify Important Alleles
  significant <- (beta_ci_lower > 0 | beta_ci_upper < 0)

  marker_info <- tryCatch({
    # Verify dimensions
    if(nrow(matrices$allele_info) != ncol(W_train)) {
      stop("allele_info rows (", nrow(matrices$allele_info), 
          ") != W_train columns (", ncol(W_train), ")")
    }

    data.frame(
      allele_id = matrices$allele_info$allele_id,
      beta_hat = as.vector(beta_hat),
      beta_abs = abs(as.vector(beta_hat)),
      ci_lower = as.vector(beta_ci_lower),
      ci_upper = as.vector(beta_ci_upper),
      component = as.integer(gamma_mode),
      component_label = factor(as.integer(gamma_mode),
                              levels = 0:3, 
                              labels = c("zero", "small", "medium", "large")),
      PIP = as.vector(PIP),
      significant = as.logical(significant),
      freq = matrices$allele_info$freq,
      block = matrices$allele_info$block,
      allele = matrices$allele_info$allele,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    warning("Marker info for BayesR failed: ", e$message)
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
      sigma2_e = sigma2_e_hat,
      sigma2_small = sigma2_small_hat,
      sigma2_medium = sigma2_medium_hat,
      sigma2_large = sigma2_large_hat,
      pi = pi_hat,
      sigma2_g = sigma2_g_bayesr,
      h2 = h2_bayesr
    ),
    markers = marker_info,
    samples = if(config$save_mcmc_samples && config$bayes_algo == "mcmc") {
      list(beta = res$beta_samples, gamma = res$gamma_samples)
    } else NULL,
    diagnostics = list(
      ess = ess,
      geweke_z = geweke
    ),
    n_cores = 1,
    runtime = NA
  )
}