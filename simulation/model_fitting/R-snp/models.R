# models.R
# Model training and prediction functions
train_gblup <- function(y_train, GRM_train, GRM_val = NULL, y_val = NULL, val_ids = NULL, n_threads) {
  result <- tryCatch({
    train_time <- system.time({
      n <- length(y_train)
      if(is.null(rownames(GRM_train))) {
        id_names <- paste0("ID_", 1:n)
        rownames(GRM_train) <- id_names
        colnames(GRM_train) <- id_names
      } else {
        id_names <- rownames(GRM_train)
      }
      
      pheno_df <- data.frame(
        ID = id_names,
        y = y_train,
        units = 1:n
      )
      
      fit <- mmes(
        fixed = y ~ 1,
        random = ~ vsm(ism(ID), Gu = GRM_train),
        rcov = ~ vsm(ism(units)),
        data = pheno_df,
        verbose = FALSE
      )
      
      var_comp <- fit$theta
      var_g <- var_comp[[1]][1,1]
      var_e <- var_comp[[2]][1,1]
      h2 <- var_g / (var_g + var_e)

      h2_test <- NA
      if (!is.null(GRM_val) && !is.null(y_val)) {
        n_val <- length(y_val)
        v_ids <- if (!is.null(val_ids)) as.character(val_ids) else paste0("ID_", 1:n_val)
        rownames(GRM_val) <- v_ids
        colnames(GRM_val) <- v_ids
        pheno_val_df <- data.frame(ID = v_ids, y = y_val, units = 1:n_val)
        h2_test <- tryCatch({
          fit_val <- mmes(
            fixed = y ~ 1,
            random = ~ vsm(ism(ID), Gu = GRM_val),
            rcov = ~ vsm(ism(units)),
            data = pheno_val_df,
            verbose = FALSE
          )
          vc <- fit_val$theta
          sg <- vc[[1]][1,1]; se <- vc[[2]][1,1]
          h2_raw <- sg / (sg + se)
          if (is.na(h2_raw) || h2_raw < 0 || h2_raw > 1) NA else h2_raw
        }, error = function(e) NA)
      }
      
      # Extract BLUP - handle list or matrix
      if(is.list(fit$u)) {
        blup <- fit$u[[1]]
      } else {
        blup <- fit$u
      }
      if(is.null(dim(blup))) blup <- matrix(blup, ncol = 1)
      
      # Extract IDs from BLUP
      if(!is.null(rownames(blup))) {
        blup_ids <- rownames(blup)
      } else if(!is.null(dimnames(blup)[[1]])) {
        blup_ids <- dimnames(blup)[[1]]
      } else {
        blup_ids <- id_names
      }
      
      # Match and reorder to phenotype ID order
      matched_idx <- match(id_names, blup_ids)
      if(any(is.na(matched_idx))) {
        stop("ID mismatch between id_names and BLUP IDs")
      }
      u <- as.vector(blup[matched_idx])
    })
    
    list(
      success = TRUE,
      fit = fit,
      u = u,
      id_names = id_names,
      var_g = var_g,
      var_e = var_e,
      train_time = train_time[3],
      h2 = h2,
      h2_test = h2_test
    )
  }, error = function(e) {
    cat("GBLUP FAILED:", e$message, "\n")
    list(success = FALSE, train_time = NA, h2 = NA)
  })
  
  return(result)
}

predict_gblup_train <- function(gblup_fit, GRM_train, y_train) {
  result <- tryCatch({
    fit <- gblup_fit$fit
    id_names <- gblup_fit$id_names
    
    mu <- as.numeric(fit$b)
    
    # Extract BLUP - handle list or matrix
    if(is.list(fit$u)) {
      blup <- fit$u[[1]]
    } else {
      blup <- fit$u
    }
    if(is.null(dim(blup))) blup <- matrix(blup, ncol = 1)
    
    # Extract IDs
    if(!is.null(rownames(blup))) {
      blup_ids <- rownames(blup)
    } else if(!is.null(dimnames(blup)[[1]])) {
      blup_ids <- dimnames(blup)[[1]]
    } else {
      blup_ids <- id_names
    }
    
    # Match and reorder
    matched_idx <- match(id_names, blup_ids)
    if(any(is.na(matched_idx))) {
      stop("ID mismatch in predict_gblup_train")
    }
    u <- as.vector(blup[matched_idx])
    
    pred <- mu + u
    list(pred = pred)
  }, error = function(e) {
    cat("ERROR in predict_gblup_train:", e$message, "\n")
    stop(e)
  })
  
  return(result)
}

predict_gblup <- function(gblup_fit, GRM_val_train, y_train, GRM_train) {
  result <- tryCatch({
    pred_time <- system.time({
      fit <- gblup_fit$fit
      id_names <- gblup_fit$id_names
      
      mu <- as.numeric(fit$b)
      
      # Extract BLUP - handle list or matrix
      if(is.list(fit$u)) {
        blup <- fit$u[[1]]
      } else {
        blup <- fit$u
      }
      if(is.null(dim(blup))) blup <- matrix(blup, ncol = 1)
      
      # Extract IDs
      if(!is.null(rownames(blup))) {
        blup_ids <- rownames(blup)
      } else if(!is.null(dimnames(blup)[[1]])) {
        blup_ids <- dimnames(blup)[[1]]
      } else {
        blup_ids <- id_names
      }
      
      # Match and reorder
      matched_idx <- match(id_names, blup_ids)
      if(any(is.na(matched_idx))) {
        stop("ID mismatch in predict_gblup")
      }
      u_train <- as.vector(blup[matched_idx])
      
      # Prediction with lambda
      var_comp <- fit$theta
      var_g <- var_comp[[1]][1,1]
      var_e <- var_comp[[2]][1,1]
      lambda <- var_e / var_g
      
      A_ref_inv <- solve(GRM_train + diag(lambda, nrow(GRM_train)))
      pred <- mu + drop(GRM_val_train %*% A_ref_inv %*% u_train)
    })
    
    list(success = TRUE, pred = pred, pred_time = pred_time[3])
  }, error = function(e) {
    cat("ERROR in predict_gblup:", e$message, "\n")
    list(success = FALSE, pred = NA, pred_time = NA)
  })
  
  return(result)
}

train_bayesA <- function(W_train, y, sigma2_g, sigma2_e, config) {
  result <- tryCatch({
    train_time <- system.time({
      wtw <- colSums(W_train^2)
      wty <- as.vector(crossprod(W_train, y))
      res <- masbayes::run_bayesa(
        w             = W_train,
        y             = y,
        wtw_diag      = wtw,
        wty           = wty,
        nu            = config$bayesA$nu,
        sigma2_g      = sigma2_g,
        sigma2_e_init = sigma2_e,
        prior_params  = list(a0_e = config$bayesA$prior_df_residual),
        mcmc_params   = list(
          n_iter = as.integer(config$n_iter),
          n_burn = as.integer(config$n_burn),
          n_thin = as.integer(config$n_thin),
          seed   = as.integer(config$seed)
        ),
        method        = "mcmc",
        fold_id       = 0L
      )
    })

    beta_hat      <- res$beta_hat
    mu_hat        <- res$mu_hat
    sigma2_e_hat  <- res$sigma2_e_hat
    sigma2_g_hat  <- res$sigma2_g
    h2            <- res$h2
    pred_train    <- as.vector(W_train %*% beta_hat) + mu_hat

    beta_ci_lower <- apply(res$beta_samples, 2, quantile, 0.025)
    beta_ci_upper <- apply(res$beta_samples, 2, quantile, 0.975)
    significant   <- (beta_ci_lower > 0 | beta_ci_upper < 0)

    marker_info <- data.frame(
      beta_hat    = as.vector(beta_hat),
      beta_abs    = abs(as.vector(beta_hat)),
      ci_lower    = as.vector(beta_ci_lower),
      ci_upper    = as.vector(beta_ci_upper),
      sigma2_j    = as.vector(res$sigma2_j_hat),
      significant = as.logical(significant),
      stringsAsFactors = FALSE
    )

    list(
      success      = TRUE,
      beta_hat     = beta_hat,
      mu_hat       = mu_hat,
      sigma2_e_hat = sigma2_e_hat,
      sigma2_g_hat = sigma2_g_hat,
      h2           = h2,
      pred_train   = pred_train,
      markers      = marker_info,
      samples      = if (isTRUE(config$save_mcmc_samples))
                       list(beta = res$beta_samples) else NULL,
      train_time   = train_time[3]
    )
  }, error = function(e) {
    cat("BayesA FAILED:", e$message, "\n")
    list(success = FALSE, train_time = NA, h2 = NA)
  })
  return(result)
}

train_bayesR <- function(W_train, y, sigma2_g, sigma2_e, config) {
  result <- tryCatch({
    train_time <- system.time({
      wtw <- colSums(W_train^2)
      wty <- as.vector(crossprod(W_train, y))
      res <- masbayes::run_bayesr(
        w             = W_train,
        y             = y,
        wtw_diag      = wtw,
        wty           = wty,
        pi_vec        = config$bayesR$pi,
        sigma2_e_init = sigma2_e,
        sigma2_ah     = sigma2_g,
        prior_params  = list(
          a0_e           = config$bayesR$prior_df$residual,
          a0_g           = config$bayesR$prior_df$genetic,
          variance_class = config$bayesR$var_class
        ),
        mcmc_params   = list(
          n_iter = as.integer(config$n_iter),
          n_burn = as.integer(config$n_burn),
          n_thin = as.integer(config$n_thin),
          seed   = as.integer(config$seed)
        ),
        method        = "mcmc",
        fold_id       = 0L
      )
    })

    beta_hat      <- res$beta_hat
    mu_hat        <- res$mu_hat
    sigma2_e_hat  <- res$sigma2_e_hat
    sigma2_g_hat  <- res$sigma2_g
    h2            <- res$h2
    pred_train    <- as.vector(W_train %*% beta_hat) + mu_hat

    beta_ci_lower <- apply(res$beta_samples, 2, quantile, 0.025)
    beta_ci_upper <- apply(res$beta_samples, 2, quantile, 0.975)
    PIP           <- colMeans(res$gamma_samples != 0)
    gamma_mode    <- apply(res$gamma_samples, 2, function(x) {
      tab <- table(x); as.numeric(names(sort(tab, decreasing=TRUE)[1]))
    })
    significant   <- (beta_ci_lower > 0 | beta_ci_upper < 0)

    marker_info <- data.frame(
      beta_hat        = as.vector(beta_hat),
      beta_abs        = abs(as.vector(beta_hat)),
      ci_lower        = as.vector(beta_ci_lower),
      ci_upper        = as.vector(beta_ci_upper),
      component       = as.integer(gamma_mode),
      PIP             = as.vector(PIP),
      significant     = as.logical(significant),
      stringsAsFactors = FALSE
    )

    list(
      success      = TRUE,
      beta_hat     = beta_hat,
      mu_hat       = mu_hat,
      sigma2_e_hat = sigma2_e_hat,
      sigma2_g_hat = sigma2_g_hat,
      h2           = h2,
      pred_train   = pred_train,
      markers      = marker_info,
      samples      = if (isTRUE(config$save_mcmc_samples))
                       list(beta = res$beta_samples, gamma = res$gamma_samples) else NULL,
      train_time   = train_time[3]
    )
  }, error = function(e) {
    cat("BayesR FAILED:", e$message, "\n")
    list(success = FALSE, train_time = NA, h2 = NA)
  })
  return(result)
}

predict_bayes_marker <- function(W_new, beta_hat, mu_hat) {
  as.vector(W_new %*% beta_hat) + mu_hat
}