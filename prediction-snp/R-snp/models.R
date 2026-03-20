# models.R
# Model training and prediction functions

suppressPackageStartupMessages({
  library(sommer)
  library(hibayes)
})

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

train_bayesA <- function(pheno_data, geno_train_filtered, train_ids, config, n_threads = 1, y_val = NULL, geno_val = NULL) {
  result <- tryCatch({
    train_time <- system.time({
      fit <- ibrm(
        formula = phenotype ~ 1,
        data = pheno_data,
        M = geno_train_filtered,
        M.id = train_ids,
        method = "BayesA",
        niter = config$bayesian$niter,
        nburn = config$bayesian$nburn,
        thin = config$bayesian$thin,
        threads = n_threads,
        verbose = FALSE
      )
    })

    h2_test <- NA
    if (!is.null(y_val) && !is.null(geno_val)) {
      pred_val <- as.vector(fit$mu + geno_val %*% fit$alpha)
      ve_post <- as.numeric(fit$Ve)
      h2_test <- tryCatch({
        h2_raw <- var(pred_val) / (var(pred_val) + ve_post)
        if (!is.finite(h2_raw) || h2_raw < 0 || h2_raw > 1) NA_real_ else h2_raw
      }, error = function(e) NA_real_)
    }
    
    list(
      success = TRUE,
      fit = fit,
      train_time = train_time[3],
      h2 = fit$h2,
      h2_test = h2_test
    )
  }, error = function(e) {
    cat("BayesA FAILED:", e$message, "\n")
    list(success = FALSE, train_time = NA, h2 = NA)
  })
  
  return(result)
}

predict_bayesA_train <- function(fit, geno_train_filtered) {
  mu <- as.numeric(fit$mu)
  alpha <- as.numeric(fit$alpha)
  pred <- mu + as.vector(geno_train_filtered %*% alpha)
  list(pred = pred)
}

predict_bayesA <- function(fit, geno_val_filtered) {
  result <- tryCatch({
    pred_time <- system.time({
      pred <- as.vector(fit$mu + geno_val_filtered %*% fit$alpha)
    })
    
    list(success = TRUE, pred = pred, pred_time = pred_time[3])
  }, error = function(e) {
    cat("ERROR in predict_bayesA:", e$message, "\n")
    list(success = FALSE, pred = NA, pred_time = NA)
  })
  
  return(result)
}

train_bayesR <- function(pheno_data, geno_train_filtered, train_ids, vg, ve, config, n_threads = 1, y_val = NULL, geno_val = NULL) {
  result <- tryCatch({
    train_time <- system.time({
      fit <- ibrm(
        formula = phenotype ~ 1,
        data = pheno_data,
        M = geno_train_filtered,
        M.id = train_ids,
        method = "BayesR",
        fold = config$bayesR$var_class,
        vg = vg,
        ve = ve,
        niter = config$bayesian$niter,
        nburn = config$bayesian$nburn,
        thin = config$bayesian$thin,
        threads = n_threads,
        verbose = FALSE
      )
    })

    h2_test <- NA
    if (!is.null(y_val) && !is.null(geno_val)) {
      pred_val <- as.vector(fit$mu + geno_val %*% fit$alpha)
      ve_post <- as.numeric(fit$Ve)
      h2_test <- tryCatch({
        h2_raw <- var(pred_val) / (var(pred_val) + ve_post)
        if (!is.finite(h2_raw) || h2_raw < 0 || h2_raw > 1) NA_real_ else h2_raw
      }, error = function(e) NA_real_)
    }
    
    list(
      success = TRUE,
      fit = fit,
      train_time = train_time[3],
      h2 = fit$h2,
      h2_test = h2_test
    )
  }, error = function(e) {
    cat("BayesR FAILED:", e$message, "\n")
    list(success = FALSE, train_time = NA, h2 = NA)
  })
  
  return(result)
}

predict_bayesR_train <- function(fit, geno_train_filtered) {
  mu <- as.numeric(fit$mu)
  alpha <- as.numeric(fit$alpha)
  pred <- mu + as.vector(geno_train_filtered %*% alpha)
  list(pred = pred)
}

predict_bayesR <- function(fit, geno_val_filtered) {
  result <- tryCatch({
    pred_time <- system.time({
      pred <- as.vector(fit$mu + geno_val_filtered %*% fit$alpha)
    })
    
    list(success = TRUE, pred = pred, pred_time = pred_time[3])
  }, error = function(e) {
    list(success = FALSE, pred = NA, pred_time = NA)
  })
  
  return(result)
}