# W Matrix and Genomic Relationship Matrix Construction
construct_W_matrix <- function(hap_cols, allele_freq_filtered, reference_structure = NULL, drop_baseline = TRUE) {
  
  # Convert to matrix format
  hap_matrix <- as.matrix(hap_cols)
  storage.mode(hap_matrix) <- "integer"
  colnames_vec <- colnames(hap_cols)
  
  # Prepare allele_freq for Rust (training set only)
  allele_freq_list <- if(is.null(reference_structure)) {
    list(
      haplotype = as.character(allele_freq_filtered$haplotype),
      allele = as.integer(allele_freq_filtered$allele),
      freq = as.numeric(allele_freq_filtered$freq)
    )
  } else {
    NULL
  }
  
  # Prepare reference structure for Rust (test set only)
  ref_struct <- if(!is.null(reference_structure)) {
    list(
      allele_info = list(
        allele_id = as.character(reference_structure$allele_info$allele_id),
        freq = as.numeric(reference_structure$allele_info$freq)
      ),
      dropped_alleles = if(!is.null(reference_structure$dropped_alleles) && 
                           nrow(reference_structure$dropped_alleles) > 0) {
        list(
          block = as.character(reference_structure$dropped_alleles$block),
          allele = as.integer(reference_structure$dropped_alleles$allele),
          freq = as.numeric(reference_structure$dropped_alleles$freq)
        )
      } else {
        list()
      }
    )
  } else {
    NULL
  }
  
  # Call Rust implementation
  result <- masbayes::construct_wah_matrix(
    hap_matrix = hap_matrix,
    colnames = colnames_vec,
    allele_freq_filtered = allele_freq_list,
    reference_structure = ref_struct,
    drop_baseline = drop_baseline
  )
  
  # Ensure column names are set
  if(is.null(colnames(result$W_ah)) && 
     !is.null(result$allele_info) && 
     "allele_id" %in% names(result$allele_info)) {
    colnames(result$W_ah) <- result$allele_info$allele_id
  }
  
  # Print summary
  cat("W_ah final columns:", ncol(result$W_ah), "\n")
  cat("Metadata rows:", length(result$allele_info$allele_id), "\n")
  cat("Dropped alleles:", length(result$dropped_alleles$block), "\n")
  
  result
}

construct_grm <- function(W_train, W_test = NULL) {

  n_train <- nrow(W_train)

  # Calculate k_αh from reference population
  G_train <- W_train %*% t(W_train)
  k_ah <- sum(diag(G_train)) / n_train

  # Construct A_gh for reference
  A_train <- G_train / k_ah

  cat("A_gh_ref dimension:", dim(A_train), "\n\n")
  cat("Mean diagonal:", mean(diag(A_train)), "\n")
  cat("Mean off-diagonal:", mean(A_train[upper.tri(A_train)]), "\n\n")

  if(is.null(W_test)) {
    return(list(
      A_train = A_train,
      A_combined = A_train,
      k_ah = k_ah
    ))
  }

  # Construct combined relationship matrix
  W_combined <- rbind(W_train, W_test)
  G_combined <- W_combined %*% t(W_combined)
  A_combined <- G_combined / k_ah
  
  cat("A_gh_combined dimension:", dim(A_combined), "\n")
  cat("Mean diagonal:", mean(diag(A_combined)), "\n")
  cat("Mean off-diagonal:", mean(A_combined[upper.tri(A_combined)]), "\n\n")

  list(
    A_train = A_train,
    A_combined = A_combined,
    k_ah = k_ah
  )
}

construct_matrices <- function(train_data, test_data, allele_freq, config) {

  # Extract haplotype columns only (exclude ID)
  hap_train <- train_data$haplotypes %>% select(-individual_id)
  hap_test <- test_data$haplotypes %>% select(-individual_id)
  
  # Construct W matrix for reference and test
  W_result_train <- construct_W_matrix(
    hap_train, 
    allele_freq, 
    drop_baseline = config$drop_baseline
  )

  cat("\n=== W_ah Matrix Preview (Training) ===\n")
  cat("Dimensions:", dim(W_result_train$W_ah), "\n")
  print(head(W_result_train$W_ah[, 1:min(10, ncol(W_result_train$W_ah))]))
  cat("\nFirst 10 column names:\n")
  print(colnames(W_result_train$W_ah)[1:10])
  cat("\n")
  
  W_result_test <- construct_W_matrix(
    hap_test, 
    allele_freq_filtered = NULL,
    reference_structure = W_result_train,
    drop_baseline = config$drop_baseline
  )
  
  grm <- construct_grm(W_result_train$W_ah, W_result_test$W_ah)

  all_ids <- c(as.character(train_data$phenotypes$individual_id), 
             as.character(test_data$phenotypes$individual_id))
  rownames(grm$A_combined) <- all_ids
  colnames(grm$A_combined) <- all_ids
  
  list(
    W_train = W_result_train$W_ah,
    W_test = W_result_test$W_ah,
    allele_info = W_result_train$allele_info,
    dropped_alleles = W_result_train$dropped_alleles,
    A_train = grm$A_train,
    A_combined = grm$A_combined,
    k_ah = grm$k_ah,
    WtW_diag = colSums(W_result_train$W_ah^2),
    Wty = crossprod(W_result_train$W_ah, train_data$phenotypes$phenotype)
  )
}