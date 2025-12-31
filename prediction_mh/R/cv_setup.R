# Cross-Validation Setup Functions

check_family_structure <- function(phenotypes, min_offspring = 5) {
  
  # Check family structure for stratification feasibility
  family_stats <- phenotypes %>%
    count(father_id) %>%
    summarise(
      n_sires = n(),
      min_offspring = min(n),
      median_offspring = median(n),
      max_offspring = max(n),
      sires_with_enough = sum(n >= min_offspring),
      prop_adequate = sires_with_enough / n_sires
    )
  
  list(
    stats = family_stats,
    feasible = family_stats$n_sires > 1 && 
               family_stats$median_offspring >= min_offspring &&
               family_stats$prop_adequate >= 0.8
  )
}

setup_stratified_cv <- function(data, config) {
  
  set.seed(config$cv_seed)
  
  pheno <- data$phenotypes
  k <- config$k_folds
  
  family_check <- check_family_structure(
    pheno, 
    min_offspring = config$min_offspring_per_sire
  )
  
  if(config$stratify_by_sire && family_check$feasible) {
    
    fold_assignments <- createFolds(
      y = pheno$father_id,
      k = k,
      list = FALSE,
      returnTrain = FALSE
    )
    
    stratification_check <- table(pheno$father_id, fold_assignments)
    
    list(
      folds = fold_assignments,
      method = "stratified",
      family_stats = family_check$stats,
      sires_per_fold = colSums(stratification_check > 0),
      avg_folds_per_sire = mean(rowSums(stratification_check > 0))
    )
    
  } else {
    
    fold_assignments <- sample(rep(1:k, length.out = nrow(pheno)))
    
    list(
      folds = fold_assignments,
      method = "random",
      family_stats = family_check$stats,
      reason = if(!family_check$feasible) {
        if(family_check$stats$n_sires <= 1) "only_one_sire"
        else if(family_check$stats$median_offspring < config$min_offspring_per_sire) "insufficient_offspring"
        else "too_few_adequate_sires"
      } else "stratification_disabled"
    )
  }
}

get_fold_split_info <- function(data, cv_setup, fold_id) {
  
  test_idx <- which(cv_setup$folds == fold_id)
  train_idx <- which(cv_setup$folds != fold_id)
  
  list(
    fold = fold_id,
    n_train = length(train_idx),
    n_test = length(test_idx),
    prop_train = length(train_idx) / nrow(data$phenotypes),
    prop_test = length(test_idx) / nrow(data$phenotypes)
  )
}