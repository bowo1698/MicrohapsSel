# Data Loading and Preparation Functions

load_genomic_data <- function(config) {
  
  # Load phenotypes
  pheno <- read_csv(
    file.path(config$data_dir, config$phenotype_file),
    show_col_types = FALSE
  ) %>%
    filter(generation == config$target_generation)

  pedigree <- read_csv(
    file.path(config$data_dir, config$pedigree_file),
    show_col_types = FALSE
  ) %>%
    filter(
      population_type == config$population_type,
      generation == config$target_generation
    )

  pheno <- pheno %>%
    left_join(
      pedigree %>% select(individual_id, father_id, mother_id),
      by = "individual_id"
    )
  
  # Load haplotype genotypes
  hap_geno_list <- lapply(1:config$n_chromosomes, function(chr) {
    read_table(
      file.path(config$data_dir, config$block_dir, paste0("hap_geno_", chr)),
      show_col_types = FALSE
    )
  })
  
  # Merge chromosomes
  hap_geno <- hap_geno_list[[1]]
  for(chr in 2:config$n_chromosomes) {
    hap_geno <- bind_cols(
      hap_geno,
      hap_geno_list[[chr]] %>% select(-ID)
    )
  }

  hap_geno <- hap_geno %>% rename(individual_id = ID)
  
  common_ids <- intersect(pheno$individual_id, hap_geno$individual_id)
  if (length(common_ids) < nrow(pheno)) {
    cat("Individuals in phenotype but not genotype:", nrow(pheno) - length(common_ids), "\n")
  }
  if (length(common_ids) < nrow(hap_geno)) {
    cat("Individuals in genotype but not phenotype:", nrow(hap_geno) - length(common_ids), "\n")
  }
  pheno_ids <- pheno$individual_id[pheno$individual_id %in% common_ids]
  pheno    <- pheno    %>% filter(individual_id %in% common_ids) %>% arrange(match(individual_id, pheno_ids))
  hap_geno <- hap_geno %>% filter(individual_id %in% common_ids) %>% arrange(match(individual_id, pheno_ids))
  cat("Common individuals used:", length(common_ids), "\n\n")

  cat("Haplotypes after filter:", nrow(hap_geno), "\n")
  cat("Order matches:", identical(hap_geno$individual_id, pheno$individual_id), "\n\n")
  
  list(
    phenotypes = pheno,
    haplotypes = hap_geno,
    pedigree = pedigree
  )
}

split_fold_data <- function(data, fold_assignments, fold_id) {
  
  test_idx <- which(fold_assignments == fold_id)
  train_idx <- which(fold_assignments != fold_id)
  
  pheno_train <- data$phenotypes[train_idx, ]
  pheno_test <- data$phenotypes[test_idx, ]
  
  hap_train <- data$haplotypes %>%
    filter(individual_id %in% pheno_train$individual_id) %>%
    arrange(match(individual_id, pheno_train$individual_id))

  hap_test <- data$haplotypes %>%
    filter(individual_id %in% pheno_test$individual_id) %>%
    arrange(match(individual_id, pheno_test$individual_id))
  
  cat("\n=== Genotype Data Preview ===\n")
  print(head(hap_train[, 1:min(10, ncol(hap_train))]))
  cat("\nDimensions:", dim(hap_train), "\n\n")
  
  list(
    train = list(phenotypes = pheno_train, haplotypes = hap_train),
    test = list(phenotypes = pheno_test, haplotypes = hap_test)
  )
}

compute_allele_frequencies <- function(hap_cols, config) {

  # Pool both copy per blok with strip suffix _1
  col_names <- colnames(hap_cols)
  block_names <- sub("_1$", "", col_names)
  hap_cols_pooled <- hap_cols
  colnames(hap_cols_pooled) <- block_names
  
  allele_freq <- hap_cols_pooled %>%
    pivot_longer(everything(), names_to = "haplotype", values_to = "allele") %>%
    count(haplotype, allele) %>%
    group_by(haplotype) %>%
    mutate(freq = n / sum(n)) %>%
    select(-n) 
  
  allele_freq_filtered <- allele_freq %>%
    filter(freq >= config$min_allele_freq & freq <= config$max_allele_freq) %>%
    ungroup()

  cat("Total alleles before filter:", nrow(allele_freq), "\n")
  cat("Total alleles after filter:", nrow(allele_freq_filtered), "\n")
  cat("Alleles removed:", nrow(allele_freq) - nrow(allele_freq_filtered), "\n\n")
  
  allele_freq_filtered
}