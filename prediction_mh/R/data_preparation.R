# Data Loading and Preparation Functions

load_genomic_data <- function(config) {
  
  # Load phenotypes
  pheno <- read_csv(
    file.path(config$data_dir, config$phenotype_file),
    show_col_types = FALSE
  )
  
  # Load pedigree
  pedigree <- read_csv(
    file.path(config$data_dir, config$pedigree_file),
    show_col_types = FALSE
  ) %>%
    filter(population_type == config$population_type)
  
  # Map phenotype IDs to pedigree IDs
  pheno <- pheno %>%
    mutate(
      ID = as.character(ID),
      ind_number = as.integer(gsub("reference_gen\\d+_ind", "", ID)),
      pedigree_id = ind_number + config$pedigree_id_offset
    ) %>%
    left_join(
      pedigree %>% select(individual_id, father_id, mother_id),
      by = c("pedigree_id" = "individual_id")
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
    filter(ID %in% pheno_train$ID) %>%
    arrange(match(ID, pheno_train$ID))
  
  hap_test <- data$haplotypes %>%
    filter(ID %in% pheno_test$ID) %>%
    arrange(match(ID, pheno_test$ID))
  
  cat("\n=== Genotype Data Preview ===\n")
  print(head(hap_train[, 1:min(10, ncol(hap_train))]))
  cat("\nDimensions:", dim(hap_train), "\n\n")
  
  list(
    train = list(phenotypes = pheno_train, haplotypes = hap_train),
    test = list(phenotypes = pheno_test, haplotypes = hap_test)
  )
}

compute_allele_frequencies <- function(hap_cols, config) {
  
  allele_freq <- hap_cols %>%
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