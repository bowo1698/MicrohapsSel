# Global Configuration for Genomic Prediction Pipeline

config <- list(
  # Paths
  base_dir = "/scratch/user/aguswibowo/Research/prediction",
  data_dir = file.path("/QRISdata/Q8514/research_data/trout"),
  output_dir = file.path("/scratch/user/aguswibowo/Research/prediction/output/trout/haplo_normal"), # haplo_mixture | mh_normal | mh_mixture # nolint
  
  # Data settings
  block_dir = "mh_genotypes_haploblock", # mh_genotypes_haploblock_125bp | haplotype_genotypes # nolint
  phenotype_file = "trout_phenotypes.csv", # trout_phenotypes_mixture.csv
  pedigree_file = "trout_pedigree.csv",
  population_type = "reference",
  n_chromosomes = 25,
  pedigree_id_offset = 3469,
  
  # W matrix construction
  min_allele_freq = 0.01,
  max_allele_freq = 0.99,
  drop_baseline = TRUE,
  
  # Cross-validation
  k_folds = 5,
  cv_seed = 123,
  stratify_by_sire = TRUE,
  min_offspring_per_sire = 5,
  
  # MCMC settings
  mcmc = list(
    n_iter = 50000,
    n_burn = 20000,
    n_thin = 10,
    seed = 123
  ),
  
  # BayesR hyperparameters
  bayesR = list(
    pi = c(0.9, 0.05, 0.03, 0.02),  # Zero, small, medium, large
    variance_scaling = c(0, 0.01, 0.1, 1),  # Relative to sigma2_ah
    prior_df = list(
      residual = 10,
      small = 5,
      medium = 5,
      large = 5
    )
  ),
  
  # BayesA hyperparameters
  bayesA = list(
    nu = 4.5,  # degrees of freedom
    prior_df_residual = 10
  ),
  
  # XGBoost hyperparameters
  xgboost = list(
    max_depth = 3,
    eta = 0.03,
    lambda = 10,
    colsample_bytree = 0.5,
    subsample = 0.8,
    nrounds = 500
  ),

  save_mcmc_samples = FALSE
)