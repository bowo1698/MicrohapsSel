# Global Configuration for Genomic Prediction Pipeline

iter <- as.integer(Sys.getenv("ITER", "1"))
gen <- as.integer(Sys.getenv("GEN", "1"))

config <- list(
  # Paths
  base_dir = "/scratch/user/aguswibowo/Research/iter_seq/model_run",
  data_dir = file.path("/scratch/user/aguswibowo/Research/simulation/output/trout", ##
                     paste0("iteration_", iter)),
  output_dir = file.path("/scratch/user/aguswibowo/Research/iter_seq/output/trout_MH_normal", 
                     paste0("iteration_", iter, "_gen", gen)),
  
  # Data settings
  block_dir = paste0("mh_genotypes_gen", gen),
  phenotype_file = paste0("trout_gen", gen, "_phenotypes_normal.csv"),
  pedigree_file = "trout_pedigree_all.csv",
  population_type = "reference",
  target_generation = gen,
  n_chromosomes = 29,
  
  # W matrix construction
  min_allele_freq = 0.01,
  max_allele_freq = 0.99,
  drop_baseline = TRUE,
  
  # Cross-validation
  k_folds = 10,
  cv_seed = 123,
  stratify_by_sire = TRUE,
  min_offspring_per_sire = 5,
  
  # Algorithm selection
  bayes_algo = "mcmc",  # "mcmc" or "em"
  
  # MCMC settings
  n_iter = 40000,      # MCMC iterations or EM max_iter
  n_burn = 20000,      # for MCMC
  n_thin = 10,         # for MCMC
  seed = 123,
  em_tol = 1e-6,       # for EM
  
  # BayesR hyperparameters
  bayesR = list(
    pi = c(0.5, 0.487, 0.01, 0.003),  # Zero, small, medium, large
    variance_scaling = c(0, 0.0001, 0.001, 0.01), # Relative to sigma2_ah
    prior_df = list(
      residual = 10,
      small = 10,
      medium = 5,
      large = 5
    )
  ),
  
  # BayesA hyperparameters
  bayesA = list(
    nu = 4.5,  # degrees of freedom
    prior_df_residual = 10
  ),

  save_mcmc_samples = FALSE
)