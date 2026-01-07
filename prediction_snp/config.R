# config.R h2 = 0.4
# Configuration file for genomic prediction pipeline
config <- list(
  # Paths
  base_dir = "/scratch/user/aguswibowo/Research/iter_seq/snp_run",
  data_dir = file.path("/scratch/user/aguswibowo/Research/simulation/output/trout_gens/iteration_1"), ##
  output_dir = file.path("/scratch/user/aguswibowo/Research/iter_seq/output/trout/comparison-snp"), ##

  # Data settings
  phenotype = "trout_gen1_phenotypes_normal.csv",
  pedigree = "trout_pedigree_all.csv",
  genotype = "trout_gen1_genotypes.csv",
  population_type = "reference",
  target_generation = 1,
  
  # QC
  min_allele_freq = 0.01,
  callrate_threshold = 0.8,
  duplicate_threshold = 0.995,

  # Fold setting
  k_folds = 10,
  seed = 123,

  # Model parameters
  bayesian = list(
    niter = 50000,
    nburn = 20000,
    thin = 10
  ),
  bayesR = list(
    fold = c(0, 0.0001, 0.001, 0.01)
  ),
  xgb = list(
    eta = 0.1,
    max_depth = 4,
    subsample = 1,
    colsample_bytree = 0.6,
    lambda = 1,
    nrounds = 200
  )
)