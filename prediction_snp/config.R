# config.R
# Configuration file for genomic prediction pipeline
config <- list(
  # Paths
  base_dir = "/scratch/user/aguswibowo/Research/prediction",
  data_dir = file.path("/scratch/user/aguswibowo/Research/simulation/output/trout"),
  output_dir = file.path("/scratch/user/aguswibowo/Research/prediction/output/trout/trout_snp_normal_parallel"),

  # Data settings
  phenotype = "trout_phenotypes_normal.csv",
  pedigree = "trout_pedigree.csv",
  genotype = "trout_genotypes.csv",
  population_type = "reference",
  target_generation = 5,

  # QC
  min_allele_freq = 0.05,
  callrate_threshold = 0.95,
  duplicate_threshold = 0.995,

  # Fold setting
  k_folds = 5,
  seed = 123,

  # Model parameters
  bayesian = list(
    niter = 30000,
    nburn = 10000,
    thin = 5
  ),
  bayesR = list(
    fold = c(0, 0.001, 0.01, 0.1)
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