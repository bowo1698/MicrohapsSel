# config.R h2 = 0.4
# Configuration file for genomic prediction pipeline
iter <- as.integer(Sys.getenv("ITER", "1"))
gen <- as.integer(Sys.getenv("GEN", "1"))

config <- list(
  # Paths
  base_dir = "/scratch/user/aguswibowo/Research/iter_seq/model_run",
  data_dir = file.path("/scratch/user/aguswibowo/Research/simulation/output/trout2", 
                     paste0("iteration_", iter)),
  output_dir = file.path("/scratch/user/aguswibowo/Research/iter_seq/output3/trout_SNP_normal-rand",
                       paste0("iteration_", iter, "_gen", gen)),

  # Data settings
  phenotype = paste0("trout_gen", gen, "_phenotypes_normal.csv"),
  pedigree = "trout_pedigree_all.csv",
  genotype = paste0("genotypes_5k_gen", gen, ".csv"),
  population_type = "reference",
  target_generation = gen,
  
  # QC
  min_allele_freq = 0.01,
  callrate_threshold = 0.8,
  duplicate_threshold = 0.995,

  # Fold setting
  k_folds = 5,
  seed = 123,

  # Model parameters
  bayesian = list(
    niter = 50000,
    nburn = 25000,
    thin = 10
  ),
  bayesR = list(
    fold = c(0, 0.0001, 0.001, 0.01)
  )
)