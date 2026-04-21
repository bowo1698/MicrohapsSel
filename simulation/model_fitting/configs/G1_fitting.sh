#!/bin/bash
#SBATCH --job-name=fit_G1
#SBATCH --array=1-20
#SBATCH --partition=general
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --account=a_zenger
#SBATCH -o logs/fit_G1_%A_%a.out
#SBATCH -e logs/fit_G1_%A_%a.err

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export FLEXIBLAS_NUM_THREADS=1

module load r/4.4.0-gfbf-2023a

ITER=$SLURM_ARRAY_TASK_ID
ITERPAD=$(printf '%02d' $ITER)
MODEL="/scratch/user/aguswibowo/Research/predict_salmon_forward/model_run"
DATA_DIR="/scratch/user/aguswibowo/Research/simulation_forward/output_G1/iter${ITERPAD}"
OUT_BASE="/scratch/user/aguswibowo/Research/predict_salmon_forward/output/fitG1/iter${ITERPAD}"

SCENARIOS='list(
    S1_equal   = "pheno_S1_G1.csv",
    S2_unequal = "pheno_S2_G1.csv",
    S3_equal   = "pheno_S3_G1.csv",
    S4_unequal = "pheno_S4_G1.csv",
    S5_equal   = "pheno_S5_G1.csv",
    S6_unequal = "pheno_S6_G1.csv"
  )'

echo "=== ITER=${ITER} Started: $(date) ==="

# ── MH ──
echo "--- MH started: $(date) ---"
TMPCONFIG=$(mktemp /tmp/config_mh_XXXXXX.R)
cat > $TMPCONFIG << EOF
config <- list(
  base_dir          = "${MODEL}",
  data_dir          = "${DATA_DIR}",
  output_dir        = "${OUT_BASE}/mh",
  block_dir         = "mh_genotypes",
  pedigree_file     = "pedigree_G1.csv",
  population_type   = "reference",
  target_generation = 1,
  n_chromosomes     = 29,
  min_allele_freq   = 0.01,
  max_allele_freq   = 0.99,
  drop_baseline     = TRUE,
  seed              = 123,
  n_iter            = 50000,
  n_burn            = 25000,
  n_thin            = 10,
  bayesR = list(
    pi        = c(0.90, 0.05, 0.03, 0.02),
    var_class = c(0, 0.001, 0.01, 0.1),
    prior_df  = list(residual = 10, genetic = 10)
  ),
  bayesA = list(nu = 4.5, prior_df_residual = 10),
  save_mcmc_samples = FALSE,
  scenarios = ${SCENARIOS}
)
EOF
export GENOMIC_CONFIG=$TMPCONFIG
Rscript ${MODEL}/main-mh-fit.R
rm -f $TMPCONFIG
echo "--- MH done: $(date) ---"

# ── SNP ──
echo "--- SNP started: $(date) ---"
TMPCONFIG=$(mktemp /tmp/config_snp_XXXXXX.R)
cat > $TMPCONFIG << EOF
config <- list(
  base_dir            = "${MODEL}",
  data_dir            = "${DATA_DIR}",
  output_dir          = "${OUT_BASE}/snp",
  genotype            = "geno_G1.csv",
  pedigree            = "pedigree_G1.csv",
  population_type     = "reference",
  target_generation   = 1,
  min_allele_freq     = 0.01,
  callrate_threshold  = 0.8,
  duplicate_threshold = 0.995,
  seed                = 123,
  n_iter            = 50000,
  n_burn            = 25000,
  n_thin            = 10,
  bayesR = list(
    pi        = c(0.90, 0.05, 0.03, 0.02),
    var_class = c(0, 0.001, 0.01, 0.1),
    prior_df  = list(residual = 10, genetic = 10)
  ),
  bayesA = list(nu = 4.5, prior_df_residual = 10),
  save_mcmc_samples = FALSE,
  scenarios = ${SCENARIOS}
)
EOF
export GENOMIC_CONFIG=$TMPCONFIG
Rscript ${MODEL}/main-snp-fit.R
rm -f $TMPCONFIG
echo "--- SNP done: $(date) ---"

echo "=== ITER=${ITER} Completed: $(date) ==="