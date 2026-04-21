#!/bin/bash
#SBATCH --job-name=pred_snp_G2
#SBATCH --array=1-40
#SBATCH --partition=general
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --time=120:00:00
#SBATCH --account=a_zenger
#SBATCH -o logs/pred_snp_G2_%A_%a.out
#SBATCH -e logs/pred_snp_G2_%A_%a.err

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export FLEXIBLAS_NUM_THREADS=1

module load r/4.4.0-gfbf-2023a

MODEL="/scratch/user/aguswibowo/Research/predict_salmon_forward/model_run"
QRIS_TAR_DIR="/QRISdata/Q8514/research_data/simulation_forward/output_G2_tar"
OUT_BASE="/scratch/user/aguswibowo/Research/predict_salmon_forward/output/predict_within_gen_G2_full"

SOURCE_COMBOS=(
  "mh_gblup_S1_equal"
  "snp_gblup_S1_equal"
)

SCENARIOS_G2='list(
  S1_equal   = "pheno_S1_equal_G2.csv",
  S2_unequal = "pheno_S2_unequal_G2.csv",
  S3_equal   = "pheno_S3_equal_G2.csv",
  S4_unequal = "pheno_S4_unequal_G2.csv",
  S5_equal   = "pheno_S5_equal_G2.csv",
  S6_unequal = "pheno_S6_unequal_G2.csv"
)'

TASK_ID=$SLURM_ARRAY_TASK_ID
ITER=$(( (TASK_ID - 1) / 2 + 1 ))
COMBO_IDX=$(( (TASK_ID - 1) % 2 ))
ITERPAD=$(printf '%02d' $ITER)

SOURCE_COMBO="${SOURCE_COMBOS[$COMBO_IDX]}"
QRIS_TAR="${QRIS_TAR_DIR}/iter${ITERPAD}_${SOURCE_COMBO}.tar.gz"
OUT_DIR="${OUT_BASE}/iter${ITERPAD}/${SOURCE_COMBO}"
G2_DIR="${TMPDIR}/iter${ITERPAD}_${SOURCE_COMBO}_snp"

echo "=== TASK=${TASK_ID} | ITER=${ITER} | SOURCE=${SOURCE_COMBO} | SNP pred started: $(date) ==="

if [ ! -f "${QRIS_TAR}" ]; then
  echo "SKIP — tar not found: ${QRIS_TAR}"
  exit 0
fi

if [ -f "${OUT_DIR}/snp/S6_unequal/bayesR/predict_cv.rds" ]; then
  echo "SKIP — already done: ${OUT_DIR}/snp"
  exit 0
fi

mkdir -p ${G2_DIR}
tar -xzf ${QRIS_TAR} -C ${G2_DIR}

TMPCONFIG=$(mktemp /tmp/config_snp_XXXXXX.R)
cat > $TMPCONFIG << EOF
config <- list(
  base_dir            = "${MODEL}",
  data_dirs           = list(G2 = "${G2_DIR}"),
  output_dir          = "${OUT_DIR}/snp",
  genotype_by_gen     = list(G2 = "geno_G2.csv"),
  pedigree            = "pedigree_G2.csv",
  population_type     = "candidate",
  min_allele_freq     = 0.01,
  callrate_threshold  = 0.8,
  duplicate_threshold = 0.995,
  seed_split          = $((ITER * 42)),
  train_prop          = 0.8,
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
  scenarios_by_gen  = list(G2 = ${SCENARIOS_G2})
)
EOF

export GENOMIC_CONFIG=$TMPCONFIG
Rscript ${MODEL}/main-snp-fit-cv-eval.R
rm -f $TMPCONFIG
rm -rf ${G2_DIR}

echo "=== TASK=${TASK_ID} | ITER=${ITER} | SOURCE=${SOURCE_COMBO} | SNP pred done: $(date) ==="