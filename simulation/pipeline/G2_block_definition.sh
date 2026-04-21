#!/bin/bash
#SBATCH --job-name=mh_G2
#SBATCH --partition=general
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=200G
#SBATCH --time=60:00:00
#SBATCH --account=a_zenger
#SBATCH --array=1-20
#SBATCH -o logs/mh_G2_%A_%a.out
#SBATCH -e logs/mh_G2_%A_%a.err

set -uo pipefail

module load anaconda3/2023.09-0
source activate microhaplotype
module load beagle/5.4.22jul22.46e-java-11

cd /scratch/user/aguswibowo/Research/simulation_forward

BASE=/scratch/user/aguswibowo/Research/simulation_forward
SCRATCH_TMP=${BASE}/tmp_mh_G2
QRIS_TAR_DIR=/QRISdata/Q8514/research_data/simulation_forward/output_G2_tar
NCORES=20
MAP=${BASE}/raw_data/salmon_60k_map_final.txt
REF_VCF=${BASE}/raw_data/salmon_60k_pairs125bp.vcf.gz

ITER=$(printf "%02d" ${SLURM_ARRAY_TASK_ID})

MARKER_TYPES=("snp" "mh")
MODELS=("gblup" "bayesR" "bayesA")
SCENARIOS=("S1_equal" "S2_unequal" "S3_equal" "S4_unequal" "S5_equal" "S6_unequal")

run_with_retry() {
    local max_attempts=3
    local delay=60
    for attempt in $(seq 1 $max_attempts); do
        "$@" && return 0
        echo "Attempt ${attempt}/${max_attempts} failed. Retrying in ${delay}s..."
        sleep $delay
    done
    echo "All ${max_attempts} attempts failed."
    return 1
}

echo "=== MH pipeline G2 iter=${ITER} started: $(date) ==="

for marker in "${MARKER_TYPES[@]}"; do
  for model in "${MODELS[@]}"; do
    for sc in "${SCENARIOS[@]}"; do

      COMBO="${marker}_${model}_${sc}"
      QRIS_TAR=${QRIS_TAR_DIR}/iter${ITER}_${COMBO}.tar.gz
      SCRATCH_COMBO=${SCRATCH_TMP}/iter${ITER}/${COMBO}

      # Cek tar ada di QRIS
      if [ ! -f "${QRIS_TAR}" ]; then
        echo "SKIP ${COMBO} — tar not found in QRIS"
        continue
      fi

      # Skip jika MH pipeline sudah selesai
      MH_DONE=$(tar -tzf ${QRIS_TAR} 2>/dev/null | grep "mh_genotypes/" | wc -l)
      if [ "${MH_DONE}" -gt 0 ]; then
        echo "SKIP ${COMBO} — MH pipeline already done"
        continue
      fi

      echo "--- Processing ${COMBO}: $(date) ---"

      # Setup scratch
      rm -rf ${SCRATCH_COMBO}
      mkdir -p ${SCRATCH_COMBO}
      
      tar -xzf ${QRIS_TAR} -C ${SCRATCH_COMBO}

      # Step 2a: convert-to-vcf
      echo "  Step 2a: convert-to-vcf"
      ./convert-to-vcf \
          -i ${SCRATCH_COMBO}/geno_G2.csv \
          -m ${MAP} \
          -o ${SCRATCH_COMBO}/geno_G2.vcf \
          --split-by-chr --parallel --ncores ${NCORES} --merge-after-split

      # Step 2a.5: remove duplicate positions
      echo "  Step 2a.5: remove duplicates"
      CONTIG_HEADER=${SCRATCH_COMBO}/contig_header.txt
      bcftools view -h ${REF_VCF} | grep "^##contig" > ${CONTIG_HEADER}

      bcftools annotate \
          --header-lines ${CONTIG_HEADER} \
          ${SCRATCH_COMBO}/geno_G2.vcf | \
          bcftools norm --rm-dup all \
          -O v -o ${SCRATCH_COMBO}/geno_G2.vcf.tmp

      mv ${SCRATCH_COMBO}/geno_G2.vcf.tmp ${SCRATCH_COMBO}/geno_G2.vcf
      rm -f ${CONTIG_HEADER}

      # Step 2b: BEAGLE phasing
      echo "  Step 2b: BEAGLE phasing"
      java -Xmx50g -jar ${EBROOTBEAGLE}/beagle.jar \
          gt=${SCRATCH_COMBO}/geno_G2.vcf \
          nthreads=${NCORES} \
          out=${SCRATCH_COMBO}/geno_G2_phased

      gunzip ${SCRATCH_COMBO}/geno_G2_phased.vcf.gz

      # Step 2c: convert-from-vcf
      echo "  Step 2c: convert-from-vcf"
      ./convert-from-vcf \
          -i ${SCRATCH_COMBO}/geno_G2_phased.vcf \
          --map ${MAP} \
          --genofolder ${SCRATCH_COMBO}/genoG2 \
          --hapfolder ${SCRATCH_COMBO}/haploG2

      # Step 2d: haplotype-hybrid MH blocks (independent G2)
      echo "  Step 2d: haplotype-hybrid MH blocks (independent)"
      mkdir -p ${SCRATCH_COMBO}/mh_info_ld_haploblock
      if ! run_with_retry ./haplotype-hybrid \
          --ncores ${NCORES} \
          -i ${SCRATCH_COMBO}/haploG2/chr* \
          -m ${MAP} \
          -o ${SCRATCH_COMBO}/mh_info_ld_haploblock \
          --generate-genotypes ${SCRATCH_COMBO}/mh_genotypes \
          -v \
          ld-haploblock micro --window-bp 125; then
        echo "FAILED Step 2d: ${COMBO}"
        rm -rf ${SCRATCH_COMBO}
        continue
      fi

      # Step 2e: haplotype-hybrid ALL blocks (independent G2)
      echo "  Step 2e: haplotype-hybrid ALL blocks (independent)"
      if ! run_with_retry ./haplotype-hybrid \
          --ncores ${NCORES} \
          -i ${SCRATCH_COMBO}/haploG2/chr* \
          -m ${MAP} \
          -o ${SCRATCH_COMBO}/all_info_haploblock \
          --generate-genotypes ${SCRATCH_COMBO}/all_genotypes \
          -v \
          fixed-kb --window-bp 100000; then
        echo "FAILED Step 2e: ${COMBO}"
        rm -rf ${SCRATCH_COMBO}
        continue
      fi
      
      # Step 2f: compress dan rsync ke QRIS
      echo "  Step 2f: compress and rsync to QRIS"
      TAR_FILE=${SCRATCH_TMP}/iter${ITER}_${COMBO}.tar.gz
      tar -czf ${TAR_FILE} -C ${SCRATCH_COMBO} .

      if [ $? -eq 0 ]; then
        tar -tzf ${TAR_FILE} > /dev/null
        if [ $? -eq 0 ]; then
          mkdir -p ${QRIS_TAR_DIR}
          rsync -a ${TAR_FILE} ${QRIS_TAR}
          rm -f ${TAR_FILE}
          echo "  Archive OK: ${COMBO}.tar.gz → QRIS"
          # Hapus file individual dari QRIS — sudah ada dalam tar
          rm -rf ${SCRATCH_COMBO}
          echo "  Done: ${COMBO}: $(date)"
        else
          echo "  ARCHIVE CORRUPT: ${COMBO} — scratch preserved"
        fi
      else
        echo "  COMPRESS FAILED: ${COMBO} — scratch preserved"
        rm -f ${TAR_FILE}
      fi

    done
  done
done

# Cleanup scratch iter
rm -rf ${SCRATCH_TMP}/iter${ITER}

echo "=== MH pipeline G2 iter=${ITER} complete: $(date) ==="