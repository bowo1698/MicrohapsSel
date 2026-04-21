#!/bin/bash
#SBATCH --job-name=mh_G1
#SBATCH --partition=general
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=22
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --account=a_zenger
#SBATCH --array=1-20
#SBATCH -o logs/mh_G1_%A_%a.out
#SBATCH -e logs/mh_G1_%A_%a.err

set -euo pipefail

module load anaconda3/2023.09-0
source activate microhaplotype
module load beagle/5.4.22jul22.46e-java-11

cd /scratch/user/aguswibowo/Research/simulation_forward

BASE=/scratch/user/aguswibowo/Research/simulation_forward
OUTDIR=${BASE}/output_G1
NCORES=20
MAP=${BASE}/raw_data/salmon_60k_map_final.txt
REF_VCF=${BASE}/raw_data/salmon_60k_pairs125bp.vcf.gz

ITER=$(printf "%02d" ${SLURM_ARRAY_TASK_ID})
ITER_DIR=${OUTDIR}/iter${ITER}

echo "Starting MH pipeline iter=${ITER}"

echo "Step 2a: convert-to-vcf iter=${ITER}"
./convert-to-vcf \
    -i ${ITER_DIR}/geno_G1.csv \
    -m ${MAP} \
    -o ${ITER_DIR}/geno_G1.vcf \
    --split-by-chr --parallel --ncores ${NCORES} --merge-after-split

echo "Step 2a.5: Remove duplicate positions iter=${ITER}"
CONTIG_HEADER=${ITER_DIR}/contig_header.txt
bcftools view -h ${REF_VCF} | grep "^##contig" > ${CONTIG_HEADER}

bcftools annotate \
    --header-lines ${CONTIG_HEADER} \
    ${ITER_DIR}/geno_G1.vcf | \
    bcftools norm --rm-dup all \
    -O v -o ${ITER_DIR}/geno_G1.vcf.tmp

mv ${ITER_DIR}/geno_G1.vcf.tmp ${ITER_DIR}/geno_G1.vcf
rm -f ${CONTIG_HEADER}

echo "Step 2b: BEAGLE phasing iter=${ITER}"
java -Xmx50g -jar ${EBROOTBEAGLE}/beagle.jar \
    gt=${ITER_DIR}/geno_G1.vcf \
    nthreads=${NCORES} \
    out=${ITER_DIR}/geno_G1_phased

gunzip ${ITER_DIR}/geno_G1_phased.vcf.gz

echo "Step 2c: convert-from-vcf iter=${ITER}"
./convert-from-vcf \
    -i ${ITER_DIR}/geno_G1_phased.vcf \
    --map ${MAP} \
    --genofolder ${ITER_DIR}/genoG1 \
    --hapfolder ${ITER_DIR}/haploG1

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

echo "Step 2d: haplotype-hybrid MH blocks iter=${ITER}"
mkdir -p ${ITER_DIR}/mh_info_ld_haploblock
run_with_retry ./haplotype-hybrid \
    --ncores ${NCORES} \
    -i ${ITER_DIR}/haploG1/chr* \
    -m ${MAP} \
    -o ${ITER_DIR}/mh_info_ld_haploblock \
    --generate-genotypes ${ITER_DIR}/mh_genotypes \
    -v \
    ld-haploblock micro --window-bp 125

echo "Step 2e: haplotype-hybrid ALL blocks iter=${ITER}"
run_with_retry ./haplotype-hybrid \
    --ncores ${NCORES} \
    -i ${ITER_DIR}/haploG1/chr* \
    -m ${MAP} \
    -o ${ITER_DIR}/all_info_haploblock \
    --generate-genotypes ${ITER_DIR}/all_genotypes \
    -v \
    fixed-kb --window-bp 100000

echo "Done MH pipeline iter=${ITER}"