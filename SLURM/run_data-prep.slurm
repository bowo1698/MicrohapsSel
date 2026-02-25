#!/bin/bash
#SBATCH --job-name=run_trout_all2
#SBATCH --partition=general
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=100:00:00
#SBATCH --account=a_zenger
#SBATCH -o logs/trout_all2_%j.out
#SBATCH -e logs/trout_all2_%j.err

OUTDIR="output/trout2"
SPECIES="trout"

# Step 1: Simulation
module load r/4.4.0-gfbf-2023a
Rscript simulation_trout_gens.R

## Output validation
echo "========================================="
echo "VALIDATING SIMULATION OUTPUT"
echo "========================================="

MAP_SNPS=$(tail -n +2 ${OUTDIR}/${SPECIES}_snp_map.txt | wc -l)
echo "Map file SNPs: ${MAP_SNPS}"

for iter in {1..5}; do
  for gen in {1..3}; do
    GENO_COLS=$(head -1 ${OUTDIR}/iteration_${iter}/${SPECIES}_gen${gen}_genotypes.csv | tr ',' '\n' | wc -l)
    GENO_SNPS=$((GENO_COLS - 3))  # Subtract ID, generation, population_type
    
    if [ ${GENO_SNPS} -ne ${MAP_SNPS} ]; then
      echo "ERROR: iteration ${iter} gen ${gen}: Genotype has ${GENO_SNPS} SNPs, map has ${MAP_SNPS}"
      exit 1
    fi
  done
done

echo "✓ All genotype files match map file (${MAP_SNPS} SNPs)"
echo "✓ Validation PASSED - safe to proceed"
echo "========================================="

module purge

# Step 2: Phasing (load Beagle)
module load beagle/5.4.22jul22.46e-java-11
module load anaconda3/2023.09-0
source activate microhaplotype

for iter in {1..5}; do
  for gen in {1..3}; do
    ./convert-to-vcf \
      -i ${OUTDIR}/iteration_${iter}/${SPECIES}_gen${gen}_genotypes.csv \
      -m ${OUTDIR}/${SPECIES}_snp_map.txt \
      -o ${OUTDIR}/iteration_${iter}/${SPECIES}_gen${gen}_genotypes.vcf \
      --split-by-chr \
      --parallel \
      --ncores 16 \
      --merge-after-split
    
    java -Xmx450g -jar ${EBROOTBEAGLE}/beagle.jar \
      gt=${OUTDIR}/iteration_${iter}/${SPECIES}_gen${gen}_genotypes.vcf \
      nthreads=16 \
      out=${OUTDIR}/iteration_${iter}/${SPECIES}_sim_phased_gen${gen}
    
    gunzip ${OUTDIR}/iteration_${iter}/${SPECIES}_sim_phased_gen${gen}.vcf.gz
    
    ./convert-from-vcf \
      -i ${OUTDIR}/iteration_${iter}/${SPECIES}_sim_phased_gen${gen}.vcf \
      --map ${OUTDIR}/${SPECIES}_snp_map.txt \
      --genofolder ${OUTDIR}/iteration_${iter}/geno_gen${gen} \
      --hapfolder ${OUTDIR}/iteration_${iter}/haplo_gen${gen} \
      --nosort
  done
done

# Step 3: Genotyping
for iter in {1..5}; do
  for gen in {1..3}; do
    ./haplotype-hybrid \
      --method ld-haploblock \
      --haplotype-type micro \
      --d-prime-threshold 0.45 \
      --window-bp 125 \
      --min-snps 2 \
      --aft 0.08 \
      --md 0.10 \
      -i ${OUTDIR}/iteration_${iter}/haplo_gen${gen}/chr{1..29} \
      -m ${OUTDIR}/${SPECIES}_snp_map.txt \
      -o ${OUTDIR}/iteration_${iter}/mh_info_ld_haploblock_gen${gen} \
      --generate-genotypes mh_genotypes_gen${gen} \
      -v
    
    python extract_selected_genotypes.py \
      -g ${OUTDIR}/iteration_${iter}/${SPECIES}_gen${gen}_genotypes.csv \
      -s ${OUTDIR}/iteration_${iter}/mh_info_ld_haploblock_gen${gen}/stats/snp_selection_detailed.csv \
      -o ${OUTDIR}/iteration_${iter}/genotypes_microhap_gen${gen}.csv \
      -v
  done
done