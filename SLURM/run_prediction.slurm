#!/bin/bash
#SBATCH --job-name=genomic_array
#SBATCH --partition=general
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --exclusive
#SBATCH --time=336:00:00             
#SBATCH --account=a_zenger
#SBATCH --array=0-7                  
#SBATCH -o logs/array_%A_%a.out     
#SBATCH -e logs/array_%A_%a.err

# Threading konsisten - otomatis sesuai SLURM allocation
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export FLEXIBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
module load r/4.4.0-gfbf-2023a

# Ambil informasi hardware spesifik node secara real-time
NODE_NAME=$(hostname)
CPU_MODEL=$(lscpu | grep "Model name" | sed 's/Model name:[[:space:]]*//')

# Daftar config
CONFIGS=(
    "config_mh_normal.R"
    "config_mh_mixture.R"
    "config_snp_normal.R"
    "config_snp_mixture.R"
    "config_snp_normal_mh.R"
    "config_snp_mixture_mh.R"
    "config_snp_normal_rand.R"
    "config_snp_mixture_rand.R"
)

echo "Running synthetic benchmark..."
BENCH_START=$(date +%s.%N)

# Menjalankan benchmark dengan proteksi internal di dalam R
Rscript -e "
if (requireNamespace('RhpcBLASctl', quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads($OMP_NUM_THREADS)
};
x <- matrix(rnorm(10000*10000), 10000, 10000); 
y <- x %*% t(x); 
cat(sum(y))
" > /dev/null

BENCH_END=$(date +%s.%N)

# Cara paling aman menghitung selisih waktu (mendukung awk jika bc tidak ada)
BENCH_TIME=$(awk -v start=$BENCH_START -v end=$BENCH_END 'BEGIN {print end - start}')

CFG=${CONFIGS[$SLURM_ARRAY_TASK_ID]}
SUMMARY_FILE="metrics_summary_${SLURM_ARRAY_TASK_ID}.txt"
DETAILED_FILE="metrics_detailed_${SLURM_ARRAY_TASK_ID}.csv"

# Inisialisasi header file
echo "config,iteration,generation,duration_seconds,node,cpu_type,bench_time" > $SUMMARY_FILE
echo "config,iter,gen,user_sec,sys_sec,elapsed_sec,max_ram_kb,node,cpu_type" > $DETAILED_FILE

run_analysis() {
    local config=$1
    local iter=$2
    local gen=$3
    local script=$4
    
    echo "=== Processing $config: iter $iter, gen $gen pada $NODE_NAME ==="
    export GENOMIC_CONFIG=$config
    export ITER=$iter
    export GEN=$gen
    
    start_time=$(date +%s)
    
    # Mencatat metrik performa + info hardware ke file CSV
    /usr/bin/time -f "$config,$iter,$gen,%U,%S,%e,%M,$NODE_NAME,$CPU_MODEL" -a -o $DETAILED_FILE Rscript $script

    end_time=$(date +%s)
    duration=$((end_time - start_time))
    
    # Mencatat ke summary file dengan bench_time
    echo "$config,$iter,$gen,$duration,$NODE_NAME,$CPU_MODEL,$BENCH_TIME" >> $SUMMARY_FILE
    
    sync
    echo "Done. Duration: $duration seconds."
    echo "------------------------------------------------"
}

# Pilih script berdasarkan nama config
if [[ $CFG == *"mh"* ]] && [[ $CFG != *"snp"* ]]; then
    RSCRIPT="main-mh.R"
else
    RSCRIPT="main-snp.R"
fi

echo "========================================"
echo "Array ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $NODE_NAME"
echo "CPU: $CPU_MODEL"
echo "Benchmark time: $BENCH_TIME seconds"
echo "Config: $CFG"
echo "========================================"

# Loop iterasi (1-5) dan generasi (1-3)
for iter in {1..5}; do  
    for gen in {1..3}; do
        run_analysis $CFG $iter $gen $RSCRIPT
    done
done

echo "Task $SLURM_ARRAY_TASK_ID (Config $CFG) completed!"
