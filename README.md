# Pipeline for genomic prediction using microhaplotype markers

## General information

We adapted and optimised the GVCHAP pipeline from Prakapenka et al to handle large genomic data and utilise microhaplotype markers in genomic prediction. The pipeline consists of four key stages: (1) phasing, (2) converting phased data to haplotypes, (3) defining microhaplotype blocks, and finally, (4) encoding them into microhaplotype genotypes. The pipeline was optimised with Rust, so it has more efficient memory and enhanced scalability due to parallelisation. For phasing, we only used Beagle as the pipeline has not been developed for other phasing methods yet (e.g., FindHap, FImpute, etc.). In addition, we encourage the use of parallelisation techniques when phasing with Beagle. 

While haplotypes refer to long-range combinations of alleles along a chromosome, microhaplotypes can be defined as several SNPs covering a short DNA segment, typically 125 - 150 bp. So in this pipeline, conducted by `haplotype-hybrid`, we provide fixed window (SNP counts) and LD-based methods to define fully haplotype blocks and microhaplotype segments, respectively. When the `--method snp-count-simple` argument is used, it simply defines haplotype blocks depending on how many SNPs are set in each block. In contrast, when the `--method ld-haploblock` along with `--haplotype-type micro` argument is used, it will discover microhaplotypes based on LD and window constraints. However, both methods will result in similar data outputs.

## Data preparation

Raw genotype data is expected in a .VCF format, as it contains genetic map positions, reference/alternative alleles, and unphased genotype calls (e.g., 0/1), which serve as the essential input for the subsequent phasing step. However, if the raw geneotype data is in PLINK formats (provided as .bed, .bim, .fam), recoding to .VCF format can be done by PLINK as well, for example:

```bash
plink --bfile GENO_DATA --chr-set 95 no-xy --recode vcf --out GENO_DATA
```

This will generate unphased genotype data in a .VCF file, for example:

```
##fileformat=VCFv4.2
##fileDate=20250503
##source=PLINKv1.9
##contig=<ID=1,length=80480704>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Ind_1	Ind_2
1	16977	snp1	C	A	.	.	PR	GT	0/0	0/1
1	33723	snp2	G	T	.	.	PR	GT	0/0	0/0
```

Alternatively, if the genotype data is in form of .CSV:

```
ID,snp_1,snp_2,snp_3,...snp_n
ind_1,0,1,2,...
ind_2,1,2,0,...
ind_3,2,0,1,...
```

Along with a genetic map (in form .txt):

```
SNPID	Chr	Position
snp_1	1	0
snp_2	1	16746
snp_3	1	217334
```

We provide `convert-to-vcf` to convert such genotype data to a VCF file. This tool can be run as follows:

```
./convert-to-vcf \
      -i genotypes.csv \
      -m map.txt \
      -o genotypes.vcf \
      --split-by-chr \
      --parallel \
      --ncores 64 \
      --merge-after-split
```

Note that the conversion is done by splitting the genotype files based on the number of chromosomes processed in parallel, then merging them into a single VCF file.

## Phasing

Phasing is the first critical step in our pipeline, where ambiguous diploid genotypes are separated into their constituent haploid sequences to allow for the accurate identification of microhaplotype blocks. We encourage to do phasing using Beagle in parallel process, and to do that, we use `bcftools` to split the VCF file per chromosome and `tabix` for indexing a compressed VCF, as an example:

```bash
for chr in {1..29}; do
    bcftools view -r ${chr} /data/genotypes.vcf.gz -O z -o /data/tmp/chr${chr}.vcf.gz
    tabix -p vcf /data/tmp/chr${chr}.vcf.gz
done

for chr in {1..29}; do
    (
        java -Xmx28g -jar ${EBROOTBEAGLE}/beagle.jar \
            gt=/data/tmp/chr${chr}.vcf.gz \
            nthreads=4 \
            out=/data/tmp/chr${chr}_phased
    ) &
    
    if (( chr % 8 == 0 )); then
        wait
    fi
done
wait

for chr in {1..29}; do
    tabix -p vcf /data/tmp/chr${chr}_phased.vcf.gz
    bcftools reheader -h <(bcftools view -h /data/genotypes.vcf.gz) /data/tmp/chr${chr}_phased.vcf.gz | \
        bcftools view -O z -o /data/tmp/chr${chr}_phased_corrected.vcf.gz
    tabix -p vcf /data/tmp/chr${chr}_phased_corrected.vcf.gz
done

bcftools concat -O z -o /data/genotypes_phased.vcf.gz \
    /data/tmp/chr{1..29}_phased_corrected.vcf.gz

tabix -p vcf /data/genotypes_phased.vcf.gz

/usr/bin/rm -rf /data/tmp
gunzip /data/genotypes_phased.vcf.gz
```

This process results in a phased VCF file, for example:

```
##fileformat=VCFv4.2
##filedate=20260108
##source="beagle.22Jul22.46e.jar"
##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated ALT Allele Frequencies">
##INFO=<ID=DR2,Number=A,Type=Float,Description="Dosage R-Squared: estimated squared correlation between estimated REF dose [P(RA) + 2*P(RR)] and true REF dose">
##INFO=<ID=IMP,Number=0,Type=Flag,Description="Imputed marker">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DS,Number=A,Type=Float,Description="estimated ALT dose [P(RA) + 2*P(AA)]">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Ind_1	Ind_2
1	16977	snp1	C	A	.	.	PR	GT	0|0	1|1	
1	33723	snp2	G	T	.	.	PR	GT	1|0	0|0
```

## Convertion from phased genotypes to haplotypes

After phasing, we want to separate the phased haplotypes into individual haploid sequences for haploblock detection. Here, we provide `convert-from-vcf` for conversion, which accomplishes two key transformations:

**1. Haplotype separation**
Each phased genotype (e.g., `0|1`) is split into two separate haplotype columns:
- Individual_A with genotype `0|1` becomes:
  - Haplotype 1: `1` (recoded from `0`)
  - Haplotype 2: `2` (recoded from `1`)

This separation is essential because LD-based haploblock detection requires calculating pairwise D' between SNPs, which depends on counting the four possible haplotype combinations (00, 01, 10, 11).

**2. Allele recoding**
VCF format uses `0` (reference) and `1` (alternate) encoding. The conversion recodes these to:
- `0` → `1` (major allele)
- `1` → `2` (minor allele)

For genotypes (diploid coding):
- `0|0` → `0` (homozygous major)
- `0|1` or `1|0` → `1` (heterozygous)
- `1|1` → `2` (homozygous minor)

To handle large datasets, individuals are processed in chunks (batches) rather than all at once. The chunk size is determined by the `--interval` parameter:
- `--interval 1`: Process all individuals in one pass (high memory)
- `--interval 10`: Process in 10 batches (lower memory, slightly slower)

Each chunk reads through the VCF file sequentially, extracts the relevant individuals, and writes output incrementally. Data is automatically split by chromosome, creating separate output files (haplotypes and genotypes), for example:

**Haplotype Files (`hap/chr*.hap`)**
Format: Each individual contributes **two columns** (one per haplotype)
```
ID              SNP1_h1  SNP1_h2  SNP2_h1  SNP2_h2  ...
Individual_A    1        2        1        1        ...
Individual_B    2        2        1        2        ...
```

This provides input for LD-based haploblock detection

**Genotype Files (`geno/chr*.geno`)**
Format: Traditional genotype matrix with header
```
ID              SNP1  SNP2  SNP3  ...
Individual_A    1     0     2     ...
Individual_B    2     1     1     ...
```

This provides input for GBLUP or other genotype-based prediction models that utilises biallelic SNPs.

### Usage Example
```bash
./convert-from-vcf \
      -i /data/genotypes_phased.vcf.gz \
      --map map.txt \
      --genofolder geno \
      --hapfolder hap \
      --nosort
```

### Key Parameters

- `-i, --input`: VCF file(s) to convert (required)
- `--interval`: Number of chunks for processing (default: 1)
- `--hapfolder`: Output directory for haplotype files (default: `hap`)
- `--genofolder`: Output directory for genotype files (default: `geno`)
- `--map`: Output path for SNP map file (default: `map_new.txt`)
- `--missing`: Code for missing genotypes (default: `-9999`)
- `-V, --verbose`: Display detailed progress information

The separated haplotype files (`hap/chr*`) then serve as direct input for the next step: defining microhaplotype segments using LD-based haploblock detection.

## Discovering microhaplotype segments and genotyping

Because our purpose is to discover microhaplotype segments, we only use outputs from `hap/chr*` and perform arguments `--method ld-haploblock` and `--haplotype-type micro`. The `haplotype-hybrid` tool works with finding all potential haplotype candidates that meet an LD criteria (e.g., D' > 0.45) as the argument of `--d-prime-threshold 0.45` is applied. As a result, there may be less haplotype blocks are defined. After that, align with the microhaplotype definition, they are selected from the haplotype block candidates, where in each segment should contain consecutive SNPs within 125 or 150 bp which then they are evaluated using Criterion-B following Jónás et al. (2017), to ensure balance between allele frequency and microhaplotype diversity, following the calculation:

$$CriterionB_{m h_i}=\sum_{k=1}^{N_i}\left(f_i-\frac{1}{H S}\right)^2-w N_i$$

In the first term, $m h_i$ denotes microhaplotype $i$, $f_i$ is the frequency of microhaplotype allele, and $HS=2^n$ is the theoretical maximum number of alleles for n SNPs. In the second term, $w$ is calculated as $(MD.N_i)/(HS.(N_i-1) )$, where $MD$ is the scaling parameter of maximum deviation to control the magnitude of microhaplotype diversity. We set the $MD$ as $0.1$ following (Jónás et al., 2017). The $N_i$ represents the number of predictable alleles for microhaplotype $i$ that have a frequency above the allele frequency threshold ($AFT≥0.08$). 

```bash
chmod +x convert-from-vcf && xattr -d com.apple.quarantine convert-from-vcf
chmod +x convert-to-vcf && xattr -d com.apple.quarantine convert-to-vcf
chmod +x haplotype-hybrid && xattr -d com.apple.quarantine haplotype-hybrid
```

