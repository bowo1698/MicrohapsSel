# Pipeline for genomic prediction using microhaplotype markers

## General information

We adapted and optimised the GVCHAP pipeline from Prakapenka et al to handle large genomic data and utilise microhaplotype markers in genomic prediction. The pipeline consists of four key stages: (1) phasing, (2) converting phased data to haplotypes, (3) defining microhaplotype blocks, and finally, (4) encoding them into genotypes. The pipeline was optimised with Rust, so it has more efficient memory and enhanced scalability due to parallelisation. For phasing, we only used Beagle as the pipeline has not been developed for other phasing methods yet (e.g., FindHap, FImpute, etc.). In addition, we encourage the use of parallelisation techniques when phasing with Beagle. 

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




```bash
chmod +x convert-from-vcf && xattr -d com.apple.quarantine convert-from-vcf
chmod +x convert-to-vcf && xattr -d com.apple.quarantine convert-to-vcf
chmod +x haplotype-hybrid && xattr -d com.apple.quarantine haplotype-hybrid
```
