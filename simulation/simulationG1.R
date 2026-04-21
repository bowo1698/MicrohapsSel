suppressPackageStartupMessages({
  library(AlphaSimR)
  library(tidyverse)
  library(vcfR)
})

# ── Config ────────────────────────────────────────────────────────────────────
BASE        <- "/scratch/user/aguswibowo/Research/simulation_forward_2/raw_data"
OUT         <- "/scratch/user/aguswibowo/Research/simulation_forward_2/output_G1"
VCF_FILE    <- file.path(BASE, "salmon_60k_pairs125bp.vcf")
MH_STATS    <- file.path(BASE, "mh_info_ld_haploblock_G0/stats")
HAP_STATS   <- file.path(BASE, "all_info_haploblock_G0/stats")
MH_INFO     <- file.path(BASE, "mh_info_ld_haploblock_G0")
HAP_INFO    <- file.path(BASE, "all_info_haploblock_G0")
MH_GENO     <- file.path(BASE, "mh_genotypes_G0")
HAP_GENO    <- file.path(BASE, "all_genotypes_G0")

args        <- commandArgs(trailingOnly=TRUE)
ITER        <- as.integer(args[1])
N_ITER      <- 1
N_OFF_FAM   <- 25
N_QTL       <- 100
H2          <- 0.3
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
RDS_DIR  <- file.path(OUT, "rds")
ITER_DIR <- file.path(OUT, sprintf("iter%02d", ITER))
dir.create(RDS_DIR,  recursive=TRUE, showWarnings=FALSE)
dir.create(ITER_DIR, recursive=TRUE, showWarnings=FALSE)

# ── 1. Read VCF → MapPop ──────────────────────────────────────────────────────
vcf       <- read.vcfR(VCF_FILE, verbose=FALSE)

# Extract haplotype matrix dari VCF phased
gt        <- extract.gt(vcf, element="GT")
haplo_mat <- do.call(cbind, lapply(1:ncol(gt), function(i) {
  m <- do.call(rbind, strsplit(gt[, i], "\\|"))
  matrix(as.integer(m), nrow=nrow(m), ncol=2)
}))

# Genetic map
gen_map   <- data.frame(
  chr = as.integer(vcf@fix[,"CHROM"]),
  pos = as.integer(vcf@fix[,"POS"]) / 1e6  # convert bp → Morgan (approx)
)

chr_list   <- split(gen_map$chr, gen_map$chr)
n_ind_vcf  <- ncol(gt)  # n ind

chrs_ordered <- sort(unique(gen_map$chr))

genMap_list <- lapply(chrs_ordered, function(c) {
  p <- gen_map$pos[gen_map$chr == c]
  p - min(p)
})


haplo_list <- lapply(chrs_ordered, function(c) {
  snp_idx <- which(gen_map$chr == c)
  t(haplo_mat[snp_idx, , drop=FALSE])  # (n_ind*2) × n_snp
})

founderPop <- newMapPop(
  genMap     = genMap_list,
  haplotypes = haplo_list,
  inbred     = FALSE
)

vcf_fix  <- as.data.frame(getFIX(vcf))
snp_meta <- tibble(snp_id = vcf_fix$ID,
                   chr    = as.integer(vcf_fix$CHROM),
                   pos    = as.integer(vcf_fix$POS)) %>%
            arrange(chr, pos)

n_snp    <- nrow(snp_meta)

# ── 2. Load block info ────────────────────────────────────────────────────────
mh_coords  <- read_csv(file.path(MH_STATS,  "microhaplotype_coordinates.csv"), show_col_types=FALSE)
hap_coords <- read_csv(file.path(HAP_STATS, "microhaplotype_coordinates.csv"), show_col_types=FALSE)
mh_snps    <- read_csv(file.path(MH_STATS,  "snp_selection_detailed.csv"),     show_col_types=FALSE)
hap_snps   <- read_csv(file.path(HAP_STATS, "snp_selection_detailed.csv"),     show_col_types=FALSE)

# ── 3. Define QTL effects once from G0 ───────────────────────────────────────
set.seed(123)

## S1/S2: SNP QTL — pick N_QTL SNPs randomly
qtl_snp_ids <- sample(snp_meta$snp_id, N_QTL)

define_snp_effects <- function(qtl_snp_ids, type="equal") {
  n <- length(qtl_snp_ids)
  eff <- if (type=="equal") rep(1/sqrt(n), n) else rgamma(n, shape=0.4, scale=1)
  eff <- eff / sqrt(sum(eff^2))  # normalize
  tibble(snp_id=qtl_snp_ids, effect=eff)
}

qtl_s1 <- define_snp_effects(qtl_snp_ids, "equal")
qtl_s2 <- define_snp_effects(qtl_snp_ids, "unequal")

## S3/S4: MH block QTL
qtl_mh_blocks <- sample(mh_coords$block_id, N_QTL)

## S5/S6: Haplotype block QTL
qtl_hap_blocks <- sample(hap_coords$block_id, N_QTL)

# Read G0 block genotypes → build allele effect lookup per block
read_block_geno_G0 <- function(geno_dir, info_dir, block_ids, block_snp_map) {
  by_chr <- block_snp_map %>% filter(block_id %in% block_ids) %>%
    group_by(block_id) %>% summarise(chr=first(chr), .groups="drop")

  # build global→local index per chr dari hap_block_* file
  local_idx_map <- map_dfr(unique(by_chr$chr), function(c_val) {
    blk_file <- file.path(info_dir, paste0("hap_block_", c_val))
    blk_ids  <- gsub("\t.*", "", readLines(blk_file))
    tibble(block_id=blk_ids, chr=c_val, local_idx=seq_along(blk_ids))
  })

  allele_effects <- list()
  for (blk in block_ids) {
    c_val     <- by_chr$chr[by_chr$block_id == blk]
    local_num <- local_idx_map$local_idx[local_idx_map$block_id == blk]
    geno_file <- file.path(geno_dir, paste0("hap_geno_", c_val))
    geno      <- read.delim(geno_file, header=FALSE, check.names=FALSE, stringsAsFactors=FALSE)
    raw_hdr   <- scan(geno_file, what="", nlines=1, quiet=TRUE)
    target  <- paste0("hap_", c_val, "_", local_num)
    matches <- which(raw_hdr == target)
    blk_pos <- matches[1]  # paternal
    mat_pos <- matches[2]  # maternal
    a1_vec  <- geno[[blk_pos]]
    a2_vec  <- geno[[mat_pos]]
    all_alleles <- unique(c(a1_vec, a2_vec))
    allele_effects[[blk]] <- setNames(rnorm(length(all_alleles)), as.character(all_alleles))
  }
  allele_effects
}

# Normalize: equal → standardize per block; unequal → divide by sqrt(n_qtl)
normalize_block_effects <- function(eff_list, type="equal") {
  n <- length(eff_list)
  lapply(eff_list, function(e) {
    if (type=="equal") e / (sd(e) * sqrt(n)) else {
      # gamma: scale by sum of squares so that magnitude comparable
      e / sqrt(sum(e^2))
    }
  })
}

qtl_s3_effects <- normalize_block_effects(
  read_block_geno_G0(MH_GENO,  MH_INFO,  qtl_mh_blocks,  mh_snps), "equal")
qtl_s4_effects <- normalize_block_effects(
  read_block_geno_G0(MH_GENO,  MH_INFO,  qtl_mh_blocks,  mh_snps), "unequal")
qtl_s5_effects <- normalize_block_effects(
  read_block_geno_G0(HAP_GENO, HAP_INFO, qtl_hap_blocks, hap_snps), "equal")
qtl_s6_effects <- normalize_block_effects(
  read_block_geno_G0(HAP_GENO, HAP_INFO, qtl_hap_blocks, hap_snps), "unequal")

# ── 4. TBV calculation functions ──────────────────────────────────────────────
calc_tbv_snp <- function(dosage_mat, qtl_effects) {
  # dosage_mat: individuals × SNPs, qtl_effects: tibble snp_id/effect
  cols <- intersect(qtl_effects$snp_id, colnames(dosage_mat))
  eff  <- qtl_effects$effect[match(cols, qtl_effects$snp_id)]
  tbv  <- as.matrix(dosage_mat[, cols]) %*% eff
  as.vector(tbv - mean(tbv))
}

calc_tbv_block <- function(geno_dir, ind_ids, block_ids, eff_list, block_snp_map, generation_haplo) {
  # generation_haplo: list per chr of matrix (n_ind*2 × n_snp) — haplotypes of G1
  # For simulation: reconstruct allele integers from haplotype data
  by_chr <- block_snp_map %>% filter(block_id %in% block_ids) %>%
    group_by(block_id) %>% summarise(chr=first(chr), snp_indices=list(snp_index), .groups="drop")
  
  n_ind <- length(ind_ids)
  tbv   <- numeric(n_ind)
  
  for (bi in seq_len(nrow(by_chr))) {
    blk      <- by_chr$block_id[bi]
    c_val    <- by_chr$chr[bi]
    snp_idx  <- unlist(by_chr$snp_indices[bi]) + 1  # 1-based
    haplo    <- generation_haplo[[c_val]]            # (n_ind*2) × n_snp
    eff      <- eff_list[[blk]]
    
    for (i in seq_len(n_ind)) {
      h1    <- haplo[(i-1)*2+1, snp_idx]
      h2    <- haplo[(i-1)*2+2, snp_idx]
      a1    <- strtoi(paste(h1, collapse=""), base=2L) + 1
      a2    <- strtoi(paste(h2, collapse=""), base=2L) + 1
      e1    <- if (!is.na(eff[as.character(a1)])) eff[as.character(a1)] else 0
      e2    <- if (!is.na(eff[as.character(a2)])) eff[as.character(a2)] else 0
      tbv[i] <- tbv[i] + e1 + e2
    }
  }
  tbv - mean(tbv)
}

# ── 5. Simulation loop ────────────────────────────────────────────────────────
SP <- SimParam$new(founderPop)
n_snp_per_chr <- snp_meta %>% count(chr) %>% pull(n)
SP$addSnpChip(nSnpPerChr=n_snp_per_chr)  # register all SNPs as chip

all_geno  <- list()
all_pheno <- list(S1=list(), S2=list(), S3=list(), S4=list(), S5=list(), S6=list())

all_pedigree <- list()

for (iter in ITER) {
  set.seed(iter * 100)
  
  # G0 population
  G0 <- newPop(founderPop, simParam=SP)
  
  # Random crosses → G1 families
  n_founders <- G0@nInd
  sires <- 1:floor(n_founders/2)
  dams  <- (floor(n_founders/2)+1):n_founders
  crosses <- data.frame(
    sire = sample(sires, min(length(sires), length(dams)), replace=FALSE),
    dam  = sample(dams,  min(length(sires), length(dams)), replace=FALSE)
  )
  
  G1_list <- lapply(seq_len(nrow(crosses)), function(f) {
    makeCross(G0, crossPlan=matrix(c(crosses$sire[f], crosses$dam[f]), nrow=1),
              nProgeny=N_OFF_FAM, simParam=SP)
  })
  G1 <- mergePops(G1_list)
  
  n_G1     <- G1@nInd
  fam_ids  <- rep(seq_len(nrow(crosses)), each=N_OFF_FAM)
  ind_ids  <- sprintf("iter%d_gen1_fam%03d_ind%02d",
                      iter, fam_ids, sequence(rep(N_OFF_FAM, nrow(crosses))))
  
  # Pedigree G1
  pedigree_g1 <- tibble(
    individual_id   = ind_ids,
    generation      = 1,
    population_type = "reference",
    mother_id       = rep(sprintf("founder_%d", crosses$dam),  each=N_OFF_FAM),
    father_id       = rep(sprintf("founder_%d", crosses$sire), each=N_OFF_FAM),
    iteration       = iter,
    family_id       = fam_ids
  )
  all_pedigree[[iter]] <- pedigree_g1

  # Save G1 pop object
  saveRDS(G1, file.path(ITER_DIR, "G1_pop.rds"))

  # Dosage genotype matrix
  dosage <- pullSnpGeno(G1, simParam=SP)
  colnames(dosage) <- snp_meta$snp_id
  rownames(dosage) <- ind_ids
  
  geno_df <- as_tibble(dosage, rownames="ID") %>%
    mutate(generation=1, population_type="reference", .after="ID")
  all_geno[[iter]] <- geno_df
  
  # Haplotype matrix G1 per chr (for block TBV)
  n_chr <- length(unique(snp_meta$chr))
  h_all <- pullSegSiteHaplo(G1, simParam=SP)
  G1_haplo <- lapply(1:n_chr, function(c) {
    snp_in_chr <- which(snp_meta$chr == c)
    h_all[, snp_in_chr, drop=FALSE]
  })
  
  # rescale TBV
  rescale_tbv <- function(tbv) (tbv - mean(tbv)) / sd(tbv)

  # SNP TBV (S1, S2)
  tbv_s1 <- rescale_tbv(calc_tbv_snp(dosage, qtl_s1))
  tbv_s2 <- rescale_tbv(calc_tbv_snp(dosage, qtl_s2))
  
  # Block TBV (S3–S6)
  tbv_s3 <- rescale_tbv(calc_tbv_block(MH_GENO,  ind_ids, qtl_mh_blocks,  qtl_s3_effects, mh_snps,  G1_haplo))
  tbv_s4 <- rescale_tbv(calc_tbv_block(MH_GENO,  ind_ids, qtl_mh_blocks,  qtl_s4_effects, mh_snps,  G1_haplo))
  tbv_s5 <- rescale_tbv(calc_tbv_block(HAP_GENO, ind_ids, qtl_hap_blocks, qtl_s5_effects, hap_snps, G1_haplo))
  tbv_s6 <- rescale_tbv(calc_tbv_block(HAP_GENO, ind_ids, qtl_hap_blocks, qtl_s6_effects, hap_snps, G1_haplo))
  
  # Simulate phenotype: y = TBV + e, e ~ N(0, Var(TBV)*(1-h2)/h2)
  sim_pheno <- function(tbv, h2=H2) {
    var_g <- var(tbv)
    var_e <- var_g * (1 - h2) / h2
    tbv + rnorm(length(tbv), 0, sqrt(var_e))
  }
  
  make_pheno_df <- function(tbv, ind_ids, iter, fam_ids, h2=H2) {
    pheno    <- sim_pheno(tbv, h2)
    h2_real  <- var(tbv) / var(pheno)
    tibble(individual_id  = ind_ids,
           tbv            = tbv,
           phenotype      = pheno,
           generation     = 1,
           population_type= "reference",
           iteration      = iter,
           family_id      = fam_ids,
           h2_realized    = h2_real)
  }
  
  all_pheno$S1[[iter]] <- make_pheno_df(tbv_s1, ind_ids, iter, fam_ids)
  all_pheno$S2[[iter]] <- make_pheno_df(tbv_s2, ind_ids, iter, fam_ids)
  all_pheno$S3[[iter]] <- make_pheno_df(tbv_s3, ind_ids, iter, fam_ids)
  all_pheno$S4[[iter]] <- make_pheno_df(tbv_s4, ind_ids, iter, fam_ids)
  all_pheno$S5[[iter]] <- make_pheno_df(tbv_s5, ind_ids, iter, fam_ids)
  all_pheno$S6[[iter]] <- make_pheno_df(tbv_s6, ind_ids, iter, fam_ids)
  
  cat("Iter", iter, "done | n_G1:", n_G1, "\n")
}

# Write output per iter
write_csv(all_geno[[ITER]], file.path(ITER_DIR, "geno_G1.csv"))
walk(names(all_pheno), ~write_csv(all_pheno[[.x]][[ITER]],
       file.path(ITER_DIR, sprintf("pheno_%s_G1.csv", .x))))
write_csv(all_pedigree[[ITER]], file.path(ITER_DIR, "pedigree_G1.csv"))

# QTL + founderPop + pedigree G0 hanya ditulis oleh iter 1
if (ITER == 1) {
  pedigree_g0 <- tibble(
    individual_id   = sprintf("founder_%d", 1:n_ind_vcf),
    generation      = 0,
    population_type = "founder",
    mother_id       = NA, father_id = NA,
    iteration       = NA, family_id = NA
  )
  write_csv(pedigree_g0, file.path(OUT, "pedigree_G0.csv"))
  saveRDS(founderPop,     file.path(RDS_DIR, "founderPop.rds"))
  saveRDS(SP,             file.path(RDS_DIR, "SP.rds"))
  saveRDS(qtl_s1,         file.path(RDS_DIR, "qtl_s1.rds"))
  saveRDS(qtl_s2,         file.path(RDS_DIR, "qtl_s2.rds"))
  saveRDS(qtl_s3_effects, file.path(RDS_DIR, "qtl_s3_effects.rds"))
  saveRDS(qtl_s4_effects, file.path(RDS_DIR, "qtl_s4_effects.rds"))
  saveRDS(qtl_s5_effects, file.path(RDS_DIR, "qtl_s5_effects.rds"))
  saveRDS(qtl_s6_effects, file.path(RDS_DIR, "qtl_s6_effects.rds"))
  saveRDS(qtl_snp_ids,    file.path(RDS_DIR, "qtl_snp_ids.rds"))
  saveRDS(qtl_mh_blocks,  file.path(RDS_DIR, "qtl_mh_blocks.rds"))
  saveRDS(qtl_hap_blocks, file.path(RDS_DIR, "qtl_hap_blocks.rds"))
  saveRDS(snp_meta,       file.path(RDS_DIR, "snp_meta.rds"))
}

# Summary per iter
pheno_s1 <- all_pheno$S1[[ITER]]
summary_lines <- c(
  "=== Simulation Summary ===",
  sprintf("Date       : %s", Sys.time()),
  sprintf("ITER       : %d", ITER),
  sprintf("N_QTL      : %d", N_QTL),
  sprintf("H2_target  : %.2f", H2),
  sprintf("N_founders : %d", founderPop@nInd),
  sprintf("N_chr      : %d", n_chr),
  sprintf("N_SNP_total: %d", n_snp),
  sprintf("N_MH_blocks: %d", nrow(mh_coords)),
  sprintf("N_HAP_blocks: %d", nrow(hap_coords)),
  sprintf("n_G1       : %d", n_G1),
  "",
  "=== h2_realized per scenario ===",
  sprintf("S1: %.4f", unique(all_pheno$S1[[ITER]]$h2_realized)),
  sprintf("S2: %.4f", unique(all_pheno$S2[[ITER]]$h2_realized)),
  sprintf("S3: %.4f", unique(all_pheno$S3[[ITER]]$h2_realized)),
  sprintf("S4: %.4f", unique(all_pheno$S4[[ITER]]$h2_realized)),
  sprintf("S5: %.4f", unique(all_pheno$S5[[ITER]]$h2_realized)),
  sprintf("S6: %.4f", unique(all_pheno$S6[[ITER]]$h2_realized))
)
writeLines(summary_lines, file.path(ITER_DIR, "simulation_summary.txt"))

cat("Output written to", OUT, "\n")