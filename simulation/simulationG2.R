suppressPackageStartupMessages({
  library(AlphaSimR)
  library(optiSel)
  library(tidyverse)
})

# ── Config ────────────────────────────────────────────────────────────────────
BASE_G1     <- "/scratch/user/aguswibowo/Research/simulation_forward/output_G1"
FIT_DIR     <- "/scratch/user/aguswibowo/Research/predict_salmon_forward/output/fitG1"
OUT         <- "/QRISdata/Q8514/research_data/simulation_forward/output_G2"
RDS_DIR     <- file.path(BASE_G1, "rds")
H2          <- 0.3
SEL_PROP    <- 0.30
N_PROG      <- 25
N_PAIRS     <- 120
DELTA_F     <- 0.01
BREED_NAME  <- "salmon"

MARKER_TYPES <- c("snp", "mh")
MODELS       <- c("gblup", "bayesR", "bayesA")
SCENARIOS    <- c("S1_equal", "S2_unequal", "S3_equal", "S4_unequal", "S5_equal", "S6_unequal")

args    <- commandArgs(trailingOnly = TRUE)
ITER    <- as.integer(args[1])
ITERPAD <- sprintf("%02d", ITER)

# ── Load shared objects ───────────────────────────────────────────────────────
SP             <- readRDS(file.path(RDS_DIR, "SP.rds"))
snp_meta       <- readRDS(file.path(RDS_DIR, "snp_meta.rds"))
qtl_s1         <- readRDS(file.path(RDS_DIR, "qtl_s1.rds"))
qtl_s2         <- readRDS(file.path(RDS_DIR, "qtl_s2.rds"))
qtl_s3_effects <- readRDS(file.path(RDS_DIR, "qtl_s3_effects.rds"))
qtl_s4_effects <- readRDS(file.path(RDS_DIR, "qtl_s4_effects.rds"))
qtl_s5_effects <- readRDS(file.path(RDS_DIR, "qtl_s5_effects.rds"))
qtl_s6_effects <- readRDS(file.path(RDS_DIR, "qtl_s6_effects.rds"))
qtl_mh_blocks  <- readRDS(file.path(RDS_DIR, "qtl_mh_blocks.rds"))
qtl_hap_blocks <- readRDS(file.path(RDS_DIR, "qtl_hap_blocks.rds"))

G1_pop <- readRDS(file.path(BASE_G1, sprintf("iter%s", ITERPAD), "G1_pop.rds"))

pedigree_g1 <- read_csv(
  file.path(BASE_G1, sprintf("iter%s", ITERPAD), "pedigree_G1.csv"),
  show_col_types = FALSE
)
g1_ind_ids <- pedigree_g1$individual_id

# ── TBV functions ─────────────────────────────────────────────────────────────
calc_tbv_snp <- function(dosage_mat, qtl_effects) {
  cols <- intersect(qtl_effects$snp_id, colnames(dosage_mat))
  eff  <- qtl_effects$effect[match(cols, qtl_effects$snp_id)]
  tbv  <- as.matrix(dosage_mat[, cols]) %*% eff
  as.vector(tbv - mean(tbv))
}

calc_tbv_block <- function(ind_ids, block_ids, eff_list, block_snp_map, G2_haplo) {
  by_chr <- block_snp_map %>%
    filter(block_id %in% block_ids) %>%
    group_by(block_id) %>%
    summarise(chr = first(chr), snp_indices = list(snp_index), .groups = "drop")

  n_ind <- length(ind_ids)
  tbv   <- numeric(n_ind)

  for (bi in seq_len(nrow(by_chr))) {
    blk     <- by_chr$block_id[bi]
    c_val   <- by_chr$chr[bi]
    snp_idx <- unlist(by_chr$snp_indices[bi]) + 1
    haplo   <- G2_haplo[[c_val]]
    eff     <- eff_list[[blk]]

    for (i in seq_len(n_ind)) {
      h1 <- haplo[(i - 1) * 2 + 1, snp_idx]
      h2 <- haplo[(i - 1) * 2 + 2, snp_idx]
      a1 <- strtoi(paste(h1, collapse = ""), base = 2L) + 1
      a2 <- strtoi(paste(h2, collapse = ""), base = 2L) + 1
      e1 <- if (!is.na(eff[as.character(a1)])) eff[as.character(a1)] else 0
      e2 <- if (!is.na(eff[as.character(a2)])) eff[as.character(a2)] else 0
      tbv[i] <- tbv[i] + e1 + e2
    }
  }
  tbv - mean(tbv)
}

rescale_tbv <- function(tbv) (tbv - mean(tbv)) / sd(tbv)

sim_pheno <- function(tbv, h2 = H2) {
  var_e <- var(tbv) * (1 - h2) / h2
  tbv + rnorm(length(tbv), 0, sqrt(var_e))
}

make_pheno_df <- function(tbv, ind_ids, iter, fam_ids, h2 = H2) {
  pheno   <- sim_pheno(tbv, h2)
  h2_real <- var(tbv) / var(pheno)
  tibble(
    individual_id   = ind_ids,
    tbv             = tbv,
    phenotype       = pheno,
    generation      = 2,
    population_type = "candidate",
    iteration       = iter,
    family_id       = fam_ids,
    h2_realized     = h2_real
  )
}

# ── Block SNP map ─────────────────────────────────────────────────────────────
BASE_RAW <- "/scratch/user/aguswibowo/Research/simulation_forward/raw_data"
mh_snps  <- read_csv(file.path(BASE_RAW, "mh_info_ld_haploblock_G0/stats/snp_selection_detailed.csv"),
                     show_col_types = FALSE)
hap_snps <- read_csv(file.path(BASE_RAW, "all_info_haploblock_G0/stats/snp_selection_detailed.csv"),
                     show_col_types = FALSE)

n_chr <- length(unique(snp_meta$chr))

# GRM
dosage_g1        <- pullSnpGeno(G1_pop, simParam = SP)
rownames(dosage_g1) <- g1_ind_ids
freq_g1          <- colMeans(dosage_g1) / 2
keep_snp_g1      <- freq_g1 > 0.01 & freq_g1 < 0.99
Z_g1             <- sweep(dosage_g1[, keep_snp_g1], 2, 2 * freq_g1[keep_snp_g1], "-")
k_g1             <- 2 * sum(freq_g1[keep_snp_g1] * (1 - freq_g1[keep_snp_g1]))
Gmat_g1          <- tcrossprod(Z_g1) / k_g1
Kin_g1           <- Gmat_g1 / 2
rownames(Kin_g1) <- colnames(Kin_g1) <- g1_ind_ids

set.seed(ITER * 999)
fam_sex <- pedigree_g1 %>%
        distinct(family_id) %>%
        mutate(fam_sex = sample(rep(c("male", "female"), length.out = n())))

# ── Main loop ─────────────────────────────────────────────────────────────────
for (marker in MARKER_TYPES) {
  for (model in MODELS) {
    for (sc in SCENARIOS) {

      combo_tag <- sprintf("%s_%s_%s", marker, model, sc)
      fit_path  <- file.path(FIT_DIR, sprintf("iter%s", ITERPAD),
                             marker, sc, model, "gebv_ranked.rds")

      if (!file.exists(fit_path)) {
        cat(sprintf("SKIP %s — gebv_ranked.rds not found\n", combo_tag))
        next
      }

      set.seed(ITER * 1000 + which(MARKER_TYPES == marker) * 100 +
               which(MODELS == model) * 10 + which(SCENARIOS == sc))

      out_dir <- file.path(OUT, sprintf("iter%s", ITERPAD), combo_tag)
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

      # ── OCS: optimum contribution selection ──
      gebv_df <- readRDS(fit_path)

      # All cand
      cand <- pedigree_g1 %>%
        left_join(fam_sex, by = "family_id") %>%
        mutate(
            GEBV  = gebv_df$GEBV[match(individual_id, gebv_df$individual_id)],
            Breed = BREED_NAME
        ) %>%
        filter(!is.na(GEBV)) %>%
        select(Indiv = individual_id, Sex = fam_sex, Breed, GEBV) %>%
        as.data.frame()

      Kin_cand         <- Kin_g1[cand$Indiv, cand$Indiv]
      n_cand           <- nrow(Kin_cand)
      mean_kin_current <- (sum(Kin_cand) - sum(diag(Kin_cand))) / (n_cand * (n_cand - 1))
      ub_kin           <- mean_kin_current + DELTA_F
    
      ocs <- tryCatch({
        cand_obj <- candes(phen = cand, sKin = Kin_cand, quiet = TRUE)
        opticont("max.GEBV", cand_obj,
            con   = list(ub.sKin = ub_kin),
            quiet = TRUE)
      }, error = function(e) {
        cat(sprintf("  OCS failed (%s), fallback ke truncation\n", e$message))
        NULL
      })
    
      if (!is.null(ocs)) {
        candidate <- ocs$parent

        # Stage 1: Hard cutoff top N_PAIRS based on oc
        N_POOL  <- N_PAIRS * 2 
        sires <- candidate %>% filter(Sex == "male")   %>%
                 arrange(desc(oc)) %>% slice_head(n = N_POOL)
        dams  <- candidate %>% filter(Sex == "female") %>%
                 arrange(desc(oc)) %>% slice_head(n = N_POOL)

        # Stage 2: matings() minimise kinship within selected pool
        pool       <- bind_rows(sires, dams)
        Kin_pool   <- Kin_cand[pool$Indiv, pool$Indiv]
        pool_df    <- data.frame(Indiv = pool$Indiv,
                                 Sex   = pool$Sex,
                                 Breed = BREED_NAME,
                                 oc    = pool$oc)
        pool_df$n  <- noffspring(pool_df, N = N_PAIRS * N_PROG)$nOff
        cross_plan <- matings(pool_df, Kin_pool)

        if (nrow(cross_plan) < N_PAIRS)
          warning(sprintf("  matings() returned %d pairs < N_PAIRS=%d — n_G2 akan < 3750",
                          nrow(cross_plan), N_PAIRS))

        cross_plan <- cross_plan %>% slice_head(n = N_PAIRS)

        # Log 1: info parents
        cat(sprintf("  Pool: %d sires | %d dams | Selected top %d each\n",
                    nrow(sires), nrow(dams), N_PAIRS))
        cat(sprintf("  Mean GEBV sires=%.4f | dams=%.4f\n",
                    mean(sires$GEBV[1:N_PAIRS], na.rm=TRUE),
                    mean(dams$GEBV[1:N_PAIRS],  na.rm=TRUE)))

        # Log 2: kinship 
        pair_kin <- mapply(function(s, d) Kin_pool[s, d],
                           cross_plan$Sire, cross_plan$Dam)
        cat(sprintf("  Pair kinship: mean=%.4f | max=%.4f\n",
                    mean(pair_kin), max(pair_kin)))

        sire_idx  <- match(cross_plan$Sire, g1_ind_ids)
        dam_idx   <- match(cross_plan$Dam,  g1_ind_ids)
        n_pairs   <- nrow(cross_plan)
        n_off_vec <- rep(N_PROG, n_pairs)

        cross_matrix <- cbind(sire_idx, dam_idx)
        G2_list <- lapply(seq_len(n_pairs), function(i) {
          makeCross(G1_pop,
                    crossPlan = matrix(cross_matrix[i, ], nrow = 1),
                    nProgeny  = N_PROG,
                    simParam  = SP)
        })
        G2 <- mergePops(G2_list)
        cat(sprintf("  OCS: %d crosses | mean_kin_before=%.4f | ub_kin=%.4f\n",
                    n_pairs, mean_kin_current, ub_kin))
        cat(sprintf("  G2 produced: %d families x %d offspring = %d individuals\n",
                    n_pairs, N_PROG, G2@nInd))
        } else {
          # Fallback: truncation selection + randCross2
          n_sel    <- max(2, round(length(g1_ind_ids) * SEL_PROP))
          sel_ids  <- gebv_df %>% arrange(rank) %>% slice_head(n = n_sel) %>% pull(individual_id)
          sel_idx  <- which(g1_ind_ids %in% sel_ids)
          G1_sel   <- G1_pop[sel_idx]
          n_half   <- floor(G1_sel@nInd / 2)
          sires    <- G1_sel[1:n_half]
          dams     <- G1_sel[(n_half + 1):G1_sel@nInd]
          n_pairs  <- n_half
          n_off_vec <- rep(N_PROG, n_pairs)
          G2       <- randCross2(sires, dams,
                                nCrosses  = n_pairs,
                                nProgeny  = N_PROG,
                                simParam  = SP)
        }

      n_G2    <- G2@nInd
      fam_ids <- rep(seq_len(n_pairs), times = n_off_vec)
      ind_ids <- sprintf("iter%d_gen2_%s_fam%03d_ind%02d",
                         ITER, combo_tag,
                         fam_ids, sequence(n_off_vec))

      cat(sprintf("iter%s | %s | n_pairs=%d | n_G2=%d\n",
                  ITERPAD, combo_tag, n_pairs, n_G2))

      # ── Pedigree G2 ──
      pedigree_g2 <- tibble(
        individual_id   = ind_ids,
        generation      = 2,
        population_type = "candidate",
        father_id       = NA_character_,
        mother_id       = NA_character_,
        iteration       = ITER,
        family_id       = fam_ids,
        marker_type     = marker,
        model           = model,
        scenario_fit    = sc
      )
      write_csv(pedigree_g2, file.path(out_dir, "pedigree_G2.csv"))

      # ── SNP dosage ──
      dosage <- pullSnpGeno(G2, simParam = SP)
      colnames(dosage) <- snp_meta$snp_id
      rownames(dosage) <- ind_ids

      geno_df <- as_tibble(dosage, rownames = "ID") %>%
        mutate(generation = 2, population_type = "candidate", .after = "ID")
      write_csv(geno_df, file.path(out_dir, "geno_G2.csv"))

      # ── Haplotype G2 per chr ──
      h        <- pullSegSiteHaplo(G2, simParam = SP)
      G2_haplo <- lapply(1:n_chr, function(c) {
        snp_in_chr <- which(snp_meta$chr == c)
        h[, snp_in_chr, drop = FALSE]
      })

      # ── VCF for phasing ──
      vcf_path <- file.path(out_dir, "geno_G2.vcf")
      writeLines(
        c("##fileformat=VCFv4.2",
          paste0("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t",
                 paste(ind_ids, collapse="\t"))),
        vcf_path
      )

      # ── Haplotype output per chr ──
      haplo_dir <- file.path(out_dir, "haploG2")
      dir.create(haplo_dir, recursive = TRUE, showWarnings = FALSE)
      for (c in 1:n_chr) {
        haplo_c <- G2_haplo[[c]]
        write.table(haplo_c, file.path(haplo_dir, sprintf("chr%d", c)),
                    row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
      }

      # ── TBV all scenario ──
      tbv_list <- list(
        S1_equal   = rescale_tbv(calc_tbv_snp(dosage, qtl_s1)),
        S2_unequal = rescale_tbv(calc_tbv_snp(dosage, qtl_s2)),
        S3_equal   = rescale_tbv(calc_tbv_block(ind_ids, qtl_mh_blocks,  qtl_s3_effects, mh_snps,  G2_haplo)),
        S4_unequal = rescale_tbv(calc_tbv_block(ind_ids, qtl_mh_blocks,  qtl_s4_effects, mh_snps,  G2_haplo)),
        S5_equal   = rescale_tbv(calc_tbv_block(ind_ids, qtl_hap_blocks, qtl_s5_effects, hap_snps, G2_haplo)),
        S6_unequal = rescale_tbv(calc_tbv_block(ind_ids, qtl_hap_blocks, qtl_s6_effects, hap_snps, G2_haplo))
      )

      pheno_list <- lapply(tbv_list, function(tbv) make_pheno_df(tbv, ind_ids, ITER, fam_ids))
      for (sc_pheno in names(pheno_list))
        write_csv(pheno_list[[sc_pheno]],
                  file.path(out_dir, sprintf("pheno_%s_G2.csv", sc_pheno)))
      
      # Log 3: h2_realized per scenario
      cat(sprintf("  h2_realized: %s\n",
                  paste(sprintf("%s=%.3f", names(pheno_list),
                                sapply(pheno_list, function(df) unique(df$h2_realized))),
                        collapse=" | ")))

      # Log 4: realized gain 
      if (!is.null(ocs)) {
        mean_gebv_parents <- mean(c(
          head(sires$GEBV, N_PAIRS),
          head(dams$GEBV,  N_PAIRS)
        ), na.rm = TRUE)
        mean_gebv_all <- mean(cand$GEBV, na.rm = TRUE)
        cat(sprintf("  Mean GEBV all G1=%.4f | selected parents=%.4f | selection diff=%.4f\n",
                    mean_gebv_all, mean_gebv_parents,
                    mean_gebv_parents - mean_gebv_all))
      }

      saveRDS(G2, file.path(out_dir, "G2_pop.rds"))

      writeLines(c(
        sprintf("Date        : %s", Sys.time()),
        sprintf("ITER        : %d", ITER),
        sprintf("Combo       : %s", combo_tag),
        sprintf("N_pairs     : %d", n_pairs),
        sprintf("N_prog      : %d", N_PROG),
        sprintf("N_G2        : %d", n_G2),
        paste(sprintf("h2_realized %s: %.4f",
                      names(pheno_list),
                      sapply(pheno_list, function(df) unique(df$h2_realized))))
      ), file.path(out_dir, "simulation_summary.txt"))

    }
  }
}

cat(sprintf("=== ITER %d complete ===\n", ITER))