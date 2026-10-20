## =====================================================================
##  lib_data.R  --  load and derive every analysis object used by figures
## =====================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })
source(file.path(dirname(sys.frame(1)$ofile %||% "."), "lib_theme.R"), local = TRUE)

`%||%` <- function(a, b) if (is.null(a)) b else a

rd <- function(f) readr::read_tsv(file.path(DATA, f), show_col_types = FALSE,
                                  progress = FALSE)
num <- function(x) suppressWarnings(as.numeric(x))

## ---------------------------------------------------------------- 1. QC
manifest <- rd("sample_manifest.tsv") %>%
  mutate(group = factor(group, levels = GRP_LAB),
         sample = factor(sample, levels = SAMPLES))

mapped_long <- manifest %>%
  select(sample, group, total_reads, mapped_reads, uniq_reads, multi_reads,
         mapped_pct, uniq_pct, multi_pct) %>%
  pivot_longer(-c(sample, group),
               names_to = "metric", values_to = "value")

## ------------------------------------------------------- 2. expression
counts <- rd("counts_matrix.tsv")
fpkm   <- rd("fpkm_matrix.tsv")
genes_all    <- counts$gene
expr_mat     <- as.matrix(counts[, SAMPLES]);  rownames(expr_mat) <- genes_all
fpkm_mat     <- as.matrix(fpkm[, SAMPLES]);    rownames(fpkm_mat) <- genes_all

## ------------------------------------------------------------- 3. DEGs
deg_all    <- rd("deg_all_genes.tsv") %>% mutate(across(c(Pvalue, FDR, log2FC), num))
deg_str    <- rd("deg_stringent.tsv") %>% mutate(across(c(Pvalue, FDR, log2FC), num))
DEG_UP     <- deg_str %>% filter(regulated == "up")   %>% pull(gene)
DEG_DN     <- deg_str %>% filter(regulated == "down") %>% pull(gene)
DEG_ALL    <- deg_str$gene
## relaxed sensitivity criterion from the manuscript
deg_relax  <- deg_all %>% filter(FDR < 0.05, abs(log2FC) >= 1)
DEG_RELAX  <- deg_relax$gene

gse <- lapply(c(KEGG = "gsea_KEGG.tsv", BP = "gsea_BP.tsv",
                CC = "gsea_CC.tsv", MF = "gsea_MF.tsv"), function(f) {
  rd(f) %>% mutate(across(c(setSize, enrichmentScore, NES, pvalue, qvalue, rank), num))
})

## --------------------------------------- 4. alternative splicing (rMATS)
as_ev <- rd("as_events.tsv")
names(as_ev) <- sub("^psi_(Control-|AA-)", "PSI_\\1", names(as_ev))
as_ev <- as_ev %>%
  mutate(across(c(exon_len, inc_form_len, skip_form_len, pvalue, fdr, dpsi,
                  psi_ctrl_mean, psi_ova_mean, max_supp_ctrl, max_supp_ova), num),
         across(all_of(paste0("PSI_", SAMPLES)), num),
         class = factor(class, levels = CLS_LEVELS),
         gene  = gsub('"', '', gene)) %>%
  mutate(neg_log10_fdr = -log10(pmax(fdr, 1e-300)),
         abs_dpsi      = abs(dpsi))

as_psi <- rd("as_psi_long.tsv") %>%
  mutate(class = factor(class, levels = CLS_LEVELS),
         gene  = gsub('"', '', gene),
         group = factor(group, levels = GRP_LAB),
         sample = factor(sample, levels = SAMPLES),
         psi = num(psi), ijc = num(ijc), sjc = num(sjc), support = num(support))

## --- the validated prioritisation cascade -----------------------------
EV_TESTED   <- as_ev
EV_FDR      <- as_ev %>% filter(fdr < 0.05)
EV_DSI      <- EV_FDR %>% filter(abs_dpsi >= 0.10)
EV_HC       <- EV_DSI %>% filter(n_ge10_ctrl >= 3, n_ge10_ova >= 3)
EV_HC4      <- EV_DSI %>% filter(n_ge10_ctrl >= 4, n_ge10_ova >= 4)
DAS_GENES   <- unique(EV_HC$gene)
DAS_GENES_CAND <- unique(EV_DSI$gene)

## one representative event per gene = lowest FDR
rep_ev <- EV_HC %>% arrange(fdr) %>% group_by(gene) %>% slice(1) %>% ungroup()

cascade <- tibble::tibble(
  step = factor(c("Events tested", "FDR < 0.05", "|dPSI| >= 0.10",
                  "Junction-read support", ">=10 reads in 4/4 both groups"),
                levels = c("Events tested", "FDR < 0.05", "|dPSI| >= 0.10",
                           "Junction-read support", ">=10 reads in 4/4 both groups")),
  n = c(nrow(EV_TESTED), nrow(EV_FDR), nrow(EV_DSI), nrow(EV_HC), nrow(EV_HC4)))

## ------------------------------------ 5. sensitivity / robustness -----
## NOTE: reconstructed from the archived per-replicate IncLevel values, which
## are exactly the quantities rMATS uses for the group comparison.
psi_wide <- as_psi %>%
  select(event, sample, psi) %>%
  pivot_wider(names_from = sample, values_from = psi)

sens_core <- function(ev) {
  m  <- as.matrix(ev[, SAMPLES])
  ok_ctrl <- rowSums(!is.na(m[, CTRL])) > 0
  ok_ova  <- rowSums(!is.na(m[, OVAGRP])) > 0
  list(m = m, ok = ok_ctrl & ok_ova)
}

sens_loo <- function(ev) {
  m <- as.matrix(ev[, SAMPLES])
  out <- lapply(SAMPLES, function(s) {
    keep_c <- setdiff(CTRL, s); keep_o <- setdiff(OVAGRP, s)
    rowMeans(m[, keep_o, drop = FALSE], na.rm = TRUE) -
      rowMeans(m[, keep_c, drop = FALSE], na.rm = TRUE)
  })
  do.call(cbind, out)
}
sens_paired <- function(ev) {
  m <- as.matrix(ev[, SAMPLES]); res <- list()
  for (c1 in CTRL) for (o1 in OVAGRP) {
    keep_c <- setdiff(CTRL, c1); keep_o <- setdiff(OVAGRP, o1)
    res[[paste(c1, o1, sep = "|")]] <-
      rowMeans(m[, keep_o, drop = FALSE], na.rm = TRUE) -
      rowMeans(m[, keep_c, drop = FALSE], na.rm = TRUE)
  }
  do.call(cbind, res)
}

hc_psi <- EV_HC %>% select(event, class, gene, dpsi, fdr, abs_dpsi,
                           psi_ctrl_mean, psi_ova_mean,
                           n_ge10_ctrl, n_ge10_ova,
                           max_supp_ctrl, max_supp_ova,
                           all_of(paste0("PSI_", SAMPLES)))
names(hc_psi) <- sub("^PSI_", "", names(hc_psi))
loo  <- sens_loo(hc_psi);  colnames(loo)  <- paste0("drop_",  c(CTRL, OVAGRP))
pair <- sens_paired(hc_psi)

hc_psi <- hc_psi %>%
  mutate(
    loo_min   = apply(loo, 1, min, na.rm = TRUE),
    loo_max   = apply(loo, 1, max, na.rm = TRUE),
    loo_same  = sign(loo) == sign(dpsi) & abs(loo) >= 0.05,
    loo_robust = rowSums(loo_same) == ncol(loo),
    loo_range = loo_max - loo_min,
    pair_min  = apply(pair, 1, min, na.rm = TRUE),
    pair_max  = apply(pair, 1, max, na.rm = TRUE),
    pair_same = sign(pair) == sign(dpsi) & abs(pair) >= 0.05,
    pair_robust = rowSums(pair_same) == ncol(pair),
    pair_range  = pair_max - pair_min,
    loo_shift   = (loo_min + loo_max) / 2 - dpsi
  )
N_LOO_OK  <- sum(hc_psi$loo_robust)
N_PAIR_OK <- sum(hc_psi$pair_robust)

## ------------------------------------------- 6. representation quality
## "measurable" = event had at least one non-missing PSI in every replicate
psi_meas <- as_psi %>% filter(!is.na(psi)) %>% count(event, sample) %>% count(event)
EV_MEASURABLE <- psi_meas$event[psi_meas$n == 8]

## ----------------------------------------- 7. splicing factor / RBP set
SF_RBP <- c(
  ## core spliceosome / snRNP
  "Snrpa","Snrpa1","Snrpb","Snrpb2","Snrpc","Snrpd1","Snrpd2","Snrpd3","Snrpe",
  "Snrpf","Snrpg","Snrnp70","Snrnp200","Snrnp25","Snrnp27","Snrnp35","Snrnp40",
  "Snrnp48","Snrnp70","Sart1","Sart3","Prpf3","Prpf4","Prpf4b","Prpf6","Prpf8",
  "Prpf19","Prpf31","Prpf38a","Prpf38b","Prpf39","Prpf40a","Prpf40b","PrpF4",
  "Phf5a","Bud31","Cdc40","Cwc15","Cwc22","Cwc25","Cwc27","Dhx8","Dhx15","Dhx16",
  "Dhx38","Ddx5","Ddx17","Ddx23","Ddx39b","Ddx41","Ddx42","Ddx46","Eftud2","Eif4a3",
  "Magoh","Magohb","Rbm8a","Casc3","Acin1","Ncbp1","Ncbp2","Ncbp3","NcBP",
  ## SR proteins and SR-related
  "Srsf1","Srsf2","Srsf3","Srsf4","Srsf5","Srsf6","Srsf7","Srsf9","Srsf10",
  "Srsf11","Srsf12","Srsf18","SrrM1","SrrM2","Tra2a","Tra2b","Sf1","U2af1","U2af2",
  "U2af1l4","U2surp","Rbfox1","Rbfox2","Rbfox3",
  ## hnRNP family
  "Hnrnpa0","Hnrnpa1","Hnrnpa2b1","Hnrnpa3","Hnrnpab","Hnrnpac","Hnrnpc","Hnrnph1",
  "Hnrnph2","Hnrnph3","Hnrnpk","Hnrnpl","Hnrnpll","HnrnpM","Hnrnpm","HnrnpU","Hnrnpr",
  "Hnrnpu","Hnrnpul1","Hnrnpul2","Hnrnpc",
  ## CELF / ELAV / KHDRBS / PTBP / others
  "Celf1","Celf2","Celf3","Celf4","Celf5","Celf6","Elavl1","Elavl2","Elavl3","Elavl4",
  "Khdrbs1","Khdrbs2","Khdrbs3","Ptbp1","Ptbp2","Ptbp3","Nova1","Nova2","Qki","Mbnl1",
  "Mbnl2","Mbnl3","Esrp1","Esrp2","Rbm3","Rbm4","Rbm4b","Rbm5","Rbm6","Rbm7","Rbm10",
  "Rbm11","Rbm12","Rbm14","Rbm15","Rbm15b","Rbm17","Rbm18","Rbm19","Rbm20","Rbm22",
  "Rbm23","Rbm24","Rbm25","Rbm26","Rbm27","Rbm28","Rbm31","Rbm33","Rbm34","Rbm38",
  "Rbm39","Rbm41","Rbm42","Rbm43","Rbm44","Rbm45","Rbm46","Rbm47","Rbm48","Rbmx",
  "Rbmxl1","Rbmxl2","Rbpms","Rbpms2","Rc3h1","Rc3h2",
  ## ZFP36 / AU-rich / stability
  "Zfp36","Zfp36l1","Zfp36l2","Zc3h7a","Zc3h7b","Zc3h11a","Zc3h12a","Zc3h13",
  "Zc3h14","Zc3h15","Zc3h18","Zc3hav1",
  ## DDX / DHX helicases with splicing roles
  "Ddx1","Ddx3x","Ddx3y","Ddx6","Ddx20","Ddx24","Ddx27","Ddx39a","Ddx50","Ddx51",
  "Ddx54","Ddx55","Ddx56","Ddx59","Dhx9","Dhx29","Dhx30","Dhx32","Dhx33","Dhx34",
  "Dhx35","Dhx36","Dhx37","Dhx40","Dhx57","Ddx19a","Ddx19b","Eif4a1","Eif4a2",
  ## exon-junction complex & NMD
  "Upf1","Upf2","Upf3a","Upf3b","Smg1","Smg5","Smg6","Smg7","Smg8","Smg9",
  ## cleavage / polyadenylation
  "Cpsf1","Cpsf2","Cpsf3","Cpsf4","Cpsf6","Cpsf7","Cstf1","Cstf2","Cstf3","Wdr33",
  "Fip1l1","Nudt21","Pcf11","Clp1","Sympk","Cpsf4l",
  ## editing / modification
  "Adar","Adarb1","Adarb2","Jmjd6","Mettl3","Mettl14","Wtap","Ythdc1","Ythdf1",
  "Ythdf2","Ythdf3","Fto","Alkbh5","Nsun2","Rbm15b",
  ## RNA helicase / chaperone
  "Hspa8","Ybx1","Ybx3","Csde1","Igf2bp1","Igf2bp2","Igf2bp3","Lin28a","Lin28b",
  ## specific pulmonary / inflammation-relevant RNA regulators
  "Zc3h12a","Mbnl1","Ncl","Npm1","Hnrnpd","Sfpq","Nono","Pspc1","Matr3","Tia1",
  "Tial1","Pabpc1","Pabpn1","Pum1","Pum2","Stau1","Stau2","Fmr1","Fxr1","Fxr2",
  "Dazap1","Dazap2","Raver1","Raver2","Cpeb1","Cpeb2","Cpeb3","Cpeb4","Zfp598",
  "Srek1","Srek1ip1","Son","Rnps1","Srsf8","Srsf7","Lsm1","Lsm2","Lsm3","Lsm4",
  "Lsm5","Lsm6","Lsm7","Lsm8","Lsm10","Lsm11","Lsm12","Lsm14a","Lsm14b","Gemin2",
  "GemIn5","Gemin5","Gemin6","Gemin7","Gemin8","Ddx20","Smn1","Fmr1"
)
SF_RBP <- unique(SF_RBP)

sf_state <- tibble(gene = SF_RBP) %>%
  mutate(
    in_expr   = gene %in% genes_all,
    is_DE     = gene %in% DEG_ALL,
    is_DE_rel = gene %in% DEG_RELAX,
    n_das     = sapply(gene, function(g) sum(EV_HC$gene == g)),
    is_DAS    = gene %in% DAS_GENES,
    n_cand    = sapply(gene, function(g) sum(EV_DSI$gene == g)),
    is_cand   = gene %in% DAS_GENES_CAND
  ) %>%
  left_join(deg_all %>% select(gene, log2FC, FDR) %>% rename(deg_log2FC = log2FC,
                                                             deg_FDR = FDR),
            by = "gene") %>%
  left_join(EV_HC %>% group_by(gene) %>% arrange(fdr) %>% slice(1) %>%
              ungroup() %>% select(gene, class, dpsi, fdr),
            by = "gene") %>%
  rename(das_class = class, das_dpsi = dpsi, das_fdr = fdr)

## ------------------------------------------------ 8. TF annotation -----
tf_class <- rd("tf_classification.tsv") %>%
  mutate(gene = gsub('"', '', gene),
         is_DE  = gene %in% DEG_ALL,
         is_DAS = gene %in% DAS_GENES,
         is_DE_rel = gene %in% DEG_RELAX) %>%
  left_join(deg_all %>% select(gene, log2FC, FDR) %>%
              rename(deg_log2FC = log2FC, deg_FDR = FDR), by = "gene") %>%
  left_join(EV_HC %>% group_by(gene) %>% arrange(fdr) %>% slice(1) %>%
              ungroup() %>% select(gene, class, dpsi, fdr) %>%
              rename(das_class = class, das_dpsi = dpsi, das_fdr = fdr),
            by = "gene")

## --------------------------------------------- 9. variants ------------
snp_anno    <- rd("snp_anno_stat.tsv")
indel_anno  <- rd("indel_anno_stat.tsv")
snp_sample  <- rd("snp_per_sample.tsv") %>%
  mutate(across(-sample, num), group = ifelse(grepl("^AA", sample), "OVA", "Control"))
snp_type    <- rd("snp_type_stat.tsv")

## --------------------------------------- 10. novel genes / structure ---
novel_genes <- rd("novel_genes.tsv")
gene_opt    <- rd("gene_structure_optimize.tsv")

## ------------------------------------------------ 11. PPI -------------
ppi_edges <- rd("ppi_edges_all.tsv")

## ---------------------------------------- 12. external validation -----
ext_path <- file.path(DATA, "external", "GSE295506_gene_counts.txt.gz")
have_external <- file.exists(ext_path)
if (have_external) {
  ext_raw <- readr::read_tsv(ext_path, show_col_types = FALSE, progress = FALSE) %>%
    select(-any_of(c(" ", ""))) %>% filter(!is.na(Gene), Gene != "")
  ext_cols <- grep("\\.Count$", colnames(ext_raw), value = TRUE)
  ext_raw <- ext_raw %>% filter(!duplicated(Gene))
  ext_mat <- as.matrix(ext_raw[, ext_cols]); rownames(ext_mat) <- ext_raw$Gene
  colnames(ext_mat) <- sub("\\.Count$", "", ext_cols)
  ext_grp <- ifelse(grepl("^sham", colnames(ext_mat)), "Control", "OVA")
  ## median-ratio normalisation (as described in the manuscript)
  lg <- log(pmax(ext_mat, 1))
  ref <- apply(lg, 1, function(v) { v <- v[is.finite(v)]; if (!length(v)) NA else
                                    exp(mean(v)) })
  ref <- ref[is.finite(ref) & ref > 0]
  sf  <- apply(sweep(ext_mat[names(ref), , drop = FALSE], 1, ref, "/"),
               2, function(v) median(v, na.rm = TRUE))
  ext_norm <- sweep(ext_mat, 2, sf, "/")
  ext_l2   <- log2(ext_norm + 1)
  ext_mean <- data.frame(
    gene = rownames(ext_l2),
    ext_log2FC = rowMeans(ext_l2[, ext_grp == "OVA", drop = FALSE], na.rm = TRUE) -
                 rowMeans(ext_l2[, ext_grp == "Control", drop = FALSE], na.rm = TRUE),
    ext_mean_expr = rowMeans(ext_l2, na.rm = TRUE), stringsAsFactors = FALSE)

  ext_cmp <- deg_str %>%
    select(gene, log2FC, FDR, regulated) %>%
    inner_join(ext_mean, by = "gene") %>%
    mutate(concordant = sign(log2FC) == sign(ext_log2FC),
           ext_abs    = abs(ext_log2FC),
           disc       = abs(log2FC) - abs(ext_log2FC))
  EXT_N_MEAS <- nrow(ext_cmp)
  EXT_N_CONC <- sum(ext_cmp$concordant)
  EXT_SUB    <- ext_cmp %>% filter(ext_abs >= 0.5)
  EXT_N_SUB  <- nrow(EXT_SUB)
  EXT_N_SUBC <- sum(EXT_SUB$concordant)
  ext_rho <- suppressWarnings(cor.test(ext_cmp$log2FC, ext_cmp$ext_log2FC,
                                       method = "spearman"))
  EXT_RHO <- unname(ext_rho$estimate); EXT_P <- ext_rho$p.value
}

cat("lib_data.R loaded:",
    nrow(deg_str), "DEGs;", nrow(EV_HC), "high-confidence DAS events;",
    length(DAS_GENES), "DAS genes;", N_LOO_OK, "LOO-robust;",
    N_PAIR_OK, "paired-robust\n")
