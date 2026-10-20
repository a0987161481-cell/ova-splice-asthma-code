## =====================================================================
##  30_human_analysis.R
##  Cross-species mapping and human asthma validation of the mouse findings.
##  Caches every derived object to data/human/human_results.rds
## =====================================================================
.libPaths(c("D:/Rlibs", .libPaths()))
setwd("D:/学校工作/新建文件夹 (2)/OVA/manuscript_extension/R")
source("lib_theme.R"); source("lib_data.R")
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(purrr) })
source("lib_human.R")

CACHE <- file.path(HDIR, "human_results.rds")

if (file.exists(CACHE)) {
  HR <- readRDS(CACHE)
} else {
  orth <- build_orthologs()

  ## -------------------------------------------------- load human cohorts
  h1 <- load_gse152004()      # nasal epithelium, 441 asthma vs 254 control
  h2 <- load_gse158752()      # bronchial epithelium, graded severity
  g67 <- load_micro("GSE67472", "asthma", "healthy")
  g43 <- load_micro("GSE43696", "asthma", "control")

  ## ------------------------------------------- mouse -> human effect sizes
  mdeg <- deg_str %>% select(mouse = gene, log2FC_mouse = log2FC, FDR_mouse = FDR)
  mdasa <- tibble(mouse = DAS_GENES) %>%
    left_join(EV_HC %>% group_by(gene) %>% arrange(fdr) %>% slice(1) %>%
                ungroup() %>% select(mouse = gene, dpsi = dpsi, fdr_das = fdr),
              by = "mouse") %>%
    left_join(deg_all %>% select(mouse = gene, log2FC_all = log2FC,
                                 FDR_all = FDR, P_all = Pvalue), by = "mouse") %>%
    ## the reference direction for the splicing-layer comparison is the gene-level
    ## effect of the same gene, so that human expression can be compared with it
    mutate(log2FC_mouse = log2FC_all, is_das = TRUE)
  msf <- tibble(mouse = sf_state$gene[sf_state$in_expr],
                is_sf = TRUE, sf_de = sf_state$is_DE[sf_state$in_expr],
                sf_das = sf_state$is_DAS[sf_state$in_expr])
  mtf <- tibble(mouse = tf_class$gene, is_tf = TRUE)

  ## human differential expression per cohort
  hum <- list(
    GSE152004 = h1$res %>% transmute(human = gene, log2FC_h = log2FoldChange,
                                     pval_h = pvalue, padj_h = padj),
    GSE158752 = h2$res %>% transmute(human = gene, log2FC_h = log2FoldChange,
                                     pval_h = pvalue, padj_h = padj),
    GSE67472  = g67$res %>% transmute(human = gene, log2FC_h = log2FC,
                                      pval_h = pvalue, padj_h = FDR),
    GSE43696  = g43$res %>% transmute(human = gene, log2FC_h = log2FC,
                                      pval_h = pvalue, padj_h = FDR))

  ## ------------------------------------ concordance of the mouse DEG layer
  concord <- function(mtab, htab, label) {
    d <- mtab %>%
      inner_join(orth, by = "mouse") %>%
      inner_join(htab, by = "human") %>%
      filter(!is.na(log2FC_h), !is.na(log2FC_mouse))
    if (!nrow(d)) return(NULL)
    ct <- suppressWarnings(cor.test(d$log2FC_mouse, d$log2FC_h, method = "spearman"))
    ci <- boot_ci(sign(d$log2FC_mouse) == sign(d$log2FC_h), n = 5000)
    list(tab = d, n = nrow(d),
         concordant = sum(sign(d$log2FC_mouse) == sign(d$log2FC_h)),
         pct = 100 * ci["est"], lo = 100 * ci["lo"], hi = 100 * ci["hi"],
         rho = unname(ct$estimate), p = ct$p.value, label = label)
  }

  CONC <- list()
  CONC$deg_152004 <- concord(mdeg, hum$GSE152004, "Stringent DEGs / GSE152004")
  CONC$deg_158752 <- concord(mdeg, hum$GSE158752, "Stringent DEGs / GSE158752")
  CONC$deg_67472  <- concord(mdeg, hum$GSE67472,  "Stringent DEGs / GSE67472")
  CONC$deg_43696  <- concord(mdeg, hum$GSE43696,  "Stringent DEGs / GSE43696")
  CONC$das_152004 <- concord(mdasa, hum$GSE152004, "High-confidence DAS genes / GSE152004")
  CONC$das_158752 <- concord(mdasa, hum$GSE158752, "High-confidence DAS genes / GSE158752")
  CONC$das_67472  <- concord(mdasa, hum$GSE67472,  "High-confidence DAS genes / GSE67472")
  CONC$das_43696  <- concord(mdasa, hum$GSE43696,  "High-confidence DAS genes / GSE43696")

  ## ------------------------------- splicing-factor module in human asthma
  sf_h <- msf %>% inner_join(orth, by = "mouse") %>% inner_join(hum$GSE152004, by = "human")
  sf_h <- sf_h %>% filter(!is.na(log2FC_h), !is.na(padj_h))
  sf_summary <- list(
    n_mapped = nrow(sf_h),
    n_sig = sum(sf_h$padj_h < 0.05),
    n_sig_up = sum(sf_h$padj_h < 0.05 & sf_h$log2FC_h > 0),
    n_sig_dn = sum(sf_h$padj_h < 0.05 & sf_h$log2FC_h < 0),
    mean_lfc = mean(sf_h$log2FC_h, na.rm = TRUE),
    wilcox = suppressWarnings(wilcox.test(sf_h$log2FC_h, mu = 0)$p.value))
  ## background: all mapped genes
  bg_h <- orth %>% inner_join(hum$GSE152004, by = "human") %>% filter(!is.na(padj_h))
  sf_summary$bg_n_sig_pct <- 100 * mean(bg_h$padj_h < 0.05, na.rm = TRUE)
  sf_summary$sf_n_sig_pct <- 100 * mean(sf_h$padj_h < 0.05, na.rm = TRUE)
  fisher_sf <- fisher.test(matrix(c(sum(sf_h$padj_h < 0.05), sum(sf_h$padj_h >= 0.05),
                                    sum(bg_h$padj_h < 0.05), sum(bg_h$padj_h >= 0.05)),
                                  nrow = 2))
  sf_summary$fisher_p <- fisher_sf$p.value
  sf_summary$fisher_or <- unname(fisher_sf$estimate)

  ## --------------------------- severity gradient (GSE158752 and GSE43696)
  grad <- function(expr, grp, genes, target_level) {
    genes <- intersect(genes, rownames(expr))
    if (!length(genes)) return(NULL)
    sc <- rowMeans(scale(t(expr[genes, , drop = FALSE])))
    tibble(sample = names(sc), score = sc, group = grp) %>%
      group_by(group) %>%
      summarise(mean = mean(score), sd = sd(score), n = dplyr::n(), .groups = "drop") %>%
      mutate(target = target_level)
  }
  g158_genes <- msf %>% inner_join(orth, by = "mouse") %>%
    pull(human) %>% intersect(rownames(h2$counts))
  g158_l2 <- log2(sweep(h2$counts, 2, colSums(h2$counts), "/") * 1e6 + 1)
  sev_158 <- list(
    sf   = grad(g158_l2, as.character(h2$group), g158_genes, "SF/RBP"),
    deg  = grad(g158_l2, as.character(h2$group),
                mdeg %>% inner_join(orth, by = "mouse") %>% pull(human), "Mouse DEGs"),
    das  = grad(g158_l2, as.character(h2$group),
                mdasa %>% inner_join(orth, by = "mouse") %>% pull(human), "Mouse DAS genes"))
  ## trend test across the ordered severity levels
  sev_trend <- map_dfr(sev_158, function(x) {
    if (is.null(x)) return(NULL)
    x <- x %>% mutate(rank = match(group, c("Healthy", "MildModerate", "Severe")))
    ct <- suppressWarnings(cor.test(x$rank, x$mean, method = "spearman"))
    tibble(target = x$target[1], rho = unname(ct$estimate), p = ct$p.value)
  })

  ## GSE43696 severity (control / moderate / severe)
  sev_43696 <- NULL
  if (nrow(g43$res)) {
    grp43 <- read_tsv(file.path(HDIR, "GSE43696_meta.tsv"), show_col_types = FALSE)
    lab43 <- ifelse(grepl("Severe", grp43$group), "Severe",
             ifelse(grepl("Moderate", grp43$group), "MildModerate",
             ifelse(grepl("Control", grp43$group), "Healthy", NA)))
    ok <- !is.na(lab43)
    E43 <- g43$expr[, ok, drop = FALSE]
    lab43 <- lab43[ok]
    sf43 <- msf %>% inner_join(orth, by = "mouse") %>% pull(human) %>%
      intersect(rownames(E43))
    sc <- rowMeans(scale(t(E43[sf43, , drop = FALSE])))
    sev_43696 <- tibble(sample = names(sc), score = sc, group = lab43) %>%
      group_by(group) %>% summarise(mean = mean(score), sd = sd(score),
                                    n = dplyr::n(), .groups = "drop") %>%
      mutate(target = "SF/RBP")
  }

  ## ---------------------------------------------------- ROC / biomarker
  ## A signature score must come from per-sample expression; a per-gene log2FC
  ## is constant across samples and cannot be used as a classifier.
  cpm1 <- log2(sweep(h1$counts, 2, colSums(h1$counts), "/") * 1e6 + 1)
  sig_score <- function(weights, expr) {
    g <- intersect(names(weights), rownames(expr))
    if (length(g) < 5) return(NULL)
    w <- weights[g]
    list(genes = g, w = w,
         score = as.numeric(crossprod(w, expr[g, , drop = FALSE]) / sum(abs(w))))
  }
  suppressPackageStartupMessages(library(pROC))
  mk_w <- function(h, v = NULL) setNames(if (is.null(v)) rep(1, length(h)) else v, h)

  sigsets <- list(
    "Mouse stringent DEGs" = {
      d <- mdeg %>% inner_join(orth, by = "mouse")
      mk_w(d$human, d$log2FC_mouse)
    },
    "Mouse high-confidence DAS genes" = {
      d <- mdasa %>% inner_join(orth, by = "mouse") %>% filter(!is.na(log2FC_mouse))
      mk_w(d$human, sign(d$log2FC_mouse))
    },
    "Mouse splicing-factor / RBP module" = {
      d <- msf %>% inner_join(orth, by = "mouse") %>% distinct(human)
      mk_w(d$human)
    })
  set.seed(99)
  sigsets[["Random gene set (n = 300)"]] <-
    mk_w(sample(orth$human, 300))

  y <- as.integer(h1$meta$group == "asthma")
  ROC <- list()
  for (nm in names(sigsets)) {
    ss <- sig_score(sigsets[[nm]], cpm1)
    if (is.null(ss)) next
    r  <- pROC::roc(y, ss$score, quiet = TRUE)
    ci <- pROC::ci.auc(r)
    ROC[[nm]] <- list(label = nm, n_genes = length(ss$genes),
                      auc = as.numeric(pROC::auc(r)), lo = ci[1], hi = ci[3],
                      roc = tibble(fpr = 1 - r$specificities,
                                   tpr = r$sensitivities, label = nm),
                      score = tibble(sample = h1$meta$sample, score = ss$score,
                                     group = h1$meta$group, label = nm))
  }
  ## individual-gene AUC for candidate markers (mouse concordant + human significant)
  cand_pool <- mdeg %>% inner_join(orth, by = "mouse") %>%
    left_join(hum$GSE152004, by = "human") %>%
    filter(!is.na(padj_h), padj_h < 0.01,
           sign(log2FC_h) == sign(log2FC_mouse)) %>%
    arrange(desc(abs(log2FC_mouse)))
  gene_auc <- map_dfr(intersect(cand_pool$human, rownames(cpm1)), function(g) {
    r <- pROC::roc(y, cpm1[g, ], quiet = TRUE)
    tibble(human = g, auc = as.numeric(pROC::auc(r)))
  })
  gene_auc$mouse <- orth$mouse[match(gene_auc$human, orth$human)]
  gene_auc$log2FC_h <- hum$GSE152004$log2FC_h[match(gene_auc$human, hum$GSE152004$human)]
  gene_auc <- gene_auc %>% arrange(desc(auc))

  ## ------------------------------------------ cell-type association (mouse)
  LUNG_MARKERS <- list(
    "Alveolar type II" = c("Sftpc","Sftpb","Sftpa1","Abca3","Lamp3","Napsa","Slc34a2"),
    "Alveolar type I"  = c("Ager","Pdpn","Hopx","Aqp5","Cav1"),
    "Ciliated"         = c("Foxj1","Pcdh1","Tmem212","Dynlrb2","Ccdc153","Hydin"),
    "Club / secretory" = c("Scgb1a1","Scgb3a2","Bpifa1","Muc5b","Foxj1"),
    "Basal"            = c("Krt5","Krt14","Trp63","Krt15"),
    "Endothelial"      = c("Pecam1","Cdh5","Tek","Col4a1","Emcn","Vwf"),
    "Fibroblast"       = c("Col1a1","Col1a2","Col3a1","Pdgfra","Dcn","Lum"),
    "Smooth muscle"    = c("Acta2","Myh11","Tagln","Des","Cnn1"),
    "Macrophage"       = c("Adgre1","Cd68","Lyz2","Csf1r","Itgam","Mrc1"),
    "T cell"           = c("Cd3e","Cd3d","Cd3g","Cd4","Cd8a","Lck"),
    "B cell"           = c("Cd79a","Cd79b","Ms4a1","Cd19","Igkc"),
    "Neutrophil"       = c("S100a8","S100a9","Ly6g","Retnlg","Mpo"),
    "Eosinophil"       = c("Prg2","Epx","Mcpt8","Il5ra","Ccr3"),
    "Dendritic cell"   = c("Itgax","Cd209a","H2-Ab1","Cd74","Batf3"))
  cm_scores <- sapply(LUNG_MARKERS, function(g) {
    g <- intersect(g, rownames(fpkm_mat))
    if (length(g) < 3) return(rep(NA_real_, ncol(fpkm_mat)))
    rowMeans(scale(t(log2(fpkm_mat[g, , drop = FALSE] + 1))))
  })
  rownames(cm_scores) <- SAMPLES
  ## correlation of each cell-type score with the expression of DAS / SF genes
  das_h <- intersect(DAS_GENES, rownames(fpkm_mat))
  sf_h2 <- intersect(sf_state$gene[sf_state$in_expr], rownames(fpkm_mat))
  deg_h <- intersect(DEG_ALL, rownames(fpkm_mat))
  cell_assoc <- map_dfr(colnames(cm_scores), function(ct) {
    s <- cm_scores[, ct]
    f <- function(genes) {
      if (length(genes) < 3) return(NA_real_)
      suppressWarnings(cor(colMeans(log2(fpkm_mat[genes, , drop = FALSE] + 1)), s,
                           method = "spearman"))
    }
    tibble(cell_type = ct,
           DAS = f(das_h), SF = f(sf_h2), DEG = f(deg_h),
           all = suppressWarnings(cor(colMeans(log2(fpkm_mat + 1)), s,
                                      method = "spearman")))
  })

  ## ------------------------------------------- drug-target / druggability
  ASTHMA_TARGETS <- c(
    "NR3C1",                              # inhaled corticosteroids (glucocorticoid receptor)
    "ADRB2",                              # SABA / LABA
    "CYSLTR1", "CYSLTR2", "ALOX5", "ALOX5AP",   # leukotriene pathway
    "IL4R", "IL13", "IL5", "IL5RA", "TSLP", "TSLPR" = "TSLP",
    "FCER1A", "MS4A2", "IGHE",            # IgE axis
    "IL33", "IL1RL1",                     # IL-33 axis
    "PTGDR2", "HPGDS",                    # CRTh2
    "TNF", "TNFRSF1B",                    # anti-TNF
    "PDE4A", "PDE4B", "PDE4D",            # PDE4 inhibitors
    "HRH1", "CysLT" = "CYSLTR1",
    "SIGLEC8", "CCR3", "CCR4",            # eosinophil / chemokine axis
    "STAT6", "GATA3", "FOXP3",            # Th2 transcription factors
    "IL2RA", "MST1R", "KIT", "PDGFRA", "EGFR", "VEGFA")
  ASTHMA_TARGETS <- unique(unname(ASTHMA_TARGETS))
  tgt <- orth %>% filter(human %in% ASTHMA_TARGETS) %>%
    left_join(hum$GSE152004, by = "human") %>%
    mutate(in_mouse_deg = mouse %in% DEG_ALL,
           in_mouse_das = mouse %in% DAS_GENES,
           in_mouse_sf  = mouse %in% sf_state$gene[sf_state$in_expr],
           human_sig = !is.na(padj_h) & padj_h < 0.05)

  ## candidate prioritisation score for biomarker / target ranking
  cand <- mdeg %>%
    inner_join(orth, by = "mouse") %>%
    left_join(hum$GSE152004, by = "human") %>%
    left_join(mdasa %>% select(mouse, dpsi, is_das), by = "mouse") %>%
    mutate(is_das = replace_na(is_das, FALSE),
           human_sig = !is.na(padj_h) & padj_h < 0.05,
           concordant = !is.na(log2FC_h) & sign(log2FC_h) == sign(log2FC_mouse),
           score = abs(log2FC_mouse) * pmax(-log10(pmax(FDR_mouse, 1e-300)), 1) *
                   ifelse(human_sig & concordant, 2, 1))

  HR <- list(orth = orth, CONC = CONC, sf_summary = sf_summary, sf_h = sf_h,
             bg_h = bg_h, sev_158 = sev_158, sev_trend = sev_trend,
             sev_43696 = sev_43696, ROC = ROC, cell_assoc = cell_assoc,
             cm_scores = cm_scores, targets = tgt, candidates = cand,
             gene_auc = gene_auc,
             hum = hum, n_h1 = nrow(h1$counts), n_h2 = nrow(h2$counts),
             h1_res = h1$res, h2_res = h2$res,
             groups = list(h1 = h1$group, h2 = h2$group),
             deg_counts = c(asthma = sum(h1$group == "asthma"),
                            control = sum(h1$group == "control")),
             n_sig_human = sum(h1$res$padj < 0.05, na.rm = TRUE))
  saveRDS(HR, CACHE)
}

## ---------------------------------------------------------- console report
cat("\n=== cross-species concordance of the expression layer ===\n")
for (k in grep("^deg_", names(HR$CONC), value = TRUE)) {
  x <- HR$CONC[[k]]
  if (is.null(x)) next
  cat(sprintf("  %-38s n=%4d  concordant=%4d (%.1f%%, 95%% CI %.1f-%.1f)  rho=%+.3f  P=%.2e\n",
              x$label, x$n, x$concordant, x$pct, x$lo, x$hi, x$rho, x$p))
}
cat("\n=== cross-species concordance of the splicing layer ===\n")
for (k in grep("^das_", names(HR$CONC), value = TRUE)) {
  x <- HR$CONC[[k]]
  if (is.null(x)) next
  cat(sprintf("  %-38s n=%4d  concordant=%4d (%.1f%%, 95%% CI %.1f-%.1f)  rho=%+.3f  P=%.2e\n",
              x$label, x$n, x$concordant, x$pct, x$lo, x$hi, x$rho, x$p))
}
cat("\n=== splicing-factor / RBP module in human asthma (GSE152004) ===\n")
print(unlist(HR$sf_summary))
cat("\n=== severity trend (GSE158752: Healthy -> MildModerate -> Severe) ===\n")
print(as.data.frame(HR$sev_trend), row.names = FALSE)
cat("\n=== severity score by group (GSE158752) ===\n")
print(as.data.frame(bind_rows(HR$sev_158)), row.names = FALSE)
cat("\n=== ROC (GSE152004) ===\n")
for (nm in names(HR$ROC))
  cat(sprintf("  %-24s AUC = %.3f  (95%% CI %.3f-%.3f)\n", nm,
              HR$ROC[[nm]]$auc, HR$ROC[[nm]]$lo, HR$ROC[[nm]]$hi))
cat("\n=== cell-type association (Spearman with module expression) ===\n")
print(as.data.frame(HR$cell_assoc %>% arrange(desc(abs(DAS)))), row.names = FALSE)
cat("\n=== asthma drug targets found ===\n")
print(as.data.frame(HR$targets %>% select(human, mouse, log2FC_h, padj_h,
        in_mouse_deg, in_mouse_das, in_mouse_sf, human_sig)), row.names = FALSE)
cat("\nDONE\n")
