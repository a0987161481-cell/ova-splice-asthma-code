## =====================================================================
##  lib_human.R --- cross-species mapping and human asthma validation
##  Caches everything to data/human/*.rds
## =====================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble)
})

HDIR   <- file.path(DATA, "human")
HCACHE <- file.path(HDIR, "human_analysis.rds")

## --------------------------------------------------- ortholog mapping
build_orthologs <- function() {
  f <- file.path(HDIR, "orthologs.rds")
  if (file.exists(f)) return(readRDS(f))
  suppressPackageStartupMessages(library(babelgene))
  ## mouse -> human orthologs for every mouse gene symbol we work with.
  ## babelgene rejects mixed identifier types, so keep only plausible symbols.
  q <- unique(c(DEG_ALL, DAS_GENES, DEG_RELAX, EV_DSI$gene,
                sf_state$gene, tf_class$gene, genes_all))
  q <- q[grepl("^[A-Za-z][A-Za-z0-9._-]*$", q) &
         !grepl("^(MSTRG|ENSMUS|ENSG|LOC|Gm[0-9]+$)", q)]
  bb <- babelgene::orthologs(genes = q, species = "mouse", human = FALSE)
  o  <- bb %>%
    transmute(mouse = symbol, human = human_symbol, support_n = support_n) %>%
    filter(!is.na(mouse), !is.na(human), human != "") %>%
    arrange(desc(support_n)) %>%
    distinct(mouse, .keep_all = TRUE)
  ## fallback: simple case-insensitive symbol identity for genes babelgene misses
  miss <- setdiff(q, o$mouse)
  fb <- tibble(mouse = miss, human = toupper(miss), support_n = 0)
  o  <- bind_rows(o, fb)
  saveRDS(o, f)
  cat("  orthologs:", nrow(o), "mapped (", nrow(bb), "from babelgene )\n")
  o
}

## ------------------------------------------- GSE152004 (n = 695, RNA-seq)
load_gse152004 <- function() {
  f <- file.path(HDIR, "gse152004_deg.rds")
  cnt <- read_tsv(file.path(HDIR, "GSE152004_counts.tsv"), show_col_types = FALSE)
  met <- read_tsv(file.path(HDIR, "GSE152004_meta.tsv"), show_col_types = FALSE)
  cnt <- cnt %>% filter(!is.na(gene), gene != "")
  cnt <- cnt[!duplicated(cnt$gene), ]
  M <- as.matrix(cnt[, met$sample]); rownames(M) <- cnt$gene
  keep <- rowSums(M >= 10) >= 0.2 * ncol(M)      # expressed filter
  M <- M[keep, , drop = FALSE]
  grp <- factor(met$group, levels = c("control", "asthma"))

  if (file.exists(f)) {
    res <- readRDS(f)
  } else {
    suppressPackageStartupMessages(library(DESeq2))
    cd <- DESeqDataSetFromMatrix(round(M), data.frame(g = grp), ~ g)
    cd <- DESeq(cd, quiet = TRUE)
    res <- as.data.frame(results(cd, contrast = c("g", "asthma", "control")))
    res$gene <- rownames(res)
    res <- as_tibble(res)
    saveRDS(res, f)
  }
  list(counts = M, meta = met, group = grp, res = res)
}

## -------------------------------------- GSE158752 (severity, RNA-seq)
load_gse158752 <- function() {
  f <- file.path(HDIR, "gse158752_deg.rds")
  cnt <- read_tsv(file.path(HDIR, "GSE158752_counts.tsv"), show_col_types = FALSE)
  met <- read_tsv(file.path(HDIR, "GSE158752_meta.tsv"), show_col_types = FALSE)
  cnt <- cnt %>% filter(!is.na(gene), gene != "")
  cnt <- cnt[!duplicated(cnt$gene), ]
  M <- as.matrix(cnt[, met$sample]); rownames(M) <- cnt$gene
  keep <- rowSums(M >= 10) >= 0.2 * ncol(M)
  M <- M[keep, , drop = FALSE]
  grp <- factor(met$group, levels = c("Healthy", "MildModerate", "Severe"))
  if (file.exists(f)) {
    res <- readRDS(f)
  } else {
    suppressPackageStartupMessages(library(DESeq2))
    cd <- DESeqDataSetFromMatrix(round(M), data.frame(g = grp), ~ g)
    cd <- DESeq(cd, quiet = TRUE)
    res <- as.data.frame(results(cd, contrast = c("g", "Severe", "Healthy")))
    res$gene <- rownames(res); res <- as_tibble(res)
    saveRDS(res, f)
  }
  list(counts = M, meta = met, group = grp, res = res)
}

## ------------------------------ microarray cohorts (severity / phenotype)
load_micro <- function(acc, case_pat, ctrl_pat) {
  f <- file.path(HDIR, paste0(tolower(acc), "_deg.rds"))
  ex <- read_tsv(file.path(HDIR, paste0(acc, "_expr.tsv")), show_col_types = FALSE)
  mt <- read_tsv(file.path(HDIR, paste0(acc, "_meta.tsv")), show_col_types = FALSE)
  ids <- ex[[1]]
  E <- as.matrix(ex[, -1]); rownames(E) <- ids
  E <- E[, mt$sample, drop = FALSE]
  grp <- ifelse(grepl(case_pat, mt$group, ignore.case = TRUE), "case",
         ifelse(grepl(ctrl_pat, mt$group, ignore.case = TRUE), "control", NA))
  ok <- !is.na(grp)
  E <- E[, ok, drop = FALSE]; grp <- grp[ok]
  E[E < 0] <- 0                       # microarray values may go slightly negative
  E <- log2(E + 1)
  if (file.exists(f)) {
    res <- readRDS(f)
  } else {
    r <- apply(E, 1, function(v) {
      a <- v[grp == "case"]; b <- v[grp == "control"]
      if (length(a) < 3 || length(b) < 3) return(c(NA, NA))
      tt <- suppressWarnings(t.test(a, b))
      c(mean(a) - mean(b), tt$p.value)
    })
    res <- tibble(gene = rownames(E), log2FC = r[1, ], pvalue = r[2, ]) %>%
      mutate(FDR = p.adjust(pvalue, "BH"))
    saveRDS(res, f)
  }
  list(expr = E, group = grp, res = res)
}

cat("lib_human.R loaded\n")
