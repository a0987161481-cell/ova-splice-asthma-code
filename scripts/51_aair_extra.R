## =====================================================================
##  51_aair_extra.R  (v2)
##  AAIR revision statistics:
##   1  DeLong comparison of the two mouse-derived signatures (GSE152004)
##   2  bootstrap CI of delta-AUC
##   3  replication of both AUCs in GSE158752 (independent RNA-seq cohort)
##   4  type-2-high / type-2-low stratification (GSE152004)
##  Writes data/human/aair_extra.rds and out/aair_extra_report.txt
## =====================================================================
.libPaths(c("D:/Rlibs", .libPaths()))
setwd("D:/学校工作/新建文件夹 (2)/OVA/manuscript_extension/R")
source("lib_theme.R"); source("lib_data.R")
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(purrr)
                                 library(pROC) })
source("lib_human.R")

REP <- file.path("..", "out", "aair_extra_report.txt")
say <- function(...) {
  m <- paste0(...); cat(m, "\n", sep = "")
  cat(m, "\n", sep = "", file = REP, append = TRUE)
}
try_catch <- function(expr, tag) {
  tryCatch(expr, error = function(e) say("  !! ", tag, " failed: ", conditionMessage(e)))
}
unlink(REP)

HR   <- readRDS(file.path(HDIR, "human_results.rds"))
orth <- HR$orth

## ---------------------------------------------------- signature weights
mdeg  <- deg_str %>% select(mouse = gene, log2FC_mouse = log2FC)
mdasa <- tibble(mouse = DAS_GENES) %>%
  left_join(deg_all %>% select(mouse = gene, log2FC_all = log2FC), by = "mouse") %>%
  mutate(log2FC_mouse = log2FC_all)
msf <- tibble(mouse = sf_state$gene[sf_state$in_expr], is_sf = TRUE)

w_deg <- { d <- mdeg  %>% inner_join(orth, by = "mouse"); setNames(d$log2FC_mouse, d$human) }
w_das <- { d <- mdasa %>% inner_join(orth, by = "mouse") %>% filter(!is.na(log2FC_mouse));
           setNames(sign(d$log2FC_mouse), d$human) }
w_sf  <- { d <- msf %>% inner_join(orth, by = "mouse") %>% distinct(human);
           setNames(rep(1, nrow(d)), d$human) }
set.seed(99)
w_rnd <- setNames(rep(1, 300), sample(orth$human, 300))
SETS  <- list("DEG signature" = w_deg, "DAS signature" = w_das,
              "splicing-factor/RBP module" = w_sf, "Random gene set (n=300)" = w_rnd)

score <- function(w, E) {
  g <- intersect(names(w), rownames(E)); if (length(g) < 5) return(NULL)
  ww <- w[g]
  list(genes = g, score = as.numeric(crossprod(ww, E[g, , drop = FALSE]) / sum(abs(ww))))
}

report_auc <- function(E, y, label) {
  out <- list()
  for (nm in names(SETS)) {
    ss <- score(SETS[[nm]], E); if (is.null(ss)) next
    r  <- pROC::roc(y, ss$score, quiet = TRUE)
    ci <- as.numeric(pROC::ci.auc(r))
    out[[nm]] <- list(roc = r, n_genes = length(ss$genes),
                      auc = as.numeric(pROC::auc(r)), lo = ci[1], hi = ci[3],
                      score = ss$score)
    say(sprintf("  %-26s %-11s AUC = %.3f (95%% CI %.3f-%.3f)   genes = %d",
                nm, label, out[[nm]]$auc, ci[1], ci[3], length(ss$genes)))
  }
  out
}
cmp <- function(a, b, lab) {
  if (is.null(a) || is.null(b)) { say("  DeLong ", lab, ": not available"); return(NULL) }
  d <- pROC::roc.test(a$roc, b$roc, method = "delong")
  say(sprintf("  DeLong %s: Z = %.3f, df = %.2f, P = %s",
              lab, as.numeric(unname(d$statistic)), as.numeric(unname(d$parameter)),
              format.pval(d$p.value, digits = 3, eps = 1e-4)))
  invisible(d)
}

## =========================== 1. GSE152004 ============================
say("=== GSE152004 (nasal airway epithelium, RNA-seq) ===")
h1 <- load_gse152004()
cpm1 <- log2(sweep(h1$counts, 2, colSums(h1$counts), "/") * 1e6 + 1)
y1 <- as.integer(h1$meta$group == "asthma")
say("  asthma = ", sum(y1 == 1), " ; control = ", sum(y1 == 0))
A1 <- report_auc(cpm1, y1, "GSE152004")
D1 <- cmp(A1[["DEG signature"]], A1[["DAS signature"]], "DEG vs DAS")
cmp(A1[["DEG signature"]], A1[["splicing-factor/RBP module"]], "DEG vs SF/RBP")

dlt <- A1[["DAS signature"]]$auc - A1[["DEG signature"]]$auc
set.seed(2026); B <- 2000
s_deg <- A1[["DEG signature"]]$score; s_das <- A1[["DAS signature"]]$score
i1 <- which(y1 == 1); i0 <- which(y1 == 0); dA <- numeric(B)
for (b in seq_len(B)) {
  ii <- c(sample(i1, length(i1), TRUE), sample(i0, length(i0), TRUE)); yy <- y1[ii]
  dA[b] <- as.numeric(pROC::auc(pROC::roc(yy, s_das[ii], quiet = TRUE))) -
           as.numeric(pROC::auc(pROC::roc(yy, s_deg[ii], quiet = TRUE)))
}
ci_d <- quantile(dA, c(.025, .975), na.rm = TRUE)
say(sprintf("  delta AUC (DAS - DEG) = %+.3f ; bootstrap 95%% CI %+.3f to %+.3f (B = %d)",
            dlt, ci_d[1], ci_d[2], B))

## ---- type-2 stratification within GSE152004 -------------------------
say("\n=== type-2 stratification (GSE152004) ===")
T2G <- c("CLCA1", "POSTN", "SERPINB2", "CCL26")
t2g <- intersect(T2G, rownames(cpm1))
say("  type-2 signature genes present: ", paste(t2g, collapse = ", "))
t2 <- list()
if (length(t2g) >= 3) {
  t2s <- rowMeans(scale(t(cpm1[t2g, , drop = FALSE])))
  med <- median(t2s[y1 == 1])
  hi  <- y1 == 1 & t2s >  med
  lo  <- y1 == 1 & t2s <= med
  say(sprintf("  type-2-high asthma = %d ; type-2-low asthma = %d ; controls = %d",
              sum(hi), sum(lo), sum(y1 == 0)))
  for (st in c("high", "low")) {
    sel <- c(which(if (st == "high") hi else lo), which(y1 == 0))
    yy <- ifelse(y1[sel] == 1, 1L, 0L)
    res <- list()
    for (nm in c("DEG signature", "DAS signature")) {
      ss <- score(SETS[[nm]], cpm1[, sel, drop = FALSE])
      r  <- pROC::roc(yy, ss$score, quiet = TRUE); ci <- as.numeric(pROC::ci.auc(r))
      res[[nm]] <- list(auc = as.numeric(pROC::auc(r)), lo = ci[1], hi = ci[3], roc = r)
      say(sprintf("  type-2-%-4s %-14s AUC = %.3f (95%% CI %.3f-%.3f)",
                  st, nm, res[[nm]]$auc, ci[1], ci[3]))
    }
    if (!is.null(res[["DEG signature"]]) && !is.null(res[["DAS signature"]]))
      cmp(res[["DEG signature"]], res[["DAS signature"]], paste0("DEG vs DAS in type-2-", st))
    t2[[st]] <- res
  }
}

## ============================ 2. GSE158752 ===========================
say("\n=== GSE158752 (bronchial epithelium, RNA-seq) ===")
try_catch({
  h2 <- load_gse158752()
  cpm2 <- log2(sweep(h2$counts, 2, colSums(h2$counts), "/") * 1e6 + 1)
  y2 <- as.integer(h2$group != "Healthy")
  say("  asthma = ", sum(y2 == 1), " ; healthy = ", sum(y2 == 0))
  A2 <<- report_auc(cpm2, y2, "GSE158752")
  cmp(A2[["DEG signature"]], A2[["DAS signature"]], "DEG vs DAS (GSE158752)")
}, "GSE158752")

## ============================ 3. microarray ==========================
say("\n=== microarray cohorts (probe-level matrices; reported only if mappable) ===")
A67 <- A43 <- NULL
try_catch({
  g67 <- load_micro("GSE67472", "asthma", "healthy")
  say("  GSE67472 rownames example: ", paste(head(rownames(g67$expr), 3), collapse = ", "))
  y67 <- as.integer(g67$group == "case")
  A67 <<- report_auc(g67$expr, y67, "GSE67472")
}, "GSE67472")
try_catch({
  g43 <- load_micro("GSE43696", "asthma", "control")
  say("  GSE43696 rownames example: ", paste(head(rownames(g43$expr), 3), collapse = ", "))
  y43 <- as.integer(g43$group == "case")
  A43 <<- report_auc(g43$expr, y43, "GSE43696")
}, "GSE43696")

## ============================ save ===================================
grab <- function(A) if (is.null(A)) NULL else lapply(A, function(x)
  c(auc = x$auc, lo = x$lo, hi = x$hi, n_genes = x$n_genes))
res <- list(delong_152004 = if (is.null(D1)) NULL else
              c(z = as.numeric(unname(D1$statistic)), p = D1$p.value),
            d_auc = dlt, d_lo = as.numeric(ci_d[1]), d_hi = as.numeric(ci_d[2]),
            auc_152004 = grab(A1), auc_158752 = grab(A2),
            auc_67472 = grab(A67), auc_43696 = grab(A43),
            t2 = lapply(t2, function(x) lapply(x, function(y) y[c("auc", "lo", "hi")])))
saveRDS(res, file.path(HDIR, "aair_extra.rds"))
say("\nDONE")
