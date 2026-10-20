## =====================================================================
##  52_aair_delong.R
##  Fixed-orientation ROC comparison (DeLong) + bootstrap + type-2 strata.
##  Saves score vectors so it can be re-run without reloading the matrices.
## =====================================================================
.libPaths(c("D:/Rlibs", .libPaths()))
setwd("D:/学校工作/新建文件夹 (2)/OVA/manuscript_extension/R")
source("lib_theme.R"); source("lib_data.R")
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(pROC) })
source("lib_human.R")

REP <- file.path("..", "out", "aair_delong_report.txt"); unlink(REP)
say <- function(...) { m <- paste0(...); cat(m, "\n", sep = ""); cat(m, "\n", sep = "", file = REP, append = TRUE) }
f3  <- function(x) formatC(x, format = "f", digits = 3)
f4  <- function(x) formatC(x, format = "f", digits = 4)

HR   <- readRDS(file.path(HDIR, "human_results.rds")); orth <- HR$orth
mdeg  <- deg_str %>% select(mouse = gene, log2FC_mouse = log2FC)
mdasa <- tibble(mouse = DAS_GENES) %>%
  left_join(deg_all %>% select(mouse = gene, log2FC_all = log2FC), by = "mouse") %>%
  mutate(log2FC_mouse = log2FC_all)
msf  <- tibble(mouse = sf_state$gene[sf_state$in_expr])
w_deg <- { d <- mdeg  %>% inner_join(orth, by = "mouse"); setNames(d$log2FC_mouse, d$human) }
w_das <- { d <- mdasa %>% inner_join(orth, by = "mouse") %>% filter(!is.na(log2FC_mouse)); setNames(sign(d$log2FC_mouse), d$human) }
w_sf  <- { d <- msf %>% inner_join(orth, by = "mouse") %>% distinct(human); setNames(rep(1, nrow(d)), d$human) }
set.seed(99); w_rnd <- setNames(rep(1, 300), sample(orth$human, 300))
SETS <- list(DEG = w_deg, DAS = w_das, SF = w_sf, RND = w_rnd)

score <- function(w, E) { g <- intersect(names(w), rownames(E))
  list(genes = g, score = as.numeric(crossprod(w[g], E[g, , drop = FALSE]) / sum(abs(w[g])))) }

aucs <- function(s, y) {                      # both orientations
  ra <- pROC::roc(y, s, quiet = TRUE, direction = "auto")
  rf <- pROC::roc(y, s, quiet = TRUE, direction = "<")   # a priori: cases score higher
  c(auto = as.numeric(pROC::auc(ra)), fixed = as.numeric(pROC::auc(rf)))
}

## =========================== load / cache =============================
CF <- file.path(HDIR, "aair_scores.rds")
if (file.exists(CF)) { S <- readRDS(CF) } else {
  h1 <- load_gse152004()
  cpm1 <- log2(sweep(h1$counts, 2, colSums(h1$counts), "/") * 1e6 + 1)
  h2 <- load_gse158752()
  cpm2 <- log2(sweep(h2$counts, 2, colSums(h2$counts), "/") * 1e6 + 1)
  S <- list(cpm1 = cpm1, y1 = as.integer(h1$meta$group == "asthma"),
            cpm2 = cpm2, y2 = as.integer(h2$group != "Healthy"))
  saveRDS(S, CF)
}
cpm1 <- S$cpm1; y1 <- S$y1; cpm2 <- S$cpm2; y2 <- S$y2

sc1 <- lapply(SETS, score, E = cpm1)
sc2 <- lapply(SETS, score, E = cpm2)

say("=== GSE152004 : asthma ", sum(y1 == 1), " / control ", sum(y1 == 0), " ===")
for (nm in names(SETS)) {
  a <- aucs(sc1[[nm]]$score, y1)
  say("  ", sprintf("%-4s", nm), " genes=", length(sc1[[nm]]$genes),
      "  AUC(auto)=", f3(a["auto"]), "  AUC(a priori direction)=", f3(a["fixed"]))
}
dl <- function(s_a, s_b, y, lab) {
  ra <- pROC::roc(y, s_a, quiet = TRUE, direction = "<")
  rb <- pROC::roc(y, s_b, quiet = TRUE, direction = "<")
  d  <- pROC::roc.test(ra, rb, method = "delong")
  st <- as.numeric(unname(d$statistic)); pv <- as.numeric(d$p.value)
  say("  DeLong ", lab, ": Z = ", f3(st), ", P = ", format.pval(pv, digits = 3, eps = 1e-4))
  c(z = st, p = pv)
}
D_das <- dl(sc1$DAS$score, sc1$DEG$score, y1, "DAS vs DEG")
D_sf  <- dl(sc1$SF$score,  sc1$DEG$score, y1, "SF/RBP vs DEG")

## bootstrap delta-AUC with the a-priori orientation
set.seed(2026); B <- 2000
i1 <- which(y1 == 1); i0 <- which(y1 == 0); dA <- numeric(B)
for (b in seq_len(B)) {
  ii <- c(sample(i1, length(i1), TRUE), sample(i0, length(i0), TRUE)); yy <- y1[ii]
  dA[b] <- as.numeric(pROC::auc(pROC::roc(yy, sc1$DAS$score[ii], quiet = TRUE, direction = "<"))) -
           as.numeric(pROC::auc(pROC::roc(yy, sc1$DEG$score[ii], quiet = TRUE, direction = "<")))
}
ci <- quantile(dA, c(.025, .975), na.rm = TRUE)
dlt <- as.numeric(pROC::auc(pROC::roc(y1, sc1$DAS$score, quiet = TRUE, direction = "<"))) -
       as.numeric(pROC::auc(pROC::roc(y1, sc1$DEG$score, quiet = TRUE, direction = "<")))
say("  delta AUC (DAS - DEG) = ", f3(dlt), " ; bootstrap 95% CI ",
    f3(ci[1]), " to ", f3(ci[2]), " (B = ", B, ")")

## ------------------------------ type-2 strata ------------------------
say("\n=== type-2 strata (GSE152004) ===")
T2G <- c("CLCA1", "POSTN", "SERPINB2", "CCL26")
t2g <- intersect(T2G, rownames(cpm1))
t2s <- rowMeans(scale(t(cpm1[t2g, , drop = FALSE])))
med <- median(t2s[y1 == 1])
hi  <- y1 == 1 & t2s >  med; lo <- y1 == 1 & t2s <= med
say("  type-2 genes: ", paste(t2g, collapse = ","), " ; high = ", sum(hi), " low = ", sum(lo),
    " controls = ", sum(y1 == 0))
T2OUT <- list()
for (st in c("high", "low")) {
  sel <- c(which(if (st == "high") hi else lo), which(y1 == 0))
  yy <- as.integer(y1[sel] == 1)
  a_d <- aucs(sc1$DEG$score[sel], yy); a_s <- aucs(sc1$DAS$score[sel], yy)
  say("  type-2-", st, ": DEG AUC(auto) = ", f3(a_d["auto"]), " / DAS AUC(auto) = ", f3(a_s["auto"]),
      "   [a priori: DEG ", f3(a_d["fixed"]), " ; DAS ", f3(a_s["fixed"]), "]")
  dd <- dl(sc1$DAS$score[sel], sc1$DEG$score[sel], yy, paste0("DAS vs DEG within type-2-", st))
  T2OUT[[st]] <- list(deg_auto = a_d["auto"], das_auto = a_s["auto"],
                      deg_fixed = a_d["fixed"], das_fixed = a_s["fixed"], delong = dd, n = sum(yy))
}

## ------------------------------ GSE158752 ----------------------------
say("\n=== GSE158752 : asthma ", sum(y2 == 1), " / healthy ", sum(y2 == 0), " ===")
for (nm in names(SETS)) {
  a <- aucs(sc2[[nm]]$score, y2)
  ci_a <- as.numeric(pROC::ci.auc(pROC::roc(y2, sc2[[nm]]$score, quiet = TRUE, direction = "<")))
  say("  ", sprintf("%-4s", nm), " genes=", length(sc2[[nm]]$genes),
      "  AUC(a priori)=", f3(a["fixed"]), " (95% CI ", f3(ci_a[1]), "-", f3(ci_a[3]), ")" ,
      "  AUC(auto)=", f3(a["auto"]))
}
D_2 <- dl(sc2$DAS$score, sc2$DEG$score, y2, "DAS vs DEG (GSE158752)")

saveRDS(list(delong_das_152004 = D_das, delong_sf_152004 = D_sf,
             delta_auc = c(d = dlt, lo = as.numeric(ci[1]), hi = as.numeric(ci[2])),
             t2 = T2OUT, delong_158752 = D_2,
             scores = list(cpm1_dim = dim(cpm1), cpm2_dim = dim(cpm2))),
        file.path(HDIR, "aair_delong.rds"))
say("\nDONE")
