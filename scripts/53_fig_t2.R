## =====================================================================
##  53_fig_t2.R  --  supplementary figure S11
##  ROC curves of the expression- and splice-derived signatures within
##  type-2-high and type-2-low asthma (GSE152004).
## =====================================================================
.libPaths(c("D:/Rlibs", .libPaths()))
setwd("D:/学校工作/新建文件夹 (2)/OVA/manuscript_extension/R")
source("lib_theme.R"); source("lib_data.R")
suppressPackageStartupMessages({ library(dplyr); library(pROC) })
source("lib_human.R")

HR   <- readRDS(file.path(HDIR, "human_results.rds")); orth <- HR$orth
S    <- readRDS(file.path(HDIR, "aair_scores.rds"))
cpm1 <- S$cpm1; y1 <- S$y1

mdeg  <- deg_str %>% select(mouse = gene, log2FC_mouse = log2FC)
mdasa <- tibble(mouse = DAS_GENES) %>%
  left_join(deg_all %>% select(mouse = gene, log2FC_all = log2FC), by = "mouse") %>%
  mutate(log2FC_mouse = log2FC_all)
w_deg <- { d <- mdeg  %>% inner_join(orth, by = "mouse"); setNames(d$log2FC_mouse, d$human) }
w_das <- { d <- mdasa %>% inner_join(orth, by = "mouse") %>% filter(!is.na(log2FC_mouse));
           setNames(sign(d$log2FC_mouse), d$human) }
score <- function(w, E) { g <- intersect(names(w), rownames(E))
  as.numeric(crossprod(w[g], E[g, , drop = FALSE]) / sum(abs(w[g]))) }

s_deg <- score(w_deg, cpm1); s_das <- score(w_das, cpm1)
T2G <- c("CLCA1", "POSTN", "SERPINB2", "CCL26")
t2s <- rowMeans(scale(t(cpm1[intersect(T2G, rownames(cpm1)), , drop = FALSE])))
med <- median(t2s[y1 == 1])
hi  <- y1 == 1 & t2s >  med; lo <- y1 == 1 & t2s <= med
ctl <- y1 == 0

png(file.path("..", "figures", "FigureS11_type2_ROC.png"),
    width = 7.2, height = 3.6, units = "in", res = 400)
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.6, 1.0), cex = 0.95)
draw <- function(sel, ttl, ncase) {
  yy <- as.integer(y1[sel] == 1)
  rd <- pROC::roc(yy, s_deg[sel], quiet = TRUE, direction = "<")
  rs <- pROC::roc(yy, s_das[sel], quiet = TRUE, direction = "<")
  plot(1 - rd$specificities, rd$sensitivities, type = "l", lwd = 2.2, col = "#4C72B0",
       xlim = c(0, 1), ylim = c(0, 1), xlab = "1 - specificity", ylab = "sensitivity",
       main = ttl, las = 1)
  lines(1 - rs$specificities, rs$sensitivities, lwd = 2.2, col = "#C44E52")
  abline(0, 1, lty = 3, col = "grey55")
  legend("bottomright", bty = "n", lwd = 2.2, col = c("#4C72B0", "#C44E52"),
         cex = 0.85,
         legend = c(sprintf("Expression-derived (AUC %.3f)", as.numeric(pROC::auc(rd))),
                    sprintf("Splice-derived (AUC %.3f)", as.numeric(pROC::auc(rs)))))
  mtext(sprintf("cases = %d, controls = %d", ncase, sum(y1[sel] == 0)), side = 3, line = 0.1, cex = 0.75)
}
draw(c(which(hi), which(ctl)), "Type-2-high asthma", sum(hi))
draw(c(which(lo), which(ctl)), "Type-2-low asthma",  sum(lo))
dev.off()
cat("wrote FigureS11_type2_ROC.png\n")
