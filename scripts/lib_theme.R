## =====================================================================
##  lib_theme.R  --  shared plotting theme / palettes / helpers
##  OVA chronic allergic airway inflammation : splicing manuscript
## =====================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(tibble)
  library(readr); library(stringr); library(forcats); library(purrr)
  library(patchwork); library(scales); library(RColorBrewer)
})

## ---- paths ----------------------------------------------------------
PROJ  <- "D:/学校工作/新建文件夹 (2)/OVA/manuscript_extension"
DATA  <- file.path(PROJ, "data")
FIGD  <- file.path(PROJ, "figures")
TABD  <- file.path(PROJ, "tables")
OUTD  <- file.path(PROJ, "out")
for (d in c(FIGD, TABD, OUTD)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

## ---- groups ---------------------------------------------------------
CTRL    <- c("Control-1", "Control-2", "Control-3", "Control-4")
OVAGRP  <- c("AA-1", "AA-2", "AA-3", "AA-4")
SAMPLES <- c(CTRL, OVAGRP)
GRP_LAB <- c("Control", "OVA")

## ---- palettes -------------------------------------------------------
PAL <- list(
  grp    = c(Control = "#2E6FA7", OVA = "#D1495B"),
  grp2   = c(Control = "#3C8DBC", OVA = "#E4572E"),
  cls    = c(SE = "#4C72B0", MXE = "#DD8452", RI = "#55A868",
             A3SS = "#C44E52", A5SS = "#8172B3"),
  seq    = c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8",
             "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026"),
  heat   = c("#2166AC", "#67A9CF", "#D1E5F0", "#F7F7F7",
             "#FDDBC7", "#EF8A62", "#B2182B"),
  div    = c("#0571B0", "#92C5DE", "#F7F7F7", "#F4A582", "#CA0020"),
  up     = "#C0392B", dn = "#2471A3", ns = "#B8BFC7",
  ok     = "#1B7837", bad = "#B2182B",
  accent = c("#2E6FA7", "#D1495B", "#3F8F5B", "#E0A32E", "#7B5EA7", "#4C9F9F")
)
CLS_LEVELS <- c("SE", "MXE", "RI", "A3SS", "A5SS")
ONT_LEVELS <- c("BP", "CC", "MF")

## ---- theme ----------------------------------------------------------
theme_pub <- function(base_size = 9.5, base_family = "sans",
                      legend_pos = "right", grid = TRUE) {
  th <- theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = if (grid) element_line(linewidth = 0.22, colour = "grey91")
                         else element_blank(),
      panel.border     = element_rect(linewidth = 0.55, colour = "grey25"),
      axis.line        = element_blank(),
      axis.ticks       = element_line(linewidth = 0.4, colour = "grey25"),
      axis.ticks.length = unit(1.6, "pt"),
      axis.text        = element_text(colour = "grey10", size = base_size - 0.7),
      axis.title       = element_text(colour = "black", size = base_size),
      plot.title       = element_text(face = "bold", size = base_size + 1.4,
                                      hjust = 0, margin = margin(b = 3)),
      plot.subtitle    = element_text(size = base_size - 1.2, colour = "grey30",
                                      hjust = 0, margin = margin(b = 3)),
      plot.tag         = element_text(face = "bold", size = base_size + 6),
      plot.tag.position = c(0, 1),
      legend.background = element_blank(),
      legend.key        = element_blank(),
      legend.key.size   = unit(9, "pt"),
      legend.title      = element_text(face = "bold", size = base_size - 0.3),
      legend.text       = element_text(size = base_size - 0.8),
      legend.position   = legend_pos,
      legend.margin     = margin(1, 1, 1, 1),
      strip.background  = element_rect(fill = "grey94", colour = "grey55",
                                       linewidth = 0.4),
      strip.text        = element_text(face = "bold", size = base_size - 0.5,
                                       colour = "grey10"),
      plot.margin       = margin(4, 6, 4, 4)
    )
  th
}
theme_set(theme_pub())

## panel tag helper used by every figure script.
## The bold letter is rendered as the leading run of the plot TITLE, which
## ggplot draws just above the panel -- i.e. outside the plotting area, close
## to the panel, and never covering any in-panel text.
esc_md <- function(x) gsub("([*_\\[\\]`~#])", "\\\\\\1", x)

## Panel-letter case.  MDPI uses lowercase (a, b, c ...); the R3 template uses
## uppercase (A, B, C ...).  Set PANEL_CASE=upper when rendering for R3.
PANEL_CASE <- Sys.getenv("PANEL_CASE", "lower")

ptag <- function(p, tag, size = 13, nspace = 3) {
  tag <- if (identical(PANEL_CASE, "upper")) toupper(tag) else tolower(tag)
  tt  <- p$labels$title
  gap <- strrep(" ", nspace)               # plain spaces lead the 2nd span:
                                           # gridtext keeps these (it trims NBSP)
  lab <- sprintf("<span style='font-size:%gpt'>**%s**</span>", size, tag)
  if (!is.null(tt) && nzchar(tt))
    lab <- sprintf(paste0("<span style='font-size:%gpt'>**%s**</span>",
                          "<span style='font-size:%gpt'>%s%s</span>"),
                   size, tag, size - 2.2, gap, esc_md(tt))
  p + labs(title = lab, tag = NULL) +
    theme(plot.title = ggtext::element_markdown(
      hjust = 0, face = "plain", lineheight = 1.1,
      margin = margin(b = 3, t = 1)))
}

## =====================================================================
##  MDPI / Biomedicines figure geometry
##  Full text width = 17.5 cm = 6.89 in; minimum recommended 600 dpi.
##  A figure laid out on a canvas of `w` inches with base font `fs` pt prints
##  at 6.89 in with an effective font of  fs * 6.89 / w  pt and an effective
##  resolution of  pixels / 6.89  dpi.
## =====================================================================
MDPI_W       <- 6.89
MDPI_MIN_DPI <- 600
mdpi_font <- function(fs, w) fs * MDPI_W / w        # effective printed pt
mdpi_dpi  <- function(pixels) pixels / MDPI_W       # effective printed dpi

## Render for MDPI: keep the canvas geometry the plot was designed for, but
## guarantee >= min_dpi when the image is placed at full text width.
save_fig_mdpi <- function(plot, name, width, height, base_font = 9.5,
                          formats = c("png", "pdf"), min_dpi = 700) {
  px_w  <- width * min_dpi
  dpi   <- ceiling(px_w / width)                    # == min_dpi
  save_fig(plot, name, width = width, height = height, dpi = dpi,
           formats = formats)
  cat(sprintf("      -> prints at %.1f pt effective font, %d dpi at %.2f in\n",
              mdpi_font(base_font, width), mdpi_dpi(width * dpi), MDPI_W))
  invisible(plot)
}

## a compact MDPI-sized theme: same layout, fonts sized for the canvas width
theme_mdpi <- function(canvas_w, base_pt = 6.6) {
  scale <- base_pt / 9.5                     # theme_pub() was designed at 9.5
  theme_pub(base_size = base_pt * scale * 9.5 / 6.6 * 1.0)
}

## ---------------------------------------------------------------------
##  MDPI output mode.  Set MDPI_MODE=1 to emit figures that satisfy
##  Biomedicines' requirements without changing the layout geometry:
##    * lowercase panel letters (handled by ptag())
##    * >=600 dpi when the image is placed at the 17.5 cm text width
##  The canvas is deliberately NOT rescaled.  Rescaling shrinks every
##  patchwork column, so long annotation labels and multi-column rows
##  overflow and collide.  Resolution alone is what the journal needs:
##  a 15.5-in canvas at 400 dpi is 6200 px, i.e. ~900 dpi at 17.5 cm.
## ---------------------------------------------------------------------
MDPI_MODE  <- identical(Sys.getenv("MDPI_MODE"), "1")
MDPI_DPI   <- suppressWarnings(as.numeric(Sys.getenv("MDPI_DPI", "400")))
MDPI_SUFFIX <- Sys.getenv("MDPI_SUFFIX", "")

if (MDPI_MODE) {
  cat(sprintf("  [MDPI mode] canvas unchanged, %d dpi; >=600 dpi at %.2f in\n",
              MDPI_DPI, MDPI_W))
}

## ---- inference helpers (shared by lib_gsea.R and lib_human.R) --------
boot_ci <- function(x, n = 5000, conf = 0.95, seed = 1) {
  set.seed(seed)
  x <- x[!is.na(x)]
  if (!length(x)) return(c(est = NA, lo = NA, hi = NA))
  b <- replicate(n, mean(sample(x, length(x), replace = TRUE)))
  q <- quantile(b, c((1 - conf) / 2, 1 - (1 - conf) / 2))
  c(est = mean(x), lo = unname(q[1]), hi = unname(q[2]))
}

overlap_perm <- function(setA, setB, universe, n = 20000, seed = 7) {
  set.seed(seed)
  u  <- unique(universe)
  a  <- intersect(setA, u); b <- intersect(setB, u)
  obs <- length(intersect(a, b))
  sim <- replicate(n, length(intersect(sample(u, length(a)), sample(u, length(b)))))
  tibble::tibble(observed = obs, expected = mean(sim),
                 p_greater = (sum(sim >= obs) + 1) / (n + 1),
                 p_less    = (sum(sim <= obs) + 1) / (n + 1),
                 p_two     = 2 * min((sum(sim >= obs) + 1) / (n + 1),
                                     (sum(sim <= obs) + 1) / (n + 1)),
                 n_universe = length(u), n_A = length(a), n_B = length(b))
}

## ---- small helpers --------------------------------------------------
fmt_p <- function(p, digits = 3) {
  ifelse(is.na(p), "NA",
  ifelse(p < 2.2e-16, "P < 2.2e-16",
  ifelse(p < 1e-4, sprintf("P = %.2e", p),
  ifelse(p < 0.001, sprintf("P = %.1e", p), sprintf("P = %.*f", digits, p)))))
}
fmt_q <- function(q) ifelse(q < 1e-4, sprintf("%.2e", q), sprintf("%.4f", q))
sci   <- function(x) formatC(x, format = "e", digits = 2)

## save a composed figure at high resolution, several formats
save_fig <- function(plot, name, width, height, dpi = 320, formats = c("png", "pdf")) {
  if (MDPI_MODE) {
    dpi  <- MDPI_DPI                     # geometry unchanged; only resolution
    name <- paste0(name, MDPI_SUFFIX)
  }
  for (f in formats) {
    fn <- file.path(FIGD, paste0(name, ".", f))
    if (f == "png") {
      ragg::agg_png(fn, width = width, height = height, units = "in",
                    res = dpi, background = "white")
    } else {
      grDevices::cairo_pdf(fn, width = width, height = height, bg = "white",
                           onefile = TRUE)
    }
    print(plot)
    grDevices::dev.off()
    cat("   saved", fn, "\n")
  }
  invisible(plot)
}

## a compact "stat box" for annotation
stat_box <- function(label, x = 0.98, y = 0.98, hjust = 1, vjust = 1,
                     size = 2.7, fill = "white", alpha = 0.82) {
  annotate("label", x = x, y = y, label = label, hjust = hjust, vjust = vjust,
           size = size, fill = fill, alpha = alpha, colour = "grey15",
           label.size = 0.15, label.padding = unit(0.16, "lines"),
           lineheight = 1.05)
}

## draw an exon-model schematic for one rMATS event (used in AS figures)
draw_event <- function(cls, x0 = 0, x1 = 1, y = 0, inc = 0.5, col = "grey20",
                       ymin = -0.32, ymax = 0.32) {
  ## returns a data.frame of segments for a schematic
  seg <- function(xs, xe, yy, h) data.frame(x = xs, xe = xe, y = yy, h = h)
  if (cls == "SE") {
    rbind(seg(x0, x0 + 0.14, y, 0.16), seg(x0 + 0.34, x0 + 0.66, y, 0.16),
          seg(x0 + 0.86, x1, y, 0.16), seg(x0 + 0.14, x0 + 0.86, y, 0.02))
  } else NULL
}

## load the analysis tables produced by 01_parse_data.py
read_tsv2 <- function(f) readr::read_tsv(file.path(DATA, f), show_col_types = FALSE)
