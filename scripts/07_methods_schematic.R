#!/usr/bin/env Rscript
# =============================================================================
# Methods schematic — workflow cartoon for Figure 1 (6 steps, 2 x 3 panels)
# Output: results/Methods_schematic/Methods_schematic.pdf/.png
# The PDF is vector and can be edited in Illustrator.
# =============================================================================

suppressPackageStartupMessages({
  library(grid)
  library(RColorBrewer)
})

base_dir <- "."   # run all scripts from the repository root
if (!file.exists(file.path(base_dir, "data/human_aa_percentages.tsv")))
  stop("Run this script from the repository root (the folder that contains data/).")
out_dir  <- file.path(base_dir, "results/Methods_schematic")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Palette -----------------------------------------------------------------
ink     <- "#222222"
muted   <- "#6B6B6B"
panelbg <- "#F6F6F4"
border  <- "#C9C9C4"
accent  <- "#C8960C"                     # glutamine / "dominant" highlight
viridis <- c("#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725")
ylgnbu  <- brewer.pal(9, "YlGnBu")
aa      <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]

set.seed(7)

# ---- Helpers -----------------------------------------------------------------
txt <- function(label, x, y, size = 7, col = ink, face = "plain", just = "centre",
                family = "sans", ...)
  grid.text(label, x, y, just = just,
            gp = gpar(fontsize = size, col = col, fontface = face, fontfamily = family), ...)

arrow_line <- function(x0, y0, x1, y1, col = muted, lwd = 1.2)
  grid.segments(x0, y0, x1, y1, gp = gpar(col = col, lwd = lwd, fill = col),
                arrow = arrow(length = unit(1.6, "mm"), type = "closed"))

mini_heat <- function(m, cols, x, y, w, h, na_col = "#EDEDED") {
  nr <- nrow(m); nc <- ncol(m)
  ramp <- colorRampPalette(cols)(100)
  rng  <- range(m, na.rm = TRUE)
  idx  <- round(1 + 99 * (m - rng[1]) / diff(rng))
  fill <- ifelse(is.na(m), na_col, ramp[pmax(1, pmin(100, idx))])
  xs <- x + (rep(seq_len(nc), each = nr) - 1) * w / nc
  ys <- y + h - rep(seq_len(nr), nc) * h / nr
  grid.rect(xs, ys, w / nc, h / nr, just = c("left", "bottom"),
            gp = gpar(fill = as.vector(fill), col = NA))
  grid.rect(x, y, w, h, just = c("left", "bottom"), gp = gpar(fill = NA, col = border, lwd = 0.4))
}

panel <- function(i, title, x, y, w, h, body) {
  vp <- viewport(x, y, w, h, just = c("left", "bottom"))
  pushViewport(vp)
  grid.roundrect(gp = gpar(fill = panelbg, col = border, lwd = 0.6), r = unit(2, "mm"))
  grid.circle(unit(4.2, "mm"), unit(1, "npc") - unit(4.2, "mm"), unit(2.6, "mm"),
              gp = gpar(fill = ink, col = NA))
  txt(i, unit(4.2, "mm"), unit(1, "npc") - unit(4.2, "mm"), size = 7.5, col = "white", face = "bold")
  txt(title, unit(8.5, "mm"), unit(1, "npc") - unit(4.2, "mm"), size = 8, face = "bold",
      just = "left")
  pushViewport(viewport(y = 0, height = unit(1, "npc") - unit(8.5, "mm"), just = "bottom"))
  body()
  popViewport(2)
}

# ---- Panel bodies --------------------------------------------------------------
body1 <- function() {
  txt("UniProt UP000005640 (FASTA)", 0.5, 0.92, size = 6, col = muted)
  seqs <- c(">sp|P20226|TBP_HUMAN", ">sp|Q92793|CBP_HUMAN", ">sp|Q96L91|EP400_HUMAN", ">...")
  for (k in seq_along(seqs)) {
    yy <- 0.81 - (k - 1) * 0.1
    txt(seqs[k], 0.08, yy, size = 5.6, col = ink, family = "mono", just = "left")
    if (k < 4) grid.segments(0.62, yy, 0.62 + c(0.22, 0.3, 0.26)[k], yy,
                             gp = gpar(col = viridis[2], lwd = 1.8, lineend = "round"))
  }
  # filter funnel
  fx <- 0.3; fy <- 0.33
  grid.polygon(c(fx - 0.14, fx + 0.14, fx + 0.035, fx + 0.035, fx - 0.035, fx - 0.035),
               c(fy + 0.09, fy + 0.09, fy - 0.03, fy - 0.11, fy - 0.11, fy - 0.03),
               gp = gpar(fill = "white", col = muted, lwd = 0.8))
  txt("remove empty sequences and\nnon-standard residues (U, X, …)", 0.47, fy, size = 5.8,
      col = muted, just = "left")
  txt("20,592 proteins", 0.5, 0.08, size = 8, face = "bold")
}

body2 <- function() {
  # protein as beads, Q highlighted
  s <- strsplit("MDQNQSLPQYQQGLASPQQ", "")[[1]]
  n <- length(s); xs <- seq(0.1, 0.9, length.out = n)
  ys <- 0.83 + 0.035 * sin(seq(0, 3 * pi, length.out = n))
  grid.lines(xs, ys, gp = gpar(col = border, lwd = 1))
  grid.circle(xs, ys, unit(1.55, "mm"),
              gp = gpar(fill = ifelse(s == "Q", accent, "white"), col = muted, lwd = 0.5))
  txt(s, xs, ys, size = 4.6, col = ifelse(s == "Q", "white", ink), face = "bold")
  txt(expression(paste("%", Q, " = ", frac("number of Q", "protein length"), " ", "\u00D7 100")),
      0.5, 0.6, size = 7)
  # bar chart of 20 AA percentages
  pct <- c(7, 2.3, 4.7, 7.1, 3.6, 6.6, 2.6, 4.3, 5.7, 10, 2.1, 3.6, 6.3, 4.7, 5.6, 8.3, 5.3, 6, 1.2, 2.7)
  pct[aa == "Q"] <- 11
  bx <- seq(0.1, 0.9, length.out = 20); bw <- 0.8 / 20 * 0.75
  grid.rect(bx, 0.2, bw, 0.2 * pct / 11, just = c("centre", "bottom"),
            gp = gpar(fill = ifelse(aa == "Q", accent, "grey65"), col = NA))
  grid.lines(c(0.08, 0.92), c(0.2, 0.2), gp = gpar(col = muted, lwd = 0.5))
  txt(aa, bx, 0.15, size = 4.6, col = ifelse(aa == "Q", accent, muted), face = "bold")
  txt("20 percentages per protein", 0.5, 0.06, size = 6, col = muted)
}

body3 <- function() {
  # proteins sorted into length bins (short -> long)
  nb <- 5
  bx <- seq(0.14, 0.86, length.out = nb)
  labs <- c("bin 1\nshort", "bin 2", "...", "bin 9", "bin 10\nlong")
  for (b in seq_len(nb)) {
    if (b == 3) { txt("\u2026", bx[b], 0.72, size = 10, col = muted); next }
    len <- c(0.04, 0.07, NA, 0.11, 0.14)[b]
    for (j in 1:4)
      grid.segments(bx[b] - len / 2, 0.86 - j * 0.035, bx[b] + len / 2, 0.86 - j * 0.035,
                    gp = gpar(col = viridis[2], lwd = 1.6, lineend = "round"))
    # bell curve: wider in short bins (noisier composition)
    sd <- c(0.05, 0.04, NA, 0.03, 0.025)[b]
    xx <- seq(-0.07, 0.07, length.out = 60)
    yy <- dnorm(xx, 0, sd); yy <- 0.13 * yy / max(dnorm(0, 0, 0.025))
    grid.polygon(c(bx[b] + xx, rev(bx[b] + xx)), c(0.47 + yy, rep(0.47, 60)),
                 gp = gpar(fill = adjustcolor(viridis[3], 0.35), col = viridis[3], lwd = 0.6))
    txt(labs[b], bx[b], 0.39, size = 5.3, col = muted, vjust = 1)
  }
  txt("10 equal-sized protein-length bins", 0.5, 0.95, size = 6, col = muted)
  txt(expression(Z == frac(x - mu[bin], sigma[bin])), 0.5, 0.2, size = 8)
  txt("each protein vs. proteins of similar length", 0.5, 0.045, size = 5.6, col = muted)
}

body4 <- function() {
  # Z-scores of one protein, max highlighted
  z <- c(-0.4, -0.9, 0.3, 0.8, -0.6, 0.2, 0.5, -1.0, 0.1, -0.2, -0.5, 0.9, 1.4, 3.6,
         0.2, 0.7, -0.3, -0.8, -0.7, -0.4)
  bx <- seq(0.06, 0.5, length.out = 20); bw <- 0.44 / 20 * 0.75; y0 <- 0.5; sc <- 0.09
  grid.rect(bx, y0 + pmin(z, 0) * sc, bw, abs(z) * sc, just = c("centre", "bottom"),
            gp = gpar(fill = ifelse(aa == "Q", accent, "grey65"), col = NA))
  txt("Z-scores of one protein", 0.28, 0.92, size = 5.6, col = muted)
  grid.lines(c(0.04, 0.52), c(y0, y0), gp = gpar(col = muted, lwd = 0.5))
  txt("Q", bx[aa == "Q"], y0 + 3.6 * sc + 0.04, size = 6, col = accent, face = "bold")
  txt("highest Z\n= dominant amino acid", 0.28, 0.12, size = 5.8, col = muted)
  arrow_line(0.54, 0.5, 0.62, 0.5)
  # mini Fig 1A: block-diagonal heatmap
  m <- matrix(rnorm(80 * 20, 0, 0.9), 80, 20)
  for (j in 1:20) m[(j - 1) * 4 + 1:4, j] <- seq(4, 2.4, length.out = 4)
  mini_heat(m, viridis, 0.64, 0.2, 0.32, 0.66)
  txt("20 groups (Fig. 1A)", 0.8, 0.1, size = 6, col = muted)
}

body5 <- function() {
  # 2 x 2 table
  x0 <- 0.03; y0 <- 0.46; cw <- 0.19; ch <- 0.13
  hdr <- c("in GO term", "not in term"); rows <- c("in group", "not in group")
  cells <- matrix(c("k", "n \u2212 k", "K \u2212 k", "N\u2212n\u2212K+k"), 2, byrow = TRUE)
  for (r in 1:2) for (c in 1:2) {
    grid.rect(x0 + 0.17 + (c - 1) * cw, y0 + (2 - r) * ch, cw, ch, just = c("left", "bottom"),
              gp = gpar(fill = if (r == 1 && c == 1) adjustcolor(ylgnbu[6], 0.35) else "white",
                        col = border, lwd = 0.5))
    txt(cells[r, c], x0 + 0.17 + (c - 0.5) * cw, y0 + (2.5 - r) * ch, size = 5.3)
  }
  txt(hdr, x0 + 0.17 + (1:2 - 0.5) * cw, y0 + 2 * ch + 0.04, size = 5.2, col = muted)
  txt(rows, x0 + 0.16, y0 + (2:1 - 0.5) * ch, size = 5.2, col = muted, just = "right")
  txt("one-sided hypergeometric test\nodds ratio (+0.5 per cell)\nBenjamini\u2013Hochberg FDR < 0.05",
      0.3, 0.25, size = 5.6, col = ink)
  txt("background: 17,632 GO-annotated proteins", 0.5, 0.04, size = 5.6, col = muted)
  arrow_line(0.585, 0.6, 0.645, 0.6)
  # mini Fig 1B
  m <- matrix(NA_real_, 60, 20)
  for (j in 1:20) m[(j - 1) * 3 + 1:3, j] <- runif(3, 1.5, 4.5)
  off <- cbind(sample(60, 30, TRUE), sample(20, 30, TRUE)); m[off] <- runif(30, 0.5, 2.5)
  mini_heat(m, ylgnbu[2:9], 0.665, 0.2, 0.3, 0.7)
  txt("Fig. 1B", 0.815, 0.13, size = 5.8, col = muted)
}

body6 <- function() {
  # merge related terms into a theme
  terms <- c("translation", "cytoplasmic translation", "ribosomal small subunit biogenesis")
  ty <- c(0.86, 0.76, 0.66)
  for (k in 1:3) {
    grid.roundrect(0.04, ty[k], 0.51, 0.075, just = c("left", "centre"), r = unit(0.8, "mm"),
                   gp = gpar(fill = "white", col = border, lwd = 0.5))
    txt(terms[k], 0.06, ty[k], size = 4.9, just = "left")
    grid.segments(0.55, ty[k], 0.6, 0.76, gp = gpar(col = muted, lwd = 0.6))
  }
  grid.roundrect(0.6, 0.76, 0.37, 0.14, just = c("left", "centre"), r = unit(1, "mm"),
                 gp = gpar(fill = adjustcolor(ylgnbu[6], 0.3), col = ylgnbu[7], lwd = 0.7))
  txt("theme:\nribosome/translation", 0.785, 0.76, size = 5.2, face = "bold")
  txt("FDR < 0.05, \u2265 10 proteins, OR \u2265 2;\nmerge terms sharing \u2265 50% of proteins",
      0.5, 0.52, size = 5.6, col = muted)
  # robustness check
  grid.lines(c(0.06, 0.94), c(0.41, 0.41), gp = gpar(col = border, lwd = 0.5, lty = 3))
  txt("Sensitivity analysis", 0.5, 0.34, size = 6.4, face = "bold")
  alts <- c("Z-score without length correction", "top 5% within each length bin")
  for (k in 1:2) {
    yy <- 0.24 - (k - 1) * 0.09
    grid.lines(c(0.12, 0.145, 0.185), c(yy, yy - 0.025, yy + 0.03),
               gp = gpar(col = "#2E8B57", lwd = 1.4, lineend = "round"))
    txt(alts[k], 0.22, yy, size = 5.6, just = "left")
  }
  txt("main functions reproduced (92\u201398%)", 0.5, 0.04, size = 5.6, col = muted)
}

# ---- Assemble ------------------------------------------------------------------
draw_schematic <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))
  W <- 0.315; H <- 0.43; gx <- (1 - 3 * W) / 4
  xs <- gx + (0:2) * (W + gx); ytop <- 0.52; ybot <- 0.04
  panel(1, "Human proteome",            xs[1], ytop, W, H, body1)
  panel(2, "Amino acid composition",    xs[2], ytop, W, H, body2)
  panel(3, "Length-adjusted Z-score",   xs[3], ytop, W, H, body3)
  panel(4, "Assign dominant amino acid", xs[1], ybot, W, H, body4)
  panel(5, "GO enrichment",             xs[2], ybot, W, H, body5)
  panel(6, "Themes & robustness",       xs[3], ybot, W, H, body6)
  # connectors
  for (k in 1:2) arrow_line(xs[k] + W + 0.003, ytop + H / 2, xs[k + 1] - 0.003, ytop + H / 2,
                            col = ink)
  for (k in 1:2) arrow_line(xs[k] + W + 0.003, ybot + H / 2, xs[k + 1] - 0.003, ybot + H / 2,
                            col = ink)
  # row wrap: panel 3 -> panel 4
  yy <- (ytop + ybot + H) / 2
  grid.lines(c(xs[3] + W / 2, xs[3] + W / 2, xs[1] + W / 2), c(ytop - 0.002, yy, yy),
             gp = gpar(col = ink, lwd = 1.2))
  arrow_line(xs[1] + W / 2, yy, xs[1] + W / 2, ybot + H + 0.002, col = ink)
  txt("Fig. 1C: Z-score distribution for one amino acid (e.g., glutamine)",
      0.5, yy + 0.013, size = 5.6, col = muted)
}

w <- 7.2; h <- 4.6
cairo_pdf(file.path(out_dir, "Methods_schematic.pdf"), width = w, height = h)
draw_schematic(); invisible(dev.off())
png(file.path(out_dir, "Methods_schematic.png"), width = w, height = h, units = "in",
    res = 300, type = "cairo")
draw_schematic(); invisible(dev.off())
cat("Saved to", out_dir, "\n")
