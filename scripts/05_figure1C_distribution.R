#!/usr/bin/env Rscript
# =============================================================================
# Figure 1C — Proteome-wide distribution of one amino acid's Z-scores
#             (default: glutamine), highlighting proteins assigned to that
#             amino acid in Figure 1A and labelling selected genes.
#
# Input : human_aa_zscores.tsv
#         (from uniprot_aa_percent_zscore_10142026_update.py)
# Output: Figure1C_<AA>_zscore_distribution.pdf / .png
# Run   : Rscript Figure1C_zscore_distribution.R   (or source() in RStudio)
# Needs : ggplot2, ggrepel, patchwork
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

# ---- 1. Settings ------------------------------------------------------------
# Optional command-line arguments: <zscore_tsv> <out_dir> <x-axis label prefix>
args       <- commandArgs(trailingOnly = TRUE)
base_dir   <- "."   # run all scripts from the repository root
if (!file.exists(file.path(base_dir, "data/human_aa_percentages.tsv")))
  stop("Run this script from the repository root (the folder that contains data/).")
zscore_tsv <- if (length(args) >= 1) args[1] else
  file.path(base_dir, "data/human_aa_zscores_length_adjusted.tsv")
out_dir    <- if (length(args) >= 2) args[2] else
  file.path(base_dir, "results/Fig1C")
z_desc     <- if (length(args) >= 3) args[3] else "length-adjusted"

aa_focus    <- "Q"
aa_label    <- "Glutamine"
# Gene symbols to label; names are the text shown on the plot
label_genes <- c(TBP = "TBP", EP400 = "EP400", CBP = "CREBBP")
binwidth    <- 0.25
x_max       <- 10         # x-axis limit; proteins beyond it are listed in a note
inset_from  <- 1.5        # left edge of the zoomed tail inset
hi_col      <- "#C8960C"  # proteins dominant for the focus amino acid
lo_col      <- "grey80"   # all other proteins

aa_cols <- c("A","C","D","E","F","G","H","I","K","L",
             "M","N","P","Q","R","S","T","V","W","Y")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 2. Data -----------------------------------------------------------------
zt <- read.delim(zscore_tsv, check.names = FALSE, stringsAsFactors = FALSE)
z  <- as.matrix(zt[, aa_cols]); storage.mode(z) <- "numeric"

# Same assignment as Figure 1A: amino acid with the highest Z-score
dominant <- aa_cols[max.col(z, ties.method = "first")]

dom_lab <- paste0(aa_focus, "-dominant")
df <- data.frame(
  gene  = zt$Gene_Symbol,
  z     = z[, aa_focus],
  group = factor(ifelse(dominant == aa_focus, dom_lab, "Other proteins"),
                 levels = c("Other proteins", dom_lab)))
N     <- nrow(df)
n_dom <- sum(df$group == dom_lab)
off   <- df[df$z > x_max, ]
off   <- off[order(-off$z), ]
x_min <- floor(min(df$z))

# Genes to label (first entry if a symbol occurs more than once)
lab <- do.call(rbind, lapply(names(label_genes), function(nm) {
  i <- which(df$gene == label_genes[[nm]])[1]
  if (is.na(i)) { warning("Gene not found: ", label_genes[[nm]]); return(NULL) }
  data.frame(name = nm, z = df$z[i], group = df$group[i],
             pct_rank = mean(df$z >= df$z[i]) * 100)
}))
lab$label <- sprintf("%s\nZ = %.1f", lab$name, lab$z)
print(lab[, c("name", "z", "group", "pct_rank")], row.names = FALSE)

fill_scale <- scale_fill_manual(
  values = setNames(c(lo_col, hi_col), levels(df$group)),
  labels = c(sprintf("Other proteins (n = %s)", format(N - n_dom, big.mark = ",")),
             sprintf("%s (n = %s)", dom_lab, format(n_dom, big.mark = ","))),
  name = NULL, drop = FALSE)

# ---- 3. Main panel: whole distribution (density scale) -----------------------
# Bars are stacked, so the gold part of each bar is the share of that bin
# assigned to the focus amino acid; bar areas sum to 1 over all proteins, which
# makes the standard normal reference curve directly comparable.
dens <- function(count) count / (N * binwidth)

main <- ggplot(df[df$z <= x_max, ], aes(z)) +
  geom_histogram(aes(y = after_stat(dens(count)), fill = group),
                 binwidth = binwidth, boundary = 0, position = "stack",
                 colour = NA) +
  stat_function(fun = dnorm, n = 500, linewidth = 0.6, linetype = "22",
                colour = "grey25", xlim = c(x_min, 4)) +
  annotate("rect", xmin = inset_from, xmax = x_max, ymin = 0, ymax = 0.03,
           fill = NA, colour = "grey40", linewidth = 0.3, linetype = "dotted") +
  fill_scale +
  scale_x_continuous(breaks = seq(-2, x_max, 2), limits = c(x_min, x_max),
                     expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = sprintf("%s %s Z-score", aa_label, z_desc), y = "Density",
       title = sprintf("Proteome-wide distribution of %s (%s) Z-scores",
                       aa_label, aa_focus),
       subtitle = sprintf("n = %s proteins; dashed curve, standard normal reference",
                          format(N, big.mark = ","))) +
  theme_classic(base_size = 11) +
  theme(legend.position = c(0.01, 0.99), legend.justification = c(0, 1),
        legend.key.size = unit(3.5, "mm"), legend.text = element_text(size = 8.5),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(colour = "grey35", size = 9))

# ---- 4. Inset: right tail (protein counts) with labelled genes ---------------
tail_df <- df[df$z >= inset_from & df$z <= x_max, ]
y_max_tail <- max(table(cut(tail_df$z, seq(inset_from, x_max + binwidth, binwidth),
                            right = FALSE)))

inset <- ggplot(tail_df, aes(z)) +
  geom_histogram(aes(fill = group), binwidth = binwidth, boundary = 0,
                 position = "stack", colour = NA, show.legend = FALSE) +
  geom_text_repel(data = lab, aes(x = z, y = 0, label = label),
                  inherit.aes = FALSE, direction = "x", nudge_y = y_max_tail * 0.7,
                  size = 2.7, lineheight = 0.9, fontface = "bold",
                  segment.colour = "black", segment.size = 0.35,
                  min.segment.length = 0, box.padding = 0.3) +
  fill_scale +
  scale_x_continuous(limits = c(inset_from, x_max + binwidth),
                     breaks = seq(2, x_max, 2), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                     limits = c(0, y_max_tail * 1.05)) +
  labs(x = NULL, y = "Proteins") +
  theme_classic(base_size = 8) +
  theme(plot.background = element_rect(colour = "grey40", linewidth = 0.3,
                                       linetype = "dotted"))

if (nrow(off)) {
  inset <- inset + annotate(
    "text", x = x_max, y = y_max_tail * 1.02, hjust = 1, vjust = 1, size = 2.3,
    colour = "grey30",
    label = sprintf("%d protein%s with Z > %g not shown\n(max: %s, Z = %.1f)",
                    nrow(off), if (nrow(off) > 1) "s" else "", x_max,
                    off$gene[1], off$z[1]))
}

fig <- main + inset_element(inset, left = 0.47, bottom = 0.3, right = 0.995,
                            top = 0.97)

f <- file.path(out_dir, sprintf("Figure1C_%s_zscore_distribution", aa_focus))
ggsave(paste0(f, ".pdf"), fig, width = 6.5, height = 4.2, device = cairo_pdf)
ggsave(paste0(f, ".png"), fig, width = 6.5, height = 4.2, dpi = 300)
print(fig)
cat("Saved:", paste0(f, c(".pdf", ".png")), sep = "\n")
