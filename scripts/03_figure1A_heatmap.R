#!/usr/bin/env Rscript
# =============================================================================
# Figure 1A — Proteome heatmap ordered by dominant amino acid bias (Z-scores)
#
# Optional command-line arguments: <zscore_tsv> <out_dir>
# Default input is the length-adjusted Z-score file from
# Step1_length_adjusted_zscores.R.
# =============================================================================

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# ---- 1. Input -----------------------------------------------------------------
args     <- commandArgs(trailingOnly = TRUE)
base_dir <- "."   # run all scripts from the repository root
if (!file.exists(file.path(base_dir, "data/human_aa_percentages.tsv")))
  stop("Run this script from the repository root (the folder that contains data/).")
tsv_file <- if (length(args) >= 1) args[1] else
  file.path(base_dir, "data/human_aa_zscores_length_adjusted.tsv")
out_dir  <- if (length(args) >= 2) args[2] else
  file.path(base_dir, "results/Fig1A")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
output_prefix <- file.path(out_dir, "Figure1A_zscore_heatmap")

# ---- 2. Read data ----
aa_cols <- c("A","C","D","E","F","G","H","I","K","L",
             "M","N","P","Q","R","S","T","V","W","Y")

dt <- read.delim(tsv_file, header = TRUE, check.names = FALSE)

missing_cols <- setdiff(aa_cols, names(dt))
if (length(missing_cols) > 0) {
  stop("Missing amino acid columns: ", paste(missing_cols, collapse = ", "))
}

z <- as.matrix(dt[, aa_cols])
storage.mode(z) <- "numeric"

# ---- 3. Order proteins by dominant amino acid ----
dominant_idx <- max.col(z, ties.method = "first")
dominant_z   <- z[cbind(seq_len(nrow(z)), dominant_idx)]

ord <- order(dominant_idx, -dominant_z)
z_sorted <- z[ord, ]

# ---- 4. Build heatmap ----
col_fun <- colorRamp2(c(-4, 0, 4), c("#440154", "#21908C", "#FDE725"))

ht <- Heatmap(
  z_sorted,
  name = "Z-score",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_rot = 0,
  column_names_gp = gpar(fontsize = 10),
  column_title = "Proteome heatmap ordered by dominant amino-acid bias (Z-scored)",
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_title = sprintf("Proteome (n = %s)", format(nrow(z), big.mark = ",")),
  use_raster = TRUE,
  raster_quality = 4      # avoids white seam lines in the rasterized body
)

# ---- 5. Save (draw() is required) ----
png(paste0(output_prefix, ".png"), width = 2100, height = 3000, res = 300)
draw(ht)
invisible(dev.off())

pdf(paste0(output_prefix, ".pdf"), width = 7, height = 10)
draw(ht)
invisible(dev.off())

cat("Saved:", paste0(output_prefix, ".png"), "and", paste0(output_prefix, ".pdf"), "\n")
