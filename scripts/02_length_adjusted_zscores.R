#!/usr/bin/env Rscript
# =============================================================================
# Step 2 — Length-adjusted amino acid composition Z-scores
#
# Short proteins have noisier compositions (one residue in a 50-aa protein is
# 2%), so proteome-wide Z-scores give them extreme values. Here each amino acid
# percentage is standardized against proteins of similar length instead:
#   proteins are split into 10 equal-sized protein-length deciles, and
#   Z = (x - mean_bin) / SD_bin   for each amino acid within each decile
#   (population SD, matching 01_aa_composition_zscores.py).
#
# Input : human_aa_percentages.tsv (from the Python script)
# Output: human_aa_zscores_length_adjusted.tsv  (same columns as the
#         unadjusted Z-score file, plus Length_decile)
#         length_adjustment_diagnostics.csv
# =============================================================================

base_dir <- "."   # run all scripts from the repository root
if (!file.exists(file.path(base_dir, "data/human_aa_percentages.tsv")))
  stop("Run this script from the repository root (the folder that contains data/).")
pct_tsv  <- file.path(base_dir, "data/human_aa_percentages.tsv")
raw_tsv  <- file.path(base_dir, "data/human_aa_zscores.tsv")
out_tsv  <- file.path(base_dir, "data/human_aa_zscores_length_adjusted.tsv")
diag_csv <- file.path(base_dir, "results/length_adjustment_diagnostics.csv")
n_bins   <- 10

aa_cols <- c("A","C","D","E","F","G","H","I","K","L",
             "M","N","P","Q","R","S","T","V","W","Y")
dir.create(dirname(diag_csv), showWarnings = FALSE, recursive = TRUE)

pt <- read.delim(pct_tsv, check.names = FALSE, stringsAsFactors = FALSE)
p  <- as.matrix(pt[, aa_cols]); storage.mode(p) <- "numeric"

# Equal-count length bins (ties broken by row order so bins are equal-sized)
bin <- ceiling(rank(pt$Length, ties.method = "first") / nrow(pt) * n_bins)

pop_sd <- function(x) sqrt(mean((x - mean(x))^2))
z <- p
for (b in seq_len(n_bins)) {
  i <- bin == b
  mu <- colMeans(p[i, , drop = FALSE])
  s  <- apply(p[i, , drop = FALSE], 2, pop_sd)
  z[i, ] <- sweep(sweep(p[i, , drop = FALSE], 2, mu), 2, ifelse(s == 0, 1, s), "/")
}

out <- data.frame(pt[, c("Protein_ID", "Gene_Symbol", "Length")],
                  round(z, 3), Length_decile = bin, check.names = FALSE)
write.table(out, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

# ---- Diagnostics: does the dominant Z-score still depend on length? ----------
zr <- as.matrix(read.delim(raw_tsv, check.names = FALSE)[, aa_cols])
dom_raw <- aa_cols[max.col(zr, ties.method = "first")]
dom_adj <- aa_cols[max.col(z,  ties.method = "first")]
diag <- data.frame(
  Length_decile = seq_len(n_bins),
  Length_range  = tapply(pt$Length, bin, function(x) paste(min(x), max(x), sep = "-")),
  Proteins      = as.vector(table(bin)),
  Median_dominant_Z_unadjusted = round(tapply(apply(zr, 1, max), bin, median), 2),
  Median_dominant_Z_adjusted   = round(tapply(apply(z, 1, max), bin, median), 2),
  Share_of_Zmax_gt4_unadjusted = round(as.vector(table(factor(bin[apply(zr, 1, max) > 4],
                                   levels = 1:n_bins))) / sum(apply(zr, 1, max) > 4), 3),
  Share_of_Zmax_gt4_adjusted   = round(as.vector(table(factor(bin[apply(z, 1, max) > 4],
                                   levels = 1:n_bins))) / sum(apply(z, 1, max) > 4), 3))
write.csv(diag, diag_csv, row.names = FALSE)
print(diag, row.names = FALSE)
cat(sprintf("\nProteins keeping the same dominant amino acid after adjustment: %d / %d (%.1f%%)\n",
            sum(dom_raw == dom_adj), nrow(pt), 100 * mean(dom_raw == dom_adj)))
cat("Written:", out_tsv, "\n")
