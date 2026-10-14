#!/usr/bin/env Rscript
# =============================================================================
# Sensitivity analysis — are the GO enrichments robust to how amino
# acid-enriched groups are defined?
#
# Primary     : dominant amino acid from length-adjusted Z-scores
# Alternative : dominant amino acid from unadjusted (proteome-wide) Z-scores
# Alternative : top 5% by amino acid percentage within protein-length deciles
#
# Requires the three Figure1B_GO_enrichment.R runs made by run_all.R.
# Output (results/Sensitivity/):
#   FigureS_sensitivity.pdf/.png, TableS3_*.csv
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

base_dir <- "."   # run all scripts from the repository root
if (!file.exists(file.path(base_dir, "data/human_aa_percentages.tsv")))
  stop("Run this script from the repository root (the folder that contains data/).")
res_dir  <- file.path(base_dir, "results")
out_dir  <- file.path(res_dir, "Sensitivity")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

runs <- c(primary    = "Fig1B_GO",
          unadjusted = "Sensitivity/GO_unadjusted_Z_dominant",
          top5       = "Sensitivity/GO_top5_length_decile")
run_lab <- c(unadjusted = "Unadjusted Z, dominant AA",
             top5       = "Top 5% within length deciles")
top_n   <- 10          # "main" functions = top N terms per amino acid by FDR
fdr_cut <- 0.05

aa_cols <- c("A","C","D","E","F","G","H","I","K","L",
             "M","N","P","Q","R","S","T","V","W","Y")

read_run <- function(d) as.data.frame(readxl::read_excel(
  file.path(res_dir, d, "TableS2_Fig1B_source_data.xlsx"), sheet = "All tests"))
tests <- lapply(runs, read_run)
key   <- function(t) paste(t$amino_acid, t$GOID)

# ---- 1. Group membership overlap ---------------------------------------------
zdir <- file.path(base_dir, "data")
dom_of <- function(f) {
  z <- as.matrix(read.delim(file.path(zdir, f), check.names = FALSE)[, aa_cols])
  aa_cols[max.col(z, ties.method = "first")]
}
pt  <- read.delim(file.path(zdir, "human_aa_percentages.tsv"), check.names = FALSE)
ids <- pt$Protein_ID
dom_adj <- dom_of("human_aa_zscores_length_adjusted.tsv")
dom_raw <- dom_of("human_aa_zscores.tsv")
decile  <- cut(pt$Length, quantile(pt$Length, 0:10 / 10), include.lowest = TRUE,
               labels = FALSE)
top5 <- lapply(setNames(aa_cols, aa_cols), function(a)
  unlist(lapply(split(seq_along(ids), decile), function(i)
    ids[i][pt[i, a] >= quantile(pt[i, a], 0.95)]), use.names = FALSE))
jac <- function(x, y) length(intersect(x, y)) / length(union(x, y))

# ---- 2. Per-amino acid concordance -------------------------------------------
concord <- do.call(rbind, lapply(aa_cols, function(a) {
  P   <- tests$primary[tests$primary$amino_acid == a, ]
  sig <- P[P$fdr < fdr_cut, ]
  top <- head(sig[order(sig$fdr), ], top_n)
  do.call(rbind, lapply(names(run_lab), function(r) {
    A  <- tests[[r]][tests[[r]]$amino_acid == a, ]
    sA <- A$GOID[A$fdr < fdr_cut]
    sh <- merge(P, A, by = "GOID")
    alt_grp <- if (r == "top5") top5[[a]] else ids[dom_raw == a]
    data.frame(
      amino_acid = a, comparison = run_lab[[r]],
      group_jaccard = jac(ids[dom_adj == a], alt_grp),
      primary_sig_terms = nrow(sig), alt_sig_terms = length(sA),
      pct_primary_sig_recovered = 100 * mean(sig$GOID %in% sA),
      pct_top_terms_recovered = 100 * mean(top$GOID %in% sA),
      log2OR_spearman = if (nrow(sh) >= 5)
        cor(sh$log2_odds_ratio.x, sh$log2_odds_ratio.y, method = "spearman")
        else NA_real_)
  }))
}))
write.csv(concord, file.path(out_dir, "TableS3_concordance_by_amino_acid.csv"),
          row.names = FALSE)

overall <- do.call(rbind, lapply(names(run_lab), function(r) {
  c_r <- concord[concord$comparison == run_lab[[r]], ]
  sigP <- tests$primary[tests$primary$fdr < fdr_cut, ]
  sigA <- tests[[r]][tests[[r]]$fdr < fdr_cut, ]
  sh   <- merge(tests$primary, tests[[r]], by = c("amino_acid", "GOID"))
  data.frame(
    comparison = run_lab[[r]],
    median_group_jaccard = median(c_r$group_jaccard),
    primary_sig = nrow(sigP), alt_sig = nrow(sigA),
    pct_primary_sig_recovered = 100 * mean(key(sigP) %in% key(sigA)),
    pct_top_terms_recovered = 100 * weighted.mean(
      c_r$pct_top_terms_recovered / 100, pmin(c_r$primary_sig_terms, top_n)),
    log2OR_spearman_all_shared_tests = cor(sh$log2_odds_ratio.x,
                                           sh$log2_odds_ratio.y,
                                           method = "spearman"))
}))
write.csv(overall, file.path(out_dir, "TableS3_concordance_overall.csv"),
          row.names = FALSE)

# ---- 3. Top primary terms per amino acid, side by side -----------------------
top_tab <- do.call(rbind, lapply(aa_cols, function(a) {
  P <- tests$primary[tests$primary$amino_acid == a & tests$primary$fdr < fdr_cut, ]
  P <- head(P[order(P$fdr), ], top_n)
  if (!nrow(P)) return(NULL)
  out <- data.frame(amino_acid = a, rank = seq_len(nrow(P)), GOID = P$GOID,
                    term = P$term, primary_log2OR = round(P$log2_odds_ratio, 2),
                    primary_FDR = signif(P$fdr, 3))
  for (r in names(run_lab)) {
    A <- tests[[r]][tests[[r]]$amino_acid == a, ]
    m <- match(P$GOID, A$GOID)
    out[[paste0(r, "_log2OR")]] <- round(A$log2_odds_ratio[m], 2)
    out[[paste0(r, "_FDR")]]    <- signif(A$fdr[m], 3)
  }
  out
}))
write.csv(top_tab, file.path(out_dir, "TableS3_top_terms_comparison.csv"),
          row.names = FALSE)

# ---- 4. Supplementary figure -------------------------------------------------
cols <- c("Unadjusted Z, dominant AA" = "#4C78A8",
          "Top 5% within length deciles" = "#E39C37")
concord$amino_acid <- factor(concord$amino_acid, levels = aa_cols)

pa <- ggplot(concord, aes(amino_acid, pct_top_terms_recovered, fill = comparison)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_hline(yintercept = 100, colour = "grey60", linewidth = 0.3) +
  scale_fill_manual(values = cols, name = NULL) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Amino acid",
       y = sprintf("Top-%d primary terms\nalso significant (%%)", top_n),
       title = "a  Recovery of the main enriched functions") +
  theme_classic(base_size = 10) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"))

sc <- do.call(rbind, lapply(names(run_lab), function(r) {
  sh <- merge(tests$primary, tests[[r]], by = c("amino_acid", "GOID"))
  sh <- sh[sh$fdr.x < fdr_cut | sh$fdr.y < fdr_cut, ]
  data.frame(comparison = run_lab[[r]], x = sh$log2_odds_ratio.x,
             y = sh$log2_odds_ratio.y)
}))
rho <- do.call(rbind, lapply(split(sc, sc$comparison), function(d)
  data.frame(comparison = d$comparison[1],
             lab = sprintf("rho == %.2f", cor(d$x, d$y, method = "spearman")))))

pb <- ggplot(sc, aes(x, y, colour = comparison)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_point(size = 0.6, alpha = 0.5, show.legend = FALSE) +
  geom_text(data = rho, aes(x = -Inf, y = Inf, label = lab), parse = TRUE,
            hjust = -0.2, vjust = 1.5, colour = "black", size = 3.2) +
  facet_wrap(~ comparison) +
  scale_colour_manual(values = cols) +
  labs(x = "log2 odds ratio, primary (length-adjusted Z)",
       y = "log2 odds ratio, alternative",
       title = "b  Agreement of effect sizes (pairs significant in either analysis)") +
  theme_classic(base_size = 10) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))

fig <- pa / pb + plot_layout(heights = c(1, 1.1))
ggsave(file.path(out_dir, "FigureS_sensitivity.pdf"), fig, width = 7.5, height = 7,
       device = cairo_pdf)
ggsave(file.path(out_dir, "FigureS_sensitivity.png"), fig, width = 7.5, height = 7,
       dpi = 300)

print(overall, row.names = FALSE)
cat("Outputs written to", out_dir, "\n")
