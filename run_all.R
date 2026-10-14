#!/usr/bin/env Rscript
# =============================================================================
# Reproduce Figure 1 (A, B, C), Table S2, the sensitivity analysis and the
# methods schematic from the amino acid percentage / Z-score tables in data/.
#
#   Rscript run_all.R          (run from the repository root)
#
# Step 1 (Python, optional) recreates data/ from the UniProt FASTA; see README.
# =============================================================================

if (!file.exists("data/human_aa_percentages.tsv"))
  stop("Run from the repository root (the folder that contains data/).")

rscript <- file.path(R.home("bin"), "Rscript")
z_adj   <- "data/human_aa_zscores_length_adjusted.tsv"
z_raw   <- "data/human_aa_zscores.tsv"

run <- function(script, ...) {
  cat("\n=====", script, paste(c(...), collapse = " "), "=====\n")
  status <- system2(rscript, shQuote(c(file.path("scripts", script), ...)))
  if (status != 0) stop(script, " failed")
}

run("02_length_adjusted_zscores.R")
run("03_figure1A_heatmap.R",       z_adj, "results/Fig1A")
run("04_figure1B_GO_enrichment.R", z_adj, "dominant", "results/Fig1B_GO")
run("05_figure1C_distribution.R",  z_adj, "results/Fig1C", "length-adjusted")

# Sensitivity analysis (alternative group definitions)
run("04_figure1B_GO_enrichment.R", z_raw, "dominant",
    "results/Sensitivity/GO_unadjusted_Z_dominant")
run("04_figure1B_GO_enrichment.R", z_adj, "top5_length_decile",
    "results/Sensitivity/GO_top5_length_decile")
run("06_sensitivity_compare.R")

run("07_methods_schematic.R")
cat("\nAll done. Outputs are in results/\n")
