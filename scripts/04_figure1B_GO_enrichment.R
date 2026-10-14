#!/usr/bin/env Rscript
# =============================================================================
# Figure 1B — GO Biological Process over-representation in amino acid-enriched
#             protein groups of the human proteome
#
# Input : <prefix>_aa_zscores.tsv and <prefix>_aa_percentages.tsv produced by
#         uniprot_aa_percent_zscore_10142026_update.py
#         (columns: Protein_ID, Gene_Symbol, Length, A ... Y)
# Output: Figure1B_GO_heatmap.pdf / .png
#         TableS2_Fig1B_source_data.xlsx (or .csv files if writexl is missing)
#         sessionInfo.txt  (package + annotation versions, for the Methods)
#
# Run   : Rscript Figure1B_GO_enrichment.R      (or source() it in RStudio)
# Needs : Bioconductor org.Hs.eg.db, GO.db, ComplexHeatmap; CRAN circlize,
#         RColorBrewer; optional writexl (install.packages("writexl"))
# =============================================================================

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(GO.db)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(grid)
})

# ---- 1. Settings ------------------------------------------------------------
# Optional command-line arguments: <zscore_tsv> <grouping> <out_dir>
args       <- commandArgs(trailingOnly = TRUE)
base_dir   <- "."   # run all scripts from the repository root
if (!file.exists(file.path(base_dir, "data/human_aa_percentages.tsv")))
  stop("Run this script from the repository root (the folder that contains data/).")
zscore_tsv <- if (length(args) >= 1) args[1] else
  file.path(base_dir, "data/human_aa_zscores_length_adjusted.tsv")
pct_tsv    <- file.path(base_dir, "data/human_aa_percentages.tsv")
out_dir    <- if (length(args) >= 3) args[3] else
  file.path(base_dir, "results/Fig1B_GO")

# How amino acid-enriched groups are defined:
#   "dominant"            each protein is assigned to the amino acid with its
#                         highest z-score (identical to the Figure 1A ordering);
#                         groups are mutually exclusive and differ in size.
#   "top5_length_decile"  for each amino acid, the top 5% of proteins by that
#                         amino acid's percentage within each protein-length
#                         decile (the definition used for the previous Table S2);
#                         groups are ~equal in size and may overlap.
grouping       <- if (length(args) >= 2) args[2] else "dominant"
top_fraction   <- 0.05      # used only by "top5_length_decile"

go_propagate   <- FALSE     # FALSE = direct GO annotations ("GO");
                            # TRUE  = include ancestor terms ("GOALL")
min_term_size  <- 20        # GO term size limits, counted in the background
max_term_size  <- 5000
min_overlap    <- 5         # min. group proteins annotated to the term
fdr_cutoff     <- 0.05
haldane        <- 0.5       # Haldane–Anscombe correction added to each cell
color_cap      <- 4.5       # upper limit of the log2 OR colour scale

# Representative themes per amino acid (for the text / summary table only;
# Figure 1B still shows every significant association)
rep_min_overlap <- 10       # >= this many group proteins in the term
rep_min_or      <- 2        # odds ratio >= this (log2 OR >= 1)
theme_overlap   <- 0.5      # a term joins a theme if >= this fraction of its
                            # group proteins (or of the lead term's, if smaller) are shared

aa_cols  <- c("A","C","D","E","F","G","H","I","K","L",
              "M","N","P","Q","R","S","T","V","W","Y")
aa_names <- c(A="Alanine", C="Cysteine", D="Aspartate", E="Glutamate",
              F="Phenylalanine", G="Glycine", H="Histidine", I="Isoleucine",
              K="Lysine", L="Leucine", M="Methionine", N="Asparagine",
              P="Proline", Q="Glutamine", R="Arginine", S="Serine",
              T="Threonine", V="Valine", W="Tryptophan", Y="Tyrosine")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 2. Read proteome table -------------------------------------------------
zt <- read.delim(zscore_tsv, check.names = FALSE, stringsAsFactors = FALSE)
pt <- read.delim(pct_tsv,    check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(identical(zt$Protein_ID, pt$Protein_ID), all(aa_cols %in% names(zt)))

# UniProt accession from "sp|P04637|P53_HUMAN"; plain accessions pass through
zt$Accession <- ifelse(grepl("\\|", zt$Protein_ID),
                       vapply(strsplit(zt$Protein_ID, "\\|"), `[`, "", 2),
                       zt$Protein_ID)
proteins <- zt$Accession
n_proteome <- length(proteins)
cat("Proteins in analysed proteome:", n_proteome, "\n")

# ---- 3. Define amino acid-enriched groups -----------------------------------
z <- as.matrix(zt[, aa_cols]); storage.mode(z) <- "numeric"
p <- as.matrix(pt[, aa_cols]); storage.mode(p) <- "numeric"

groups <- if (grouping == "dominant") {
  dom <- aa_cols[max.col(z, ties.method = "first")]
  setNames(lapply(aa_cols, function(a) proteins[dom == a]), aa_cols)
} else if (grouping == "top5_length_decile") {
  decile <- cut(zt$Length, quantile(zt$Length, 0:10 / 10), include.lowest = TRUE,
                labels = FALSE)
  setNames(lapply(aa_cols, function(a) {
    unlist(lapply(split(seq_len(n_proteome), decile), function(idx) {
      thr <- quantile(p[idx, a], 1 - top_fraction, type = 7)
      proteins[idx][p[idx, a] >= thr]
    }), use.names = FALSE)
  }), aa_cols)
} else stop("Unknown grouping: ", grouping)

# ---- 4. GO Biological Process annotations -----------------------------------
# Map UniProt accession -> Entrez; fall back to gene symbol when unmapped.
acc2eg <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db, keys = unique(proteins), keytype = "UNIPROT", columns = "ENTREZID"))
acc2eg <- acc2eg[!is.na(acc2eg$ENTREZID), ]
unmapped <- setdiff(proteins, acc2eg$UNIPROT)
sym <- setNames(zt$Gene_Symbol, proteins)[unmapped]
sym <- sym[!is.na(sym) & sym != "NA"]
if (length(sym)) {
  sym2eg <- suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db, keys = unique(sym), keytype = "SYMBOL", columns = "ENTREZID"))
  sym2eg <- sym2eg[!is.na(sym2eg$ENTREZID), ]
  fb <- data.frame(UNIPROT = names(sym), SYMBOL = unname(sym))
  fb <- merge(fb, sym2eg, by = "SYMBOL")[, c("UNIPROT", "ENTREZID")]
  acc2eg <- rbind(acc2eg, fb)
}

go_col  <- if (go_propagate) "GOALL" else "GO"
ont_col <- if (go_propagate) "ONTOLOGYALL" else "ONTOLOGY"
eg2go <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db, keys = unique(acc2eg$ENTREZID), keytype = "ENTREZID",
  columns = c(go_col, ont_col)))
eg2go <- eg2go[!is.na(eg2go[[go_col]]) & eg2go[[ont_col]] == "BP",
               c("ENTREZID", go_col)]
names(eg2go)[2] <- "GOID"

ann <- unique(merge(acc2eg, eg2go, by = "ENTREZID")[, c("UNIPROT", "GOID")])
go2prot <- split(ann$UNIPROT, ann$GOID)

# Background = analysed proteins with >= 1 GO BP annotation
universe <- intersect(proteins, unique(ann$UNIPROT))
N <- length(universe)
cat("Proteins with >= 1 GO BP annotation (background):", N, "\n")

term_size <- lengths(go2prot)
go2prot   <- go2prot[term_size >= min_term_size & term_size <= max_term_size]
cat("GO BP terms with", min_term_size, "-", max_term_size, "proteins:",
    length(go2prot), "\n")

# ---- 5. Hypergeometric test + Haldane–Anscombe odds ratio -------------------
res <- do.call(rbind, lapply(aa_cols, function(a) {
  grp <- intersect(groups[[a]], universe)
  n   <- length(grp)
  k   <- vapply(go2prot, function(g) sum(g %in% grp), integer(1))
  K   <- lengths(go2prot)
  keep <- k >= min_overlap
  if (!any(keep)) return(NULL)
  k <- k[keep]; K <- K[keep]
  # 2x2: [in group & term, in group & not term; not group & term, neither]
  a11 <- k + haldane;         a12 <- n - k + haldane
  a21 <- K - k + haldane;     a22 <- N - n - K + k + haldane
  or  <- (a11 * a22) / (a12 * a21)
  se  <- sqrt(1 / a11 + 1 / a12 + 1 / a21 + 1 / a22)   # Woolf SE of log OR
  data.frame(
    amino_acid      = a,
    amino_acid_name = aa_names[[a]],
    GOID            = names(k),
    group_size      = n,
    background      = N,
    term_size       = K,
    overlap         = k,
    expected        = n * K / N,
    odds_ratio      = or,
    or_ci95_low     = exp(log(or) - qnorm(0.975) * se),
    or_ci95_high    = exp(log(or) + qnorm(0.975) * se),
    log2_odds_ratio = log2(or),
    p_value         = phyper(k - 1, K, N - K, n, lower.tail = FALSE),
    stringsAsFactors = FALSE)
}))
res$fdr  <- p.adjust(res$p_value, method = "BH")   # global, all AA x GO tests
res$term <- unname(Term(res$GOID))
cat("Amino acid–GO tests performed:", nrow(res), "\n")

sig <- res[res$fdr < fdr_cutoff & res$log2_odds_ratio > 0, ]
sig <- sig[order(match(sig$amino_acid, aa_cols), sig$p_value), ]
cat("Significant associations (FDR <", fdr_cutoff, "):", nrow(sig),
    "involving", length(unique(sig$GOID)), "GO terms\n")

# Group proteins annotated to each significant term (UniProt accessions)
acc2sym <- setNames(zt$Gene_Symbol, proteins)
sig_members <- mapply(function(a, go) intersect(groups[[a]], go2prot[[go]]),
                      sig$amino_acid, sig$GOID, SIMPLIFY = FALSE)
sig$genes <- vapply(sig_members, function(x)
  paste(sort(unique(acc2sym[x])), collapse = ";"), "")

# ---- 5b. Representative themes per amino acid --------------------------------
# Eligible: significant, >= rep_min_overlap group proteins, OR >= rep_min_or.
# Ranked by FDR; the best remaining term starts a theme and absorbs every
# lower-ranked term sharing >= theme_overlap of the smaller protein set (overlap coefficient).
ovl <- function(x, y) length(intersect(x, y)) / min(length(x), length(y))
themes <- do.call(rbind, lapply(aa_cols, function(a) {
  idx <- which(sig$amino_acid == a & sig$overlap >= rep_min_overlap &
               sig$odds_ratio >= rep_min_or)
  idx <- idx[order(sig$fdr[idx], -sig$odds_ratio[idx])]
  out <- list(); theme_no <- 0
  while (length(idx)) {
    lead <- idx[1]; theme_no <- theme_no + 1
    same <- idx[vapply(idx, function(j)
      ovl(sig_members[[lead]], sig_members[[j]]) >= theme_overlap, logical(1))]
    out[[theme_no]] <- data.frame(
      amino_acid = a, amino_acid_name = aa_names[[a]], theme_rank = theme_no,
      lead_GOID = sig$GOID[lead], lead_term = sig$term[lead],
      overlap = sig$overlap[lead], expected = sig$expected[lead],
      odds_ratio = sig$odds_ratio[lead], or_ci95_low = sig$or_ci95_low[lead],
      or_ci95_high = sig$or_ci95_high[lead], fdr = sig$fdr[lead],
      n_terms_in_theme = length(same),
      related_terms = paste(sig$term[setdiff(same, lead)], collapse = "; "),
      lead_genes = sig$genes[lead], stringsAsFactors = FALSE)
    idx <- setdiff(idx, same)
  }
  do.call(rbind, out)
}))
cat("Representative themes:", nrow(themes), "across",
    length(unique(themes$amino_acid)), "amino acids\n")

# ---- 6. Matrix + row ordering -----------------------------------------------
terms <- unique(sig$GOID)
mat <- matrix(NA_real_, length(terms), length(aa_cols),
              dimnames = list(terms, aa_cols))
mat[cbind(sig$GOID, sig$amino_acid)] <- sig$log2_odds_ratio

row_info <- data.frame(
  GOID      = terms,
  term      = unname(Term(terms)),
  strongest = aa_cols[apply(mat, 1, which.max)],     # which.max ignores NA
  max_log2  = apply(mat, 1, max, na.rm = TRUE),
  n_sig     = rowSums(!is.na(mat)),
  stringsAsFactors = FALSE)
ord <- order(match(row_info$strongest, aa_cols), -row_info$n_sig,
             -row_info$max_log2, row_info$term)
row_info <- row_info[ord, ]
mat      <- mat[ord, , drop = FALSE]
row_info$row_order <- seq_len(nrow(row_info))

# ---- 7. Heatmap --------------------------------------------------------------
col_fun <- colorRamp2(seq(0, color_cap, length.out = 9),
                      brewer.pal(9, "YlGnBu"))      # values > cap drawn as cap
plot_mat <- pmin(mat, color_cap)
row_split <- factor(row_info$strongest, levels = aa_cols)

ht <- Heatmap(
  plot_mat,
  name              = "log2 OR",
  col               = col_fun,
  na_col            = "grey95",
  cluster_rows      = FALSE,
  cluster_columns   = FALSE,
  row_split         = row_split,
  row_title         = NULL,
  row_gap           = unit(0.4, "mm"),
  border            = FALSE,
  show_row_names    = FALSE,
  column_names_rot  = 0,
  column_names_centered = TRUE,
  column_names_gp   = gpar(fontsize = 11),
  column_title      = "Amino acid",
  column_title_side = "bottom",
  column_title_gp   = gpar(fontsize = 12),
  left_annotation   = NULL,
  heatmap_legend_param = list(
    title = "Enrichment\n(log2 odds ratio)", at = seq(0, color_cap, 1.5),
    labels = c(seq(0, color_cap - 1.5, 1.5), paste0("\u2265", color_cap)),
    legend_height = unit(4, "cm"), title_gp = gpar(fontsize = 10)),
  use_raster        = FALSE)

draw_fig <- function() {
  draw(ht, heatmap_legend_side = "right",
       padding = unit(c(2, 8, 2, 2), "mm"))
  decorate_heatmap_body("log2 OR", slice = 1, {
    grid.text(sprintf("GO Biological Process terms (n = %d)", nrow(mat)),
              x = unit(-4, "mm"), y = unit(1, "npc"), rot = 90,
              just = c("right", "bottom"), gp = gpar(fontsize = 12))
  })
}

pdf(file.path(out_dir, "Figure1B_GO_heatmap.pdf"), width = 6, height = 7)
draw_fig(); invisible(dev.off())
png(file.path(out_dir, "Figure1B_GO_heatmap.png"),
    width = 6, height = 7, units = "in", res = 300)
draw_fig(); invisible(dev.off())

# ---- 8. Source data (Table S2) ----------------------------------------------
heatmap_sheet <- data.frame(
  `Row order` = row_info$row_order, `GO ID` = row_info$GOID,
  `GO Biological Process term` = row_info$term,
  `Strongest amino acid association` = row_info$strongest,
  `Maximum log2 odds ratio` = row_info$max_log2,
  `Significant amino acid groups` = row_info$n_sig,
  mat, check.names = FALSE, row.names = NULL)

assoc_sheet <- data.frame(
  `Amino acid` = sig$amino_acid, `Amino acid name` = sig$amino_acid_name,
  `GO ID` = sig$GOID, `GO Biological Process term` = sig$term,
  `Group proteins (annotated)` = sig$group_size,
  `Background proteins` = sig$background,
  `GO term proteins` = sig$term_size, Overlap = sig$overlap,
  Expected = sig$expected, `Odds ratio` = sig$odds_ratio,
  `OR 95% CI low` = sig$or_ci95_low, `OR 95% CI high` = sig$or_ci95_high,
  `Log2 odds ratio` = sig$log2_odds_ratio, `P value` = sig$p_value,
  `Global FDR` = sig$fdr, `Group genes in term` = sig$genes, check.names = FALSE)

theme_sheet <- data.frame(
  `Amino acid` = themes$amino_acid, `Amino acid name` = themes$amino_acid_name,
  `Theme rank` = themes$theme_rank, `Lead GO ID` = themes$lead_GOID,
  `Lead GO term` = themes$lead_term, Overlap = themes$overlap,
  Expected = themes$expected, `Odds ratio` = themes$odds_ratio,
  `OR 95% CI low` = themes$or_ci95_low, `OR 95% CI high` = themes$or_ci95_high,
  `Global FDR` = themes$fdr, `Terms in theme` = themes$n_terms_in_theme,
  `Related terms` = themes$related_terms, `Lead term genes` = themes$lead_genes,
  check.names = FALSE)

summary_sheet <- data.frame(
  `Amino acid` = aa_cols, `Amino acid name` = unname(aa_names[aa_cols]),
  `Group proteins` = lengths(groups[aa_cols]),
  `Group proteins with GO BP` = vapply(aa_cols, function(a)
    length(intersect(groups[[a]], universe)), integer(1)),
  `Terms tested` = vapply(aa_cols, function(a) sum(res$amino_acid == a), integer(1)),
  `Significant GO BP terms` = vapply(aa_cols, function(a)
    sum(sig$amino_acid == a), integer(1)),
  check.names = FALSE, row.names = NULL)

params_sheet <- data.frame(
  Parameter = c("Z-score input", "Grouping", "Proteins analysed", "Background (GO BP-annotated)",
                "GO annotation", "Term size range", "Minimum overlap",
                "Odds-ratio correction", "Multiple testing", "Tests performed",
                "Significant associations", "Distinct GO terms",
                "Colour scale cap", "Representative theme criteria",
                "Representative themes", "org.Hs.eg.db", "GO.db",
                "GO source date", "Run date"),
  Value = c(basename(zscore_tsv), grouping, n_proteome, N,
            if (go_propagate) "propagated (GOALL)" else "direct (GO)",
            paste(min_term_size, max_term_size, sep = "-"), min_overlap,
            paste("Haldane-Anscombe +", haldane),
            paste("Benjamini-Hochberg, global; FDR <", fdr_cutoff), nrow(res),
            nrow(sig), length(terms), color_cap,
            sprintf("FDR < %g, overlap >= %d, OR >= %g; ranked by FDR; terms merged when group-protein overlap coefficient >= %g",
                    fdr_cutoff, rep_min_overlap, rep_min_or, theme_overlap),
            nrow(themes),
            as.character(packageVersion("org.Hs.eg.db")),
            as.character(packageVersion("GO.db")),
            subset(metadata(GO.db), name == "GOSOURCEDATE")$value,
            as.character(Sys.Date())))

# One row per amino acid: top 3 themes, for the main text
top_themes_sheet <- do.call(rbind, lapply(aa_cols, function(a) {
  t <- head(themes[themes$amino_acid == a, ], 3)
  data.frame(`Amino acid` = a, `Amino acid name` = aa_names[[a]],
             `Representative themes` = nrow(themes[themes$amino_acid == a, ]),
             `Top themes (lead term: overlap, OR [95% CI], FDR)` = if (nrow(t))
               paste(sprintf("%s: n = %d, OR %.1f [%.1f-%.1f], FDR %.1e", t$lead_term,
                             t$overlap, t$odds_ratio, t$or_ci95_low, t$or_ci95_high,
                             t$fdr), collapse = " | ") else "none",
             check.names = FALSE)
}))

sheets <- list(`Heatmap matrix` = heatmap_sheet, Associations = assoc_sheet,
               `Representative themes` = theme_sheet, `Top themes` = top_themes_sheet,
               `AA summary` = summary_sheet, Parameters = params_sheet,
               `All tests` = res)
# .xlsx keeps very small P values as numbers (Excel can read them from CSV as
# text, which number filters then silently drop); CSVs are for scripts.
if (requireNamespace("writexl", quietly = TRUE))
  writexl::write_xlsx(sheets, file.path(out_dir, "TableS2_Fig1B_source_data.xlsx"))
for (s in names(sheets)) {
  f <- file.path(out_dir, paste0("TableS2_", gsub(" ", "_", s), ".csv"))
  tryCatch(write.csv(sheets[[s]], f, row.names = FALSE),
           error = function(e) warning("Could not write ", basename(f),
                                       " (open in Excel?); skipped", call. = FALSE))
}

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
print(summary_sheet, row.names = FALSE)
cat("Outputs written to", out_dir, "\n")
