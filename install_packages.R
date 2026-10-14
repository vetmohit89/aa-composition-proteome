# Install the R packages used by this repository (run once).
#   Rscript install_packages.R
# Tested with R 4.5.3 and Bioconductor 3.22 (org.Hs.eg.db 3.22.0, GO.db 3.22.0).

cran <- c("ggplot2", "ggrepel", "patchwork", "circlize", "RColorBrewer",
          "readxl", "writexl", "BiocManager")
bioc <- c("AnnotationDbi", "org.Hs.eg.db", "GO.db", "ComplexHeatmap")

repos <- "https://cloud.r-project.org"
missing <- setdiff(cran, rownames(installed.packages()))
if (length(missing)) install.packages(missing, repos = repos)

missing <- setdiff(bioc, rownames(installed.packages()))
if (length(missing)) BiocManager::install(missing, update = FALSE, ask = FALSE)

cat("Installed versions:\n")
for (p in c(cran, bioc)) cat(sprintf("  %-15s %s\n", p, as.character(packageVersion(p))))
