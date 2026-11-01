#!/usr/bin/env Rscript
# 00_setup.R -- install everything the pipeline needs, record the environment.
# Per CLAUDE.md s6 step 0.

options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  Ncpus = max(1L, parallel::detectCores() - 1L),
  timeout = 3600
)

dir.create("results", showWarnings = FALSE, recursive = TRUE)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

bioc_pkgs <- c(
  "GEOquery", "oligo", "pd.hg.u133.plus.2", "hgu133plus2.db",
  "AnnotationDbi", "limma", "fgsea", "msigdbr", "BiocParallel"
)
cran_pkgs <- c("data.table", "ggplot2")

need <- function(p) !requireNamespace(p, quietly = TRUE)

to_get <- c(cran_pkgs, bioc_pkgs)
to_get <- to_get[vapply(to_get, need, logical(1))]

if (length(to_get)) {
  cat("Installing:", paste(to_get, collapse = ", "), "\n")
  BiocManager::install(to_get, update = FALSE, ask = FALSE)
}

missing <- c(cran_pkgs, bioc_pkgs)[vapply(c(cran_pkgs, bioc_pkgs), need, logical(1))]
if (length(missing)) {
  cat("STILL MISSING:", paste(missing, collapse = ", "), "\n")
} else {
  cat("All packages present.\n")
}

cat("Bioconductor version:", as.character(BiocManager::version()), "\n")
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
cat("Wrote results/sessionInfo.txt\n")
