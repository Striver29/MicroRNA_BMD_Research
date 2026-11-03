#!/usr/bin/env Rscript
# check_duplicate_arrays.R -- before merging any GEO accession, prove its arrays are new.
# Usage: Rscript scripts/check_duplicate_arrays.R <new.CEL> [more.CEL ...]
# GEO re-deposits the same physical arrays under different accessions (GSE38417 vs GSE109178).
# Metadata will not reveal it; identical raw intensities will.
suppressPackageStartupMessages(library(affxparser))
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("give one or more CEL files to check")
ours <- list.files("data/raw/CEL", pattern = "\\.CEL$", full.names = TRUE)
cat("checking", length(args), "candidate(s) against", length(ours), "arrays we already hold\n\n")
for (f in args) {
  x <- readCel(f)$intensities
  best_r <- -1; best <- ""; best_diff <- NA_integer_
  for (o in ours) {
    y <- readCel(o)$intensities
    if (length(x) != length(y)) next
    r <- cor(x, y)
    if (r > best_r) { best_r <- r; best <- basename(o); best_diff <- sum(x != y) }
  }
  verdict <- if (best_diff == 0) "*** DUPLICATE - DO NOT MERGE ***" else "independent, safe to use"
  cat(sprintf("%-28s best match %-40s r=%.6f differing=%d  %s\n",
              basename(f), best, best_r, best_diff, verdict))
}
