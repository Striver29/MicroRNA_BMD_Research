#!/usr/bin/env Rscript
# 10a_ecm_removed.R -- is the footprint just the scar genes?  CLAUDE.md s6 step 10a.
# Drop every matrisome-core gene from the ranked list and rerun fgsea.
suppressPackageStartupMessages({library(data.table); library(fgsea); library(ggplot2)})
set.seed(42)
pick <- function(f) { a <- file.path("inputs", f); if (file.exists(a)) a else file.path("inputs/auto", f) }

fs  <- readRDS("results/rds/family_sets.rds")
lim <- readRDS("results/rds/limma.rds")
fib <- fread(pick("fibrosis_sets.tsv"))
ecm <- unique(fib[set_name == "matrisome_core", gene_symbol])
cat("matrisome-core genes:", length(ecm), "\n")

contrasts <- c("BMDvsCtrl", "DMDvsCtrl")
out <- rbindlist(lapply(contrasts, function(cn) {
  tt <- lim$tables[[cn]]
  on_array <- intersect(ecm, tt$gene_symbol)
  tt2 <- tt[!gene_symbol %in% ecm]
  cat(sprintf("%s: %d genes -> %d after dropping %d matrisome-core genes present on the array\n",
              cn, nrow(tt), nrow(tt2), length(on_array)))
  # sets must also lose those genes, else sizes are wrong
  fs2 <- lapply(fs, function(g) setdiff(g, ecm))
  fs2 <- fs2[lengths(fs2) >= 15]
  rk <- sort(setNames(tt2$t, tt2$gene_symbol), decreasing = TRUE)
  r <- as.data.table(fgsea(fs2, rk, minSize = 15, maxSize = Inf, nPermSimple = 10000, eps = 0))
  r[, contrast := cn][]
}))

real <- readRDS("results/rds/fgsea_all.rds")
cmp <- merge(real[contrast %in% contrasts, .(contrast, family_id = pathway,
                                             size_full = size, NES_full = NES, padj_full = padj)],
             out[, .(contrast, family_id = pathway, size_noECM = size,
                     NES_noECM = NES, padj_noECM = padj)],
             by = c("contrast", "family_id"))
cmp[, survives := padj_noECM < 0.05 & sign(NES_noECM) == sign(NES_full)]
setorder(cmp, contrast, padj_full)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
fwrite(cmp, "results/tables/10a_ecm_removed.tsv", sep = "\t")

cat("\n=============== ECM REMOVED ===============\n")
cat("Question: if we delete the scar genes themselves, is there still a footprint?\n\n")
print(cmp[, .(n_families = .N,
              sig_full = sum(padj_full < 0.05),
              sig_noECM = sum(padj_noECM < 0.05),
              survive = sum(survives)), by = contrast])
cat("\n== miR-29, with and without the matrisome-core genes ==\n")
print(cmp[family_id == "miR-29-3p",
          .(contrast, size_full, size_noECM,
            NES_full = round(NES_full, 3), padj_full = signif(padj_full, 2),
            NES_noECM = round(NES_noECM, 3), padj_noECM = signif(padj_noECM, 2), survives)])
cat("\n== top 15 DMD families, before vs after ==\n")
print(head(cmp[contrast == "DMDvsCtrl",
               .(family_id, size_full, size_noECM, NES_full = round(NES_full, 2),
                 NES_noECM = round(NES_noECM, 2), padj_noECM = signif(padj_noECM, 2), survives)], 15))

d <- cmp[contrast == "DMDvsCtrl"]
cat(sprintf("\nNES shrinkage across families (DMD vs ctrl): median %.3f -> %.3f (%.1f%% of the original)\n",
            median(d$NES_full), median(d$NES_noECM), 100 * median(d$NES_noECM) / median(d$NES_full)))
m <- cmp[family_id == "miR-29-3p" & contrast == "DMDvsCtrl"]
cat(sprintf("miR-29 specifically: NES %.3f -> %.3f (%.1f%% retained)\n",
            m$NES_full, m$NES_noECM, 100 * m$NES_noECM / m$NES_full))
saveRDS(cmp, "results/rds/ecm_removed.rds")
