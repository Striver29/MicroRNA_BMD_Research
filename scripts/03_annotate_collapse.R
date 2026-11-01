#!/usr/bin/env Rscript
# 03_annotate_collapse.R -- probes -> HGNC symbols, one probe per gene.  CLAUDE.md s6 step 3.
suppressPackageStartupMessages({library(data.table); library(AnnotationDbi); library(hgu133plus2.db)})
set.seed(1)

x <- readRDS("data/processed/rma.rds"); mat <- x$mat; sheet <- x$sheet
cat("RMA matrix in:", nrow(mat), "probes x", ncol(mat), "samples\n")

map <- suppressMessages(AnnotationDbi::select(hgu133plus2.db, keys = rownames(mat),
                                              keytype = "PROBEID", columns = "SYMBOL"))
map <- as.data.table(map)[!is.na(SYMBOL)]
n_probe_all <- nrow(mat)
n_mapped <- uniqueN(map$PROBEID)

# a probe that hits >1 symbol is ambiguous -> drop it (CLAUDE.md step 3)
multi <- map[, .N, by = PROBEID][N > 1, PROBEID]
map <- map[!PROBEID %in% multi]
cat("probes total           :", n_probe_all, "\n")
cat("probes with a symbol   :", n_mapped, "\n")
cat("probes dropped, no sym :", n_probe_all - n_mapped, "\n")
cat("probes dropped, multi  :", length(multi), "\n")
cat("probes kept            :", uniqueN(map$PROBEID), "\n")

# collapse: per gene keep the probe with the highest mean expression (CLAUDE.md s10 default)
rm_ <- rowMeans(mat)
map[, mean_expr := rm_[PROBEID]]
setorder(map, SYMBOL, -mean_expr)
best <- map[, .SD[1], by = SYMBOL]
cat("genes kept             :", nrow(best), "\n")
cat("genes with >1 probe    :", sum(map[, .N, by = SYMBOL]$N > 1), "\n")

gmat <- mat[best$PROBEID, , drop = FALSE]
rownames(gmat) <- best$SYMBOL
gmat <- gmat[order(rownames(gmat)), , drop = FALSE]
cat("gene matrix out        :", nrow(gmat), "x", ncol(gmat), "\n")

saveRDS(list(mat = gmat, sheet = sheet, probe_map = best), "data/processed/gene_matrix.rds")
cat("Wrote data/processed/gene_matrix.rds\n")

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
fwrite(data.table(step = c("probes_total","probes_with_symbol","probes_dropped_unmapped",
                           "probes_dropped_multimap","probes_kept","genes_kept"),
                  n = c(n_probe_all, n_mapped, n_probe_all - n_mapped,
                        length(multi), uniqueN(map$PROBEID), nrow(best))),
       "results/tables/03_collapse_counts.tsv", sep = "\t")
fwrite(best[, .(gene_symbol = SYMBOL, probe_id = PROBEID, mean_expr = round(mean_expr, 4))],
       "results/tables/03_gene_probe_map.tsv", sep = "\t")
cat("Wrote results/tables/03_collapse_counts.tsv and 03_gene_probe_map.tsv\n")
