#!/usr/bin/env Rscript
# 11_utr_matched.R -- 3'UTR-length-matched null.  CLAUDE.md s6 step 11.
#
# Why: a gene with a long 3'UTR has room for many microRNA sites, so it appears on many
# target lists. Long-UTR genes also behave alike in disease. A plain random gene set is
# therefore not a fair comparison for a real target set. We instead draw random sets that
# copy the family's 3'UTR-length profile, decile by decile.
#
# Worked example with 4 deciles and a family of 10 genes:
#   family has 1 gene in D1, 2 in D2, 3 in D3, 4 in D4
#   each random set must also take exactly 1 from D1, 2 from D2, 3 from D3, 4 from D4
suppressPackageStartupMessages({library(data.table); library(fgsea); library(BiocParallel)})
set.seed(20261103)
N_DRAW    <- as.integer(Sys.getenv("N_DRAW", "1000"))
N_CORE    <- as.integer(Sys.getenv("N_CORE", "7"))
N_TOP     <- as.integer(Sys.getenv("N_TOP", "25"))
pick <- function(f) { a <- file.path("inputs", f); if (file.exists(a)) a else file.path("inputs/auto", f) }

fs   <- readRDS("results/rds/family_sets.rds")
lim  <- readRDS("results/rds/limma.rds")
real <- readRDS("results/rds/fgsea_all.rds")
utr  <- fread(pick("utr_lengths.tsv"))

gm <- readRDS("data/processed/gene_matrix.rds")$mat
genes <- rownames(gm)
u <- utr[gene_symbol %in% genes, .(gene_symbol, utr3_length_nt)]
cat("genes on the array:", length(genes), " with a 3'UTR length:", nrow(u), "\n")
u[, decile := cut(utr3_length_nt, breaks = quantile(utr3_length_nt, probs = seq(0, 1, 0.1)),
                  include.lowest = TRUE, labels = FALSE)]
cat("\n== 3'UTR length deciles ==\n")
print(u[, .(n = .N, min_nt = min(utr3_length_nt), max_nt = max(utr3_length_nt)), by = decile][order(decile)])
by_dec <- split(u$gene_symbol, u$decile)

# which families to test
top_fams <- unique(c("miR-29-3p",
                     head(real[contrast == "DMDvsCtrl"][order(padj, -abs(NES)), pathway], N_TOP)))
top_fams <- intersect(top_fams, names(fs))
cat("\ntesting", length(top_fams), "families against", N_DRAW, "UTR-matched sets each\n")

contrasts <- c("BMDvsCtrl", "DMDvsCtrl")
ranks <- lapply(contrasts, function(cn) {
  tt <- lim$tables[[cn]]; sort(setNames(tt$t, tt$gene_symbol), decreasing = TRUE) })
names(ranks) <- contrasts

# per-family decile profile, restricted to genes that have a UTR length
profiles <- lapply(fs[top_fams], function(g) {
  gg <- intersect(g, u$gene_symbol)
  list(genes = gg, counts = table(factor(u[gene_symbol %in% gg, decile], levels = 1:10)))
})

draw_sets <- function(counts, n) lapply(seq_len(n), function(i) {
  unlist(lapply(names(counts), function(d) {
    k <- counts[[d]]; if (k == 0) character(0) else sample(by_dec[[d]], k)
  }), use.names = FALSE)
})

bp <- MulticoreParam(workers = N_CORE, RNGseed = 20261103)
t0 <- Sys.time()
res <- rbindlist(bplapply(top_fams, function(f) {
  set.seed(300000L + which(top_fams == f))
  pr <- profiles[[f]]
  sets <- draw_sets(pr$counts, N_DRAW)
  names(sets) <- paste0("null", seq_along(sets))
  rbindlist(lapply(contrasts, function(cn) {
    rk <- ranks[[cn]]
    # real NES recomputed on the UTR-restricted gene set so it is comparable to the nulls
    rr <- fgsea(list(real = pr$genes), rk, minSize = 15, maxSize = Inf, nPermSimple = 10000, eps = 0)
    nn <- fgsea(sets, rk, minSize = 15, maxSize = Inf, nPermSimple = 1000, eps = 0)
    v <- nn$NES[is.finite(nn$NES)]
    data.table(family_id = f, contrast = cn, size_utr = length(pr$genes),
               NES_real = rr$NES[1], padj_fgsea = rr$padj[1],
               null_mean = mean(v), null_sd = sd(v), n_null = length(v),
               utr_matched_p = (sum(abs(v) >= abs(rr$NES[1])) + 1) / (length(v) + 1))
  }))
}, BPPARAM = bp))
cat(sprintf("done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

res[, utr_matched_q := p.adjust(utr_matched_p, "BH"), by = contrast]
setorder(res, contrast, utr_matched_p, -abs(NES_real))
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
fwrite(res, "results/tables/11_utr_matched.tsv", sep = "\t")

cat("\n=============== 3'UTR-LENGTH-MATCHED NULL ===============\n")
for (cn in contrasts) {
  cat("\n--", cn, "--\n")
  print(res[contrast == cn, .(family_id, size_utr, NES_real = round(NES_real, 2),
            padj_fgsea = signif(padj_fgsea, 2), null_mean = round(null_mean, 2),
            utr_p = signif(utr_matched_p, 3), utr_q = signif(utr_matched_q, 3))])
  cat("families surviving the matched null at q < 0.05:",
      res[contrast == cn & utr_matched_q < 0.05, .N], "of", res[contrast == cn, .N], "\n")
}
cat("\nRead null_mean like this: if a UTR-matched random set already scores NES ~ +1.5,\n")
cat("then a real family at NES +2.0 is only a little above the background, not 2.0 above zero.\n")
saveRDS(res, "results/rds/utr_matched.rds")
