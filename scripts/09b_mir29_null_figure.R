#!/usr/bin/env Rscript
# 09b_mir29_null_figure.R -- regenerate step 9's miR-29 figure, and the family
# classification, without repeating the full 128-family shuffle.
#
# Why this is legitimate: fgsea normalises each pathway against ITS OWN permutation
# distribution, so a pathway's NES is identical whether it is scored alone or alongside 127
# others. Only padj (a BH correction across pathways) would differ, and the null distribution
# does not use padj. Running miR-29 alone is therefore exact and ~100x cheaper.
suppressPackageStartupMessages({library(data.table); library(limma); library(fgsea)
                                library(BiocParallel); library(ggplot2)})
set.seed(20261102)                      # same seed as step 9
N_SHUFFLE <- as.integer(Sys.getenv("N_SHUFFLE", "500"))
N_CORE    <- as.integer(Sys.getenv("N_CORE", "4"))

fs  <- readRDS("results/rds/family_sets.rds")
gm  <- readRDS("data/processed/gene_matrix.rds")
mat <- gm$mat; sheet <- as.data.table(gm$sheet); setkey(sheet, gsm); sheet <- sheet[colnames(mat)]
pat <- sheet$group %in% c("BMD", "DMD")
mp  <- mat[, pat, drop = FALSE]; gp0 <- sheet$group[pat]
res9 <- fread("results/tables/09_bmd_vs_dmd.tsv")
rn   <- res9[family_id == "miR-29-3p", NES_BMDvsDMD]
one_fam <- fs["miR-29-3p"]

one <- function(i) {
  set.seed(200000L + i)                 # identical draws to step 9
  g <- factor(sample(gp0), levels = c("DMD", "BMD"))
  d <- model.matrix(~ 0 + g); colnames(d) <- levels(g)
  f <- eBayes(contrasts.fit(lmFit(mp, d), makeContrasts(BMDvsDMD = BMD - DMD, levels = d)))
  tt <- topTable(f, coef = "BMDvsDMD", number = Inf, sort.by = "none")
  rk <- sort(setNames(tt$t, rownames(tt)), decreasing = TRUE)
  fgsea(one_fam, rk, minSize = 15, maxSize = Inf, nPermSimple = 10000, eps = 0, nproc = 1)$NES[1]
}
cat(sprintf("Running %d miR-29-only patient-label shuffles on %d cores ...\n", N_SHUFFLE, N_CORE))
t0 <- Sys.time()
v <- unlist(bplapply(seq_len(N_SHUFFLE), one, BPPARAM = MulticoreParam(workers = N_CORE, RNGseed = 20261102)))
cat(sprintf("done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
v <- v[is.finite(v)]
p_emp <- (sum(abs(v) >= abs(rn)) + 1) / (length(v) + 1)
cat(sprintf("\nreal NES %.3f | null mean %.3f sd %.3f | empirical p = %.4f (n=%d)\n",
            rn, mean(v), sd(v), p_emp, length(v)))
cat(sprintf("step 9's full run reported p = %.4f -- these should agree closely.\n",
            res9[family_id == "miR-29-3p", shuffle_p]))
saveRDS(list(null_nes = v, real_nes = rn, p = p_emp), "results/rds/mir29_bmdvsdmd_null.rds")

p <- ggplot(data.table(v = v), aes(v)) +
  geom_histogram(bins = 40, fill = "grey75", colour = "white") +
  geom_vline(xintercept = rn, colour = "#B2182B", linewidth = 1) +
  annotate("label", x = rn, y = Inf, vjust = 1.3, size = 2.9, colour = "#B2182B",
           label = sprintf("real NES = %.2f\nempirical p = %.2f\nNOT significant", rn, p_emp)) +
  labs(x = "miR-29 NES for BMD vs DMD", y = sprintf("number of shuffles (of %d)", length(v)),
       title = "miR-29, BMD vs DMD, against shuffled PATIENT labels",
       subtitle = paste0("Grey = BMD/DMD labels shuffled among the 28 patients, keeping 11/17.\n",
                         "Red = the value we measured. It sits INSIDE the grey bulk, so the\n",
                         "BMD-vs-DMD difference is not distinguishable from relabelling noise.")) +
  theme_bw(base_size = 9)
ggsave("results/figures/09_mir29_bmd_vs_dmd_shuffle.png", p, width = 6.5, height = 4.4, dpi = 200)
cat("Wrote results/figures/09_mir29_bmd_vs_dmd_shuffle.png\n")
