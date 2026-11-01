#!/usr/bin/env Rscript
# 14_subsample_dmd.R -- is BMD's weaker signal just a sample-size effect?
# CLAUDE.md s6 step 14. Cut DMD down to 11 samples, 100 draws, and compare.
suppressPackageStartupMessages({library(data.table); library(limma); library(fgsea)
                                library(BiocParallel); library(ggplot2)})
set.seed(20261106)
N_DRAW <- as.integer(Sys.getenv("N_DRAW", "100"))
N_CORE <- as.integer(Sys.getenv("N_CORE", "7"))

fs  <- readRDS("results/rds/family_sets.rds")
gm  <- readRDS("data/processed/gene_matrix.rds")
mat <- gm$mat; sheet <- as.data.table(gm$sheet); setkey(sheet, gsm); sheet <- sheet[colnames(mat)]
dmd <- sheet[group == "DMD", gsm]; ctl <- sheet[group == "control", gsm]
real <- readRDS("results/rds/fgsea_all.rds")
nes_bmd <- real[contrast == "BMDvsCtrl", setNames(NES, pathway)][names(fs)]
nes_dmd <- real[contrast == "DMDvsCtrl", setNames(NES, pathway)][names(fs)]

one <- function(i) {
  set.seed(500000L + i)
  pick <- c(ctl, sample(dmd, 11))
  m <- mat[, pick]
  g <- factor(sheet[match(pick, gsm), group], levels = c("control", "DMD"))
  d <- model.matrix(~ 0 + g); colnames(d) <- levels(g)
  f <- eBayes(contrasts.fit(lmFit(m, d), makeContrasts(DMDvsCtrl = DMD - control, levels = d)))
  tt <- topTable(f, coef = "DMDvsCtrl", number = Inf, sort.by = "none")
  rk <- sort(setNames(tt$t, rownames(tt)), decreasing = TRUE)
  r <- as.data.table(fgsea(fs, rk, minSize = 15, maxSize = Inf, nPermSimple = 10000, eps = 0, nproc = 1))
  list(nes = setNames(r$NES, r$pathway)[names(fs)],
       n_sig = r[padj < 0.05, .N],
       m29_nes = r[pathway == "miR-29-3p", NES], m29_padj = r[pathway == "miR-29-3p", padj])
}

bp <- MulticoreParam(workers = N_CORE, RNGseed = 20261106)
cat(sprintf("Running %d draws of DMD down-sampled to 11 (controls fixed) ...\n", N_DRAW))
t0 <- Sys.time()
out <- bplapply(seq_len(N_DRAW), one, BPPARAM = bp)
cat(sprintf("done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

nm <- do.call(rbind, lapply(out, `[[`, "nes"))
summ <- data.table(draw = seq_len(N_DRAW),
                   n_sig = vapply(out, `[[`, numeric(1), "n_sig"),
                   m29_nes = vapply(out, `[[`, numeric(1), "m29_nes"),
                   m29_padj = vapply(out, `[[`, numeric(1), "m29_padj"))
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
fwrite(summ, "results/tables/14_subsample_dmd.tsv", sep = "\t")

bmd_sig <- real[contrast == "BMDvsCtrl" & padj < 0.05, .N]
dmd_sig <- real[contrast == "DMDvsCtrl" & padj < 0.05, .N]
cat("\n=============== DMD SUBSAMPLED TO n = 11 ===============\n")
cat(sprintf("families significant, full DMD (n=17)      : %d of %d\n", dmd_sig, length(fs)))
cat(sprintf("families significant, DMD cut to 11        : median %.0f (range %d-%d)\n",
            median(summ$n_sig), min(summ$n_sig), max(summ$n_sig)))
cat(sprintf("families significant, actual BMD (n=11)    : %d\n", bmd_sig))
cat(sprintf("\nmiR-29 NES, full DMD                       : %.3f\n", nes_dmd[["miR-29-3p"]]))
cat(sprintf("miR-29 NES, DMD cut to 11                  : median %.3f (range %.3f-%.3f)\n",
            median(summ$m29_nes), min(summ$m29_nes), max(summ$m29_nes)))
cat(sprintf("miR-29 NES, actual BMD (n=11)              : %.3f\n", nes_bmd[["miR-29-3p"]]))
pct <- mean(summ$m29_nes <= nes_bmd[["miR-29-3p"]])
cat(sprintf("\nShare of 11-sample DMD draws whose miR-29 NES is at or below BMD's: %.1f%%\n", 100*pct))
cat("Read it like this: if that share is large, BMD's smaller miR-29 footprint is within what\n")
cat("DMD itself produces at n = 11, so the difference may be sample size rather than biology.\n")
cat("If it is small, BMD really is weaker than DMD even at matched n.\n")

p <- ggplot(summ, aes(m29_nes)) +
  geom_histogram(bins = 30, fill = "grey75", colour = "white") +
  geom_vline(xintercept = nes_bmd[["miR-29-3p"]], colour = "#B2182B", linewidth = 1) +
  geom_vline(xintercept = nes_dmd[["miR-29-3p"]], colour = "#2166AC", linewidth = 1, linetype = 2) +
  labs(x = "miR-29 NES, DMD vs control", y = sprintf("number of draws (of %d)", N_DRAW),
       title = "Does BMD look weaker only because there are fewer BMD samples?",
       subtitle = "Grey = DMD randomly cut to 11 samples. Red solid = the real BMD value.\nBlue dashed = the full 17-sample DMD value.") +
  theme_bw(base_size = 9)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("results/figures/14_subsample_dmd.png", p, width = 6.5, height = 4.2, dpi = 200)
cat("Wrote results/figures/14_subsample_dmd.png\n")
saveRDS(list(summ = summ, nes = nm), "results/rds/subsample_dmd.rds")
