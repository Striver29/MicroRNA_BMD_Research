#!/usr/bin/env Rscript
# 09_bmd_vs_dmd.R -- the direct comparison, where the shared controls cancel out.
# CLAUDE.md s6 step 9.  Shuffles BMD/DMD labels among the 28 patients only (11/17).
suppressPackageStartupMessages({library(data.table); library(limma); library(fgsea)
                                library(BiocParallel); library(ggplot2)})
set.seed(20261102)
N_SHUFFLE  <- as.integer(Sys.getenv("N_SHUFFLE", "1000"))
NPERM_NULL <- as.integer(Sys.getenv("NPERM_NULL", "10000"))
N_CORE     <- as.integer(Sys.getenv("N_CORE", "7"))
TIME_ONLY  <- as.integer(Sys.getenv("TIME_ONLY", "0"))

fs  <- readRDS("results/rds/family_sets.rds")
gm  <- readRDS("data/processed/gene_matrix.rds")
mat <- gm$mat; sheet <- as.data.table(gm$sheet); setkey(sheet, gsm); sheet <- sheet[colnames(mat)]

pat <- sheet$group %in% c("BMD", "DMD")
mp  <- mat[, pat, drop = FALSE]
gp0 <- sheet$group[pat]
cat("patients only:", ncol(mp), " BMD:", sum(gp0 == "BMD"), " DMD:", sum(gp0 == "DMD"), "\n")
cat("NOTE: the 6 controls are not used here at all, so a fluke in the controls cannot\n")
cat("      create a BMD-DMD difference. This is the comparison that answers\n")
cat("      'is BMD different from DMD?'.\n")

real <- readRDS("results/rds/fgsea_all.rds")[contrast == "BMDvsDMD"]
real_nes <- setNames(real$NES, real$pathway)[names(fs)]

one <- function(i) {
  set.seed(200000L + i)
  g <- factor(sample(gp0), levels = c("DMD", "BMD"))
  d <- model.matrix(~ 0 + g); colnames(d) <- levels(g)
  f <- eBayes(contrasts.fit(lmFit(mp, d), makeContrasts(BMDvsDMD = BMD - DMD, levels = d)))
  tt <- topTable(f, coef = "BMDvsDMD", number = Inf, sort.by = "none")
  rk <- sort(setNames(tt$t, rownames(tt)), decreasing = TRUE)
  r <- fgsea(fs, rk, minSize = 15, maxSize = Inf, nPermSimple = NPERM_NULL, eps = 0)
  setNames(r$NES, r$pathway)[names(fs)]
}

t0 <- Sys.time(); invisible(lapply(1:5, one)); per <- as.numeric(difftime(Sys.time(), t0, units="secs"))/5
cat(sprintf("\nTiming: %.2f s per shuffle; projected %.1f min for %d on %d cores\n",
            per, N_SHUFFLE*per/N_CORE/60, N_SHUFFLE, N_CORE))
if (TIME_ONLY == 1L) quit(status = 0)

bp <- MulticoreParam(workers = N_CORE, RNGseed = 20261102)
cat(sprintf("Running %d patient-label shuffles ...\n", N_SHUFFLE))
t1 <- Sys.time()
nullm <- do.call(rbind, bplapply(seq_len(N_SHUFFLE), one, BPPARAM = bp))
cat(sprintf("Done in %.1f min\n", as.numeric(difftime(Sys.time(), t1, units = "mins"))))

res <- rbindlist(lapply(names(fs), function(f) {
  v <- nullm[, f]; v <- v[is.finite(v)]; rn <- real_nes[[f]]
  data.table(family_id = f, size = length(fs[[f]]), NES_BMDvsDMD = rn,
             fgsea_padj = real[pathway == f, padj],
             shuffle_p = (sum(abs(v) >= abs(rn)) + 1) / (length(v) + 1),
             null_mean = mean(v), null_sd = sd(v), n_null = length(v))
}))
res[, shuffle_q := p.adjust(shuffle_p, "BH")]
setorder(res, shuffle_p, -abs(NES_BMDvsDMD))
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
fwrite(res, "results/tables/09_bmd_vs_dmd.tsv", sep = "\t")

cat("\n=============== DIRECT BMD vs DMD ===============\n")
cat("Reading the sign: NES < 0 means the family's targets are LOWER in BMD than in DMD,\n")
cat("so the brake is predicted MORE active in BMD -- i.e. less weakened than in DMD.\n\n")
print(head(res[, .(family_id, size, NES = round(NES_BMDvsDMD, 3),
                   fgsea_padj = signif(fgsea_padj, 2),
                   shuffle_p = signif(shuffle_p, 3), shuffle_q = signif(shuffle_q, 3))], 20))
cat("\nfamilies with fgsea padj < 0.05        :", res[fgsea_padj < 0.05, .N], "of", nrow(res), "\n")
cat("families with label-shuffle p < 0.05   :", res[shuffle_p < 0.05, .N], "\n")
cat("families passing BOTH                  :", res[fgsea_padj < 0.05 & shuffle_p < 0.05, .N], "\n")
cat("\n== miR-29 ==\n"); print(res[family_id == "miR-29-3p"])

## ---- shared / BMD-only / DMD-only, from HERE (CLAUDE.md step 9) -----------
w <- fread("results/tables/07_fgsea_wide.tsv")
sh <- fread("results/tables/08_family_shuffle_p.tsv", showProgress = FALSE)
cls <- merge(w[, .(family_id, NES_BMDvsCtrl, padj_BMDvsCtrl, NES_DMDvsCtrl, padj_DMDvsCtrl)],
             res[, .(family_id, NES_BMDvsDMD, diff_p = shuffle_p, diff_q = shuffle_q)], by = "family_id")
cls[, sig_BMD := padj_BMDvsCtrl < 0.05]
cls[, sig_DMD := padj_DMDvsCtrl < 0.05]
cls[, differs := diff_q < 0.05]
cls[, class := fifelse(!differs & sig_BMD & sig_DMD, "shared",
                fifelse(differs & sig_DMD & !sig_BMD, "DMD-only",
                fifelse(differs & sig_BMD & !sig_DMD, "BMD-only",
                fifelse(differs & sig_BMD & sig_DMD, "both, but different size",
                        "neither"))))]
cat("\n== classification (differs = label-shuffle q < 0.05 in the DIRECT comparison) ==\n")
print(cls[, .N, by = class][order(-N)])
fwrite(cls, "results/tables/09_family_classification.tsv", sep = "\t")
cat("\n== families that differ between BMD and DMD ==\n")
print(cls[differs == TRUE][order(diff_q), .(family_id, NES_BMD = round(NES_BMDvsCtrl,2),
          NES_DMD = round(NES_DMDvsCtrl,2), NES_BvD = round(NES_BMDvsDMD,2),
          diff_q = signif(diff_q,3), class)])

saveRDS(list(res = res, nullm = nullm, cls = cls), "results/rds/bmd_vs_dmd.rds")

## ---- figure: null NES vs real NES for miR-29 ------------------------------
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
v <- nullm[, "miR-29-3p"]; v <- v[is.finite(v)]
rn <- real_nes[["miR-29-3p"]]
p <- ggplot(data.table(v = v), aes(v)) +
  geom_histogram(bins = 40, fill = "grey75", colour = "white") +
  geom_vline(xintercept = rn, colour = "#B2182B", linewidth = 1) +
  annotate("label", x = rn, y = Inf, vjust = 1.3, size = 2.8, colour = "#B2182B",
           label = sprintf("real NES = %.2f\nempirical p = %.3g", rn,
                           (sum(abs(v) >= abs(rn)) + 1) / (length(v) + 1))) +
  labs(x = "miR-29 NES for BMD vs DMD", y = sprintf("number of shuffles (of %d)", length(v)),
       title = "miR-29 in the direct BMD-vs-DMD comparison, against shuffled patient labels",
       subtitle = "Grey = BMD/DMD labels shuffled among the 28 patients keeping 11/17.\nRed = the value we measured. Controls play no part in this test.") +
  theme_bw(base_size = 9)
ggsave("results/figures/09_mir29_bmd_vs_dmd_shuffle.png", p, width = 6.5, height = 4.2, dpi = 200)
cat("Wrote results/figures/09_mir29_bmd_vs_dmd_shuffle.png\n")
