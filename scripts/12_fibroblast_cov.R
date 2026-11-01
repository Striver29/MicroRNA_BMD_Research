#!/usr/bin/env Rscript
# 12_fibroblast_cov.R -- does the footprint survive adjusting for how much fibroblast is
# in the biopsy?  CLAUDE.md s6 step 12.
suppressPackageStartupMessages({library(data.table); library(limma); library(fgsea); library(ggplot2)})
set.seed(42)
pick <- function(f) { a <- file.path("inputs", f); if (file.exists(a)) a else file.path("inputs/auto", f) }

fs  <- readRDS("results/rds/family_sets.rds")
gm  <- readRDS("data/processed/gene_matrix.rds")
mat <- gm$mat; sheet <- as.data.table(gm$sheet); setkey(sheet, gsm); sheet <- sheet[colnames(mat)]
fmk <- fread(pick("fibroblast_markers.tsv"))

## ---- the score ------------------------------------------------------------
# Only markers that are NOT in a fibrosis set, so the score is "how much fibroblast is here",
# not "how much scar protein is here" (CLAUDE.md step 12).
use <- intersect(fmk[in_fibrosis_set == FALSE, gene_symbol], rownames(mat))
cat("fibroblast markers available:", nrow(fmk), " not in a fibrosis set:",
    fmk[in_fibrosis_set == FALSE, .N], " of those on the array:", length(use), "\n")
# row z-score: for each marker gene, (value - gene mean) / gene sd, then average per sample.
# Worked example, 2 genes x 3 samples:
#   gene A: 8, 9, 10 -> mean 9, sd 1 -> z = -1, 0, +1
#   gene B: 5, 6, 7  -> mean 6, sd 1 -> z = -1, 0, +1
#   score  = -1, 0, +1.  So the score says "this sample is high/low on fibroblast genes
#   relative to the other samples", with every marker given equal weight.
z <- t(scale(t(mat[use, , drop = FALSE])))
fibro <- colMeans(z, na.rm = TRUE)
sheet[, fibro := fibro[gsm]]
cat("\n== fibroblast score by group ==\n")
print(sheet[, .(n = .N, mean = round(mean(fibro), 3), sd = round(sd(fibro), 3),
                min = round(min(fibro), 3), max = round(max(fibro), 3)), by = group])
cat("\nANOVA of the score on group: p =",
    signif(anova(lm(fibro ~ group, data = sheet))$`Pr(>F)`[1], 3), "\n")
cat("CAVEAT (CLAUDE.md step 12): scar load genuinely differs by group, so this covariate\n")
cat("absorbs real disease signal too. Read the result as 'does the footprint survive',\n")
cat("not as the main estimate.\n")
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
fwrite(sheet[, .(gsm, group, age_years, severity, fibro_score = fibro)],
       "results/tables/12_fibroblast_scores.tsv", sep = "\t")

## ---- refit with the covariate ---------------------------------------------
grp <- factor(sheet$group, levels = c("control", "BMD", "DMD"))
d <- model.matrix(~ 0 + grp + sheet$fibro)
colnames(d) <- c(levels(grp), "fibro")
fit <- eBayes(contrasts.fit(lmFit(mat, d),
       makeContrasts(BMDvsCtrl = BMD - control, DMDvsCtrl = DMD - control,
                     BMDvsDMD = BMD - DMD, levels = d)))
contrasts <- c("BMDvsCtrl", "DMDvsCtrl", "BMDvsDMD")
adj <- rbindlist(lapply(contrasts, function(cn) {
  tt <- topTable(fit, coef = cn, number = Inf, sort.by = "none")
  rk <- sort(setNames(tt$t, rownames(tt)), decreasing = TRUE)
  r <- as.data.table(fgsea(fs, rk, minSize = 15, maxSize = Inf, nPermSimple = 10000, eps = 0))
  r[, contrast := cn][]
}))
real <- readRDS("results/rds/fgsea_all.rds")
cmp <- merge(real[, .(contrast, family_id = pathway, size, NES_raw = NES, padj_raw = padj)],
             adj[,  .(contrast, family_id = pathway, NES_adj = NES, padj_adj = padj)],
             by = c("contrast", "family_id"))
cmp[, survives := padj_adj < 0.05 & sign(NES_adj) == sign(NES_raw)]
setorder(cmp, contrast, padj_raw)
fwrite(cmp, "results/tables/12_fibro_adjusted_fgsea.tsv", sep = "\t")

cat("\n=============== FIBROBLAST-ADJUSTED fgsea ===============\n")
print(cmp[, .(n = .N, sig_raw = sum(padj_raw < 0.05), sig_adj = sum(padj_adj < 0.05),
              survive = sum(survives), r_NES = round(cor(NES_raw, NES_adj), 3)), by = contrast])
cat("\n== miR-29, raw vs fibroblast-adjusted ==\n")
print(cmp[family_id == "miR-29-3p", .(contrast, NES_raw = round(NES_raw, 3),
          padj_raw = signif(padj_raw, 2), NES_adj = round(NES_adj, 3),
          padj_adj = signif(padj_adj, 2), survives)])
cat("\n== top 12 DMD families, raw vs adjusted ==\n")
print(head(cmp[contrast == "DMDvsCtrl", .(family_id, NES_raw = round(NES_raw, 2),
          NES_adj = round(NES_adj, 2), padj_adj = signif(padj_adj, 2), survives)], 12))
saveRDS(list(cmp = cmp, scores = sheet[, .(gsm, group, fibro)], markers_used = use),
        "results/rds/fibro_cov.rds")

p <- ggplot(sheet, aes(group, fibro, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = .5) + geom_jitter(width = .15, size = 1.4) +
  labs(x = NULL, y = "fibroblast score (mean row z of non-ECM fibroblast markers)",
       title = "How much fibroblast is in each biopsy",
       subtitle = sprintf("%d markers from the Rubenstein skeletal-muscle atlas that are NOT in a fibrosis set.\nHigher = more fibroblast relative to the other samples.", length(use))) +
  theme_bw(base_size = 9) + theme(legend.position = "none")
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("results/figures/12_fibroblast_score.png", p, width = 5.5, height = 4, dpi = 200)
cat("Wrote results/figures/12_fibroblast_score.png\n")
