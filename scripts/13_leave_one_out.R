#!/usr/bin/env Rscript
# 13_leave_one_out.R -- does the top-family list depend on any single sample?
# CLAUDE.md s6 step 13. Drops each BMD (11 runs) and each control (6 runs); also
# a dedicated run without GSM2934804 (the disputed-severity BMD sample).
suppressPackageStartupMessages({library(data.table); library(limma); library(fgsea)
                                library(BiocParallel); library(ggplot2)})
set.seed(20261105)
N_CORE <- as.integer(Sys.getenv("N_CORE", "7"))
N_TOP  <- as.integer(Sys.getenv("N_TOP", "20"))

fs  <- readRDS("results/rds/family_sets.rds")
gm  <- readRDS("data/processed/gene_matrix.rds")
mat <- gm$mat; sheet <- as.data.table(gm$sheet); setkey(sheet, gsm); sheet <- sheet[colnames(mat)]
real <- readRDS("results/rds/fgsea_all.rds")
top_real <- head(real[contrast == "BMDvsCtrl"][order(padj, -abs(NES)), pathway], N_TOP)

drop_set <- sheet[group %in% c("BMD", "control"), gsm]
cat("leave-one-out over", length(drop_set), "samples (", sheet[group=="BMD", .N], "BMD +",
    sheet[group=="control", .N], "control )\n")

run_one <- function(drop) {
  keep <- colnames(mat) != drop
  g <- factor(sheet$group[keep], levels = c("control", "BMD", "DMD"))
  d <- model.matrix(~ 0 + g); colnames(d) <- levels(g)
  f <- eBayes(contrasts.fit(lmFit(mat[, keep], d),
       makeContrasts(BMDvsCtrl = BMD - control, BMDvsDMD = BMD - DMD, levels = d)))
  rbindlist(lapply(c("BMDvsCtrl", "BMDvsDMD"), function(cn) {
    tt <- topTable(f, coef = cn, number = Inf, sort.by = "none")
    rk <- sort(setNames(tt$t, rownames(tt)), decreasing = TRUE)
    r <- as.data.table(fgsea(fs, rk, minSize = 15, maxSize = Inf, nPermSimple = 10000, eps = 0, nproc = 1))
    r[, absNES := abs(NES)]; setorder(r, padj, -absNES); r[, rk := .I]
    r[, .(dropped = drop, dropped_group = sheet$group[sheet$gsm == drop],
          contrast = cn, family_id = pathway, NES, padj, rank = rk)]
  }))
}

bp <- MulticoreParam(workers = N_CORE, RNGseed = 20261105)
t0 <- Sys.time()
loo <- rbindlist(bplapply(drop_set, run_one, BPPARAM = bp))
cat(sprintf("done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
fwrite(loo, "results/tables/13_leave_one_out.tsv", sep = "\t")

## ---- stability of the top-20 list -----------------------------------------
stab <- loo[contrast == "BMDvsCtrl", .(top = list(head(family_id[order(rank)], N_TOP))), by = dropped]
stab[, overlap := vapply(top, function(x) length(intersect(x, top_real)), integer(1))]
cat("\n=============== LEAVE-ONE-OUT (BMD vs control) ===============\n")
cat(sprintf("top-%d list overlap with the full-data list: min %d, median %.1f, max %d (of %d)\n",
            N_TOP, min(stab$overlap), median(stab$overlap), max(stab$overlap), N_TOP))
print(stab[order(overlap), .(dropped, overlap_with_full_top20 = overlap)])

## ---- miR-29 specifically ---------------------------------------------------
m <- loo[family_id == "miR-29-3p"]
cat("\n== miR-29 across every leave-one-out run ==\n")
for (cn in unique(m$contrast)) {
  d <- m[contrast == cn]
  cat(sprintf("%-10s NES %.3f to %.3f (full data %.3f); padj max %.2g; worst rank %d\n", cn,
              min(d$NES), max(d$NES), real[contrast == cn & pathway == "miR-29-3p", NES],
              max(d$padj), max(d$rank)))
  cat("   runs where miR-29 stays padj<0.05:", d[padj < 0.05, .N], "of", nrow(d), "\n")
}
fwrite(m, "results/tables/13_mir29_leave_one_out.tsv", sep = "\t")

## ---- the GSM2934804 question -----------------------------------------------
cat("\n== dropping GSM2934804 (disputed severity) ==\n")
print(m[dropped == "GSM2934804", .(contrast, NES = round(NES, 3), padj = signif(padj, 3), rank)])
cat("full data for comparison:\n")
print(real[pathway == "miR-29-3p" & contrast %in% c("BMDvsCtrl", "BMDvsDMD"),
           .(contrast, NES = round(NES, 3), padj = signif(padj, 3))])

p <- ggplot(m, aes(reorder(dropped, NES), NES, fill = dropped_group)) +
  geom_col() + coord_flip() + facet_wrap(~contrast, scales = "free_x") +
  geom_hline(data = data.table(contrast = c("BMDvsCtrl","BMDvsDMD"),
                               y = c(real[contrast=="BMDvsCtrl" & pathway=="miR-29-3p", NES],
                                     real[contrast=="BMDvsDMD" & pathway=="miR-29-3p", NES])),
             aes(yintercept = y), linetype = 2, colour = "#B2182B") +
  labs(x = "sample left out", y = "miR-29 NES",
       title = "miR-29 with each BMD or control sample removed in turn",
       subtitle = "Red dashed = the full-data value. Short bars = that sample was carrying the result.") +
  theme_bw(base_size = 8)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("results/figures/13_mir29_leave_one_out.png", p, width = 8, height = 5, dpi = 200)
cat("Wrote results/figures/13_mir29_leave_one_out.png\n")
saveRDS(list(loo = loo, stab = stab), "results/rds/leave_one_out.rds")
