#!/usr/bin/env Rscript
# 02_rma_qc.R -- RMA on the 34 included CELs + QC without scatter plots.  CLAUDE.md s6 step 2.
suppressPackageStartupMessages({library(data.table); library(oligo); library(ggplot2)})
set.seed(1)

sheet_path <- if (file.exists("inputs/sample_sheet.csv")) "inputs/sample_sheet.csv" else "inputs/auto/sample_sheet.csv"
cat("Using sample sheet:", sheet_path, "\n")
sheet <- fread(sheet_path)
inc <- sheet[include == TRUE]
stopifnot(nrow(inc) == 34L)
cat("Group sizes:\n"); print(table(inc$group))

cel <- file.path("data/raw/CEL", inc$cel_file)
stopifnot(all(file.exists(cel)))

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
if (file.exists("data/processed/rma.rds")) {
  cat("\nReusing data/processed/rma.rds (delete it to force a re-run of RMA).\n")
  mat <- readRDS("data/processed/rma.rds")$mat
} else {
  cat("\nReading", length(cel), "CEL files ...\n")
  raw <- read.celfiles(cel, verbose = FALSE)
  cat("Platform:", annotation(raw), "\n")
  cat("Running RMA (background correct -> quantile normalise -> median polish) ...\n")
  es <- rma(raw)
  mat <- exprs(es)
  colnames(mat) <- inc$gsm                    # CEL order == sheet order by construction
  saveRDS(list(mat = mat, sheet = inc), "data/processed/rma.rds")
  cat("Wrote data/processed/rma.rds\n")
}
cat("RMA matrix:", nrow(mat), "x", ncol(mat), "\n")
stopifnot(nrow(mat) == 54675L, ncol(mat) == 34L)

## ---- QC 1: boxplot --------------------------------------------------------
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dl <- as.data.table(mat)[, melt(.SD, measure.vars = names(.SD),
                                variable.name = "gsm", value.name = "expr")]
dl <- merge(dl, inc[, .(gsm, group)], by = "gsm", sort = FALSE)
dl[, gsm := factor(gsm, levels = inc[order(group, gsm), gsm])]
p <- ggplot(dl, aes(gsm, expr, fill = group)) +
  geom_boxplot(outlier.size = 0.2, linewidth = 0.25) +
  labs(x = NULL, y = expression(log[2]~"RMA expression"),
       title = "QC: per-sample log2 expression after RMA (n = 34)",
       subtitle = "Boxes should line up. A shifted or squashed box = a sample to question.") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave("results/figures/02_qc_boxplot.png", p, width = 9, height = 4.5, dpi = 200)
cat("Wrote results/figures/02_qc_boxplot.png\n")

## ---- QC 2: inter-sample correlation, OVERALL and WITHIN GROUP -------------
# Why both: a sample's correlation with *all* others is driven by how many samples
# share its biology. Here 28 of 34 arrays are dystrophic muscle and only 6 are healthy,
# so every control looks like an outlier on the overall metric even if its array is fine.
# The within-group metric (each sample vs the others in its own group) removes that.
cm <- cor(mat, method = "pearson"); diag(cm) <- NA_real_
grp <- setNames(inc$group, inc$gsm)
within_med <- vapply(colnames(cm), function(g)
  median(cm[names(grp)[grp == grp[[g]] & names(grp) != g], g], na.rm = TRUE), numeric(1))
qc <- data.table(gsm = colnames(cm), group = as.character(grp[colnames(cm)]),
                 median_cor_all = apply(cm, 2, median, na.rm = TRUE),
                 median_cor_within_group = as.numeric(within_med))
# flag rule: >3 MADs below the centre of the WITHIN-GROUP metric, computed per group
qc[, `:=`(grp_centre = median(median_cor_within_group),
           grp_mad    = mad(median_cor_within_group)), by = group]
qc[, mad_below_within := fifelse(grp_mad > 0, (grp_centre - median_cor_within_group) / grp_mad, 0)]
qc[, flagged := mad_below_within > 3]
setorder(qc, median_cor_within_group)

## ---- QC 3: RLE and NUSE from the probe-level model ------------------------
# These are the standard Affymetrix array-quality measures and they are computed
# within-probeset, so group biology does not drive them.
#   RLE  = Relative Log Expression. Each probeset's value minus the median across arrays.
#          A healthy array has a median RLE near 0 and a narrow spread.
#   NUSE = Normalised Unscaled Standard Error. Scaled so the median array sits at 1.
#          A NUSE median above ~1.05 means that array's probe fits are unusually noisy.
plm_file <- "data/processed/plm_qc.rds"
if (file.exists(plm_file)) {
  pq_stats <- readRDS(plm_file)
} else {
  cat("\nFitting the probe-level model for RLE/NUSE (slow, cached afterwards) ...\n")
  raw2 <- read.celfiles(cel, verbose = FALSE)
  plm <- fitProbeLevelModel(raw2)
  rle_m <- oligo::RLE(plm, type = "values"); nuse_m <- oligo::NUSE(plm, type = "values")
  pq_stats <- data.table(gsm = inc$gsm,
                         rle_median = apply(rle_m, 2, median, na.rm = TRUE),
                         rle_iqr    = apply(rle_m, 2, IQR, na.rm = TRUE),
                         nuse_median = apply(nuse_m, 2, median, na.rm = TRUE))
  saveRDS(pq_stats, plm_file)
}
qc <- merge(qc, pq_stats, by = "gsm", sort = FALSE)
# conventional thresholds
qc[, rle_bad  := abs(rle_median) > 0.1]
qc[, nuse_bad := nuse_median > 1.05]
setorder(qc, -nuse_median)

fwrite(qc, "results/tables/02_qc_metrics.tsv", sep = "\t")
cat("\n== QC table (sorted by NUSE, worst first) ==\n")
print(qc[, .(gsm, group, cor_all = round(median_cor_all, 4),
             cor_within = round(median_cor_within_group, 4),
             rle_med = round(rle_median, 3), rle_iqr = round(rle_iqr, 3),
             nuse_med = round(nuse_median, 3), flagged, rle_bad, nuse_bad)])

cat("\n== overall-correlation metric, by group (shows the group effect) ==\n")
print(qc[, .(n = .N, median_cor_all = round(median(median_cor_all), 4),
             median_cor_within = round(median(median_cor_within_group), 4)), by = group])

bad <- qc[flagged | rle_bad | nuse_bad, gsm]
cat("\nFLAGGED on the within-group correlation:",
    if (any(qc$flagged)) paste(qc[flagged == TRUE, gsm], collapse = ", ") else "none", "\n")
cat("FLAGGED on |median RLE| > 0.10      :",
    if (any(qc$rle_bad)) paste(qc[rle_bad == TRUE, gsm], collapse = ", ") else "none", "\n")
cat("FLAGGED on median NUSE > 1.05       :",
    if (any(qc$nuse_bad)) paste(qc[nuse_bad == TRUE, gsm], collapse = ", ") else "none", "\n")

pq <- ggplot(qc, aes(reorder(gsm, nuse_median), nuse_median, fill = group)) +
  geom_col() + coord_flip() +
  geom_hline(yintercept = 1.05, linetype = 2, colour = "red") +
  labs(x = NULL, y = "median NUSE (1.00 = the typical array here)",
       title = "QC: per-array noise (NUSE)",
       subtitle = "Red dashed line = 1.05, the usual flag threshold") +
  theme_bw(base_size = 8)
ggsave("results/figures/02_qc_nuse.png", pq, width = 6, height = 6, dpi = 200)

pc <- ggplot(qc, aes(reorder(gsm, median_cor_within_group), median_cor_within_group, fill = group)) +
  geom_col() + coord_flip(ylim = c(min(qc$median_cor_within_group) - 0.01, 1)) +
  labs(x = NULL, y = "median Pearson r with the other samples in its own group",
       title = "QC: array-to-array agreement within group",
       subtitle = "Within group, so the control-vs-patient biology difference cannot create a false outlier") +
  theme_bw(base_size = 8)
ggsave("results/figures/02_qc_within_group_cor.png", pc, width = 6, height = 6, dpi = 200)
cat("Wrote results/figures/02_qc_nuse.png and 02_qc_within_group_cor.png\n")

qc[, flagged := gsm %in% bad]
if (any(qc$flagged)) cat("\n*** STOP (CLAUDE.md s6 step 2): Aziz decides keep or drop. ***\n")
