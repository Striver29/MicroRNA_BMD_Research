#!/usr/bin/env Rscript
# 10b2_calibrated_catch.R -- re-score step 10b's 600 saved simulations with a catch rule
# whose false-alarm rate is actually 5%. No limma or fgsea rerun; reads 10b's saved runs.
#
# Why this exists: with the spec's rule (NES > 0 and padj < 0.05) step 10b caught the
# planted set in 100/100 runs at f = 0, when nothing was planted. A random gene set copying
# miR-29's 3'UTR-length profile already scores NES ~ +1.68 (step 11 found the same, 1.67),
# and fgsea's gene-permutation padj calls that significant every time.
#
# Worked example. The 100 unplanted (f = 0) runs give NES values from 1.46 to 1.87; their
# 95th percentile is the bar, ~1.81. At f = 10% a run with NES 1.87 clears the bar (caught);
# a run with NES 1.75 does not. At f = 0 exactly 5 of 100 clear it -- by construction.
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
set.seed(20261105)
pick <- function(f) { a <- file.path("inputs", f); if (file.exists(a)) a else file.path("inputs/auto", f) }

dl  <- readRDS("results/rds/detection_limit.rds")
res <- as.data.table(dl$runs); RULER <- dl$ruler

## ---- 1. calibrated catch: beat the 95th percentile of the unplanted NES ------
null <- res[f == 0, NES]
BAR  <- unname(quantile(null, 0.95))
cat(sprintf("Unplanted (f = 0) NES over %d runs: mean %.3f, SD %.3f, range %.3f to %.3f\n",
            length(null), mean(null), sd(null), min(null), max(null)))
cat(sprintf("Calibrated bar = 95th percentile = %.3f\n", BAR))
cat(sprintf("Spec rule at f = 0: %d/%d caught; the LEAST significant unplanted set had padj %.1e\n",
            res[f == 0, sum(caught)], length(null), res[f == 0, max(padj)]))

curve <- res[, .(n = .N, planted_log2 = f[1] * RULER,
                 median_NES = median(NES), z_vs_unplanted = (mean(NES) - mean(null)) / sd(null),
                 catch_spec = mean(caught), catch_calibrated = mean(NES > BAR)), by = f]
setorder(curve, f)

# where the calibrated curve crosses 80%: straight-line interpolation between grid points
k <- which(curve$catch_calibrated >= 0.8)[1]
CROSS <- if (is.na(k)) NA_real_ else if (k == 1) curve$f[1] else
  with(curve, f[k - 1] + (0.8 - catch_calibrated[k - 1]) /
                         (catch_calibrated[k] - catch_calibrated[k - 1]) * (f[k] - f[k - 1]))

## ---- 2. how much of the ruler is 3'UTR length alone? -------------------------
# The ruler (+0.2686) is miR-29 targets minus ALL other genes, so it includes the drift that
# any long-UTR set shows. Same formula on UTR-matched random sets tells us that share.
fs  <- readRDS("results/rds/family_sets.rds"); m29 <- fs[["miR-29-3p"]]
tt  <- readRDS("results/rds/limma.rds")$tables
utr <- fread(pick("utr_lengths.tsv"))
u <- utr[gene_symbol %in% tt$DMDvsCtrl$gene_symbol, .(gene_symbol, utr3_length_nt)]
u[, decile := cut(utr3_length_nt, quantile(utr3_length_nt, seq(0, 1, .1)),
                  include.lowest = TRUE, labels = FALSE)]
by_dec <- split(u$gene_symbol, u$decile)
prof <- table(factor(u[gene_symbol %in% m29, decile], levels = 1:10))
gap <- function(t, s) t[gene_symbol %in% s, mean(logFC)] - t[!gene_symbol %in% s, mean(logFC)]
draw <- function() unlist(lapply(names(prof), function(d)
  if (prof[[d]] == 0) character(0) else sample(by_dec[[d]], prof[[d]])), use.names = FALSE)
utr_gap <- rbindlist(lapply(c("DMDvsCtrl", "BMDvsCtrl"), function(cn) {
  g <- replicate(1000, gap(tt[[cn]], draw()))
  data.table(contrast = cn, mir29_gap = gap(tt[[cn]], m29), utr_matched_gap_mean = mean(g),
             utr_matched_gap_sd = sd(g)) }))
utr_gap[, utr_share_of_mir29_gap := utr_matched_gap_mean / mir29_gap]
utr_gap[, mir29_gap_pct_of_ruler := 100 * mir29_gap / RULER]

## ---- 3. outputs ---------------------------------------------------------------
fwrite(curve, "results/tables/10b2_detection_curve_calibrated.tsv", sep = "\t")
fwrite(utr_gap, "results/tables/10b2_ruler_utr_share.tsv", sep = "\t")

cat("\n========== DETECTION LIMIT AT n = 11 BMD, TWO CATCH RULES ==========\n")
print(curve[, .(f_pct = 100 * f, planted_log2 = round(planted_log2, 4), median_NES = round(median_NES, 2),
                z = round(z_vs_unplanted, 2), catch_spec, catch_calibrated, n)])
cat(sprintf("\nSpec rule (NES > 0, padj < 0.05): false-alarm rate %.2f -> curve uninformative\n",
            curve[f == 0, catch_spec]))
cat(if (is.na(CROSS)) "Calibrated rule: 80% NOT reached on the grid\n" else
    sprintf("Calibrated rule: 80%% crossed at f ~ %.0f%% of the DMD footprint (%+.3f log2), by interpolation\n",
            100 * CROSS, CROSS * RULER))
cat("\nShare of the miR-29 ruler that UTR-matched random sets already show (1,000 draws):\n")
print(utr_gap[, .(contrast, mir29_gap = round(mir29_gap, 4), utr_matched_gap = round(utr_matched_gap_mean, 4),
                  utr_share = round(utr_share_of_mir29_gap, 3), pct_of_ruler = round(mir29_gap_pct_of_ruler, 1))])

long <- melt(curve, id.vars = "f", measure.vars = c("catch_spec", "catch_calibrated"),
             variable.name = "rule", value.name = "catch")
long[, rule := factor(rule, levels = c("catch_spec", "catch_calibrated"),
       labels = c("spec rule: NES > 0 and padj < 0.05",
                  sprintf("calibrated: NES > %.2f (95th pct of unplanted runs)", BAR)))]
p <- ggplot(long, aes(100 * f, catch, colour = rule)) +
  geom_hline(yintercept = 0.8, linetype = 2, colour = "#B2182B") +
  geom_hline(yintercept = 0.05, linetype = 3, colour = "grey40") +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  scale_colour_manual(values = c("grey55", "#2166AC"), name = NULL) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(x = "planted footprint, as a % of the DMD miR-29 footprint",
       y = "share of 100 simulations where we caught it",
       title = "Detection limit at n = 11 BMD, under two catch rules",
       subtitle = sprintf("Ruler = %+.3f log2. Red dashed = 80%%. Dotted = 5%% false-alarm level.", RULER)) +
  theme_bw(base_size = 9) + theme(legend.position = "bottom", legend.direction = "vertical")
ggsave("results/figures/10b2_detection_curve_calibrated.png", p, width = 6.5, height = 4.8, dpi = 200)
cat("Wrote results/figures/10b2_detection_curve_calibrated.png\n")
saveRDS(list(curve = curve, bar = BAR, cross_f = CROSS, utr_gap = utr_gap),
        "results/rds/detection_limit_calibrated.rds")
