#!/usr/bin/env Rscript
# 10b_detection_limit.R -- how big a footprint would we need to catch at n = 11?
# CLAUDE.md s6 step 10b.
#
# The idea in one line: we PLANT a fake microRNA footprint of known size in the 11 BMD
# samples, then see how often our own pipeline finds it.
#
# Worked example. Suppose the DMD miR-29 footprint is +0.27 log2 units (the "ruler").
# At f = 50% we add 0.5 * 0.27 = +0.135 to each of ~1158 planted genes, in the 11 BMD
# samples only. Then we rerun BMD-vs-control limma and fgsea and ask whether the planted
# set comes out with NES > 0 and padj < 0.05. Repeat 100 times to get a catch rate.
# At f = 0 nothing is added, so the catch rate there IS the false-alarm rate.
suppressPackageStartupMessages({library(data.table); library(limma); library(fgsea)
                                library(BiocParallel); library(ggplot2)})
set.seed(20261104)
N_REP   <- as.integer(Sys.getenv("N_REP", "100"))
N_CORE  <- as.integer(Sys.getenv("N_CORE", "7"))
NPERM   <- as.integer(Sys.getenv("NPERM", "10000"))
GRID    <- as.numeric(strsplit(Sys.getenv("GRID", "0,10,25,50,75,100"), ",")[[1]]) / 100
TIME_ONLY <- as.integer(Sys.getenv("TIME_ONLY", "0"))
pick <- function(f) { a <- file.path("inputs", f); if (file.exists(a)) a else file.path("inputs/auto", f) }

fs  <- readRDS("results/rds/family_sets.rds")
gm  <- readRDS("data/processed/gene_matrix.rds")
mat <- gm$mat; sheet <- as.data.table(gm$sheet); setkey(sheet, gsm); sheet <- sheet[colnames(mat)]
grp <- factor(sheet$group, levels = c("control", "BMD", "DMD"))
is_bmd <- which(grp == "BMD")
utr <- fread(pick("utr_lengths.tsv"))

## ---- 1. the ruler: the size of the DMD miR-29 footprint -------------------
tt_dmd <- readRDS("results/rds/limma.rds")$tables$DMDvsCtrl
m29 <- fs[["miR-29-3p"]]
RULER <- tt_dmd[gene_symbol %in% m29, mean(logFC)] - tt_dmd[!gene_symbol %in% m29, mean(logFC)]
cat(sprintf("RULER = mean DMD-vs-ctrl log2FC of %d miR-29 targets minus the other %d genes\n",
            length(m29), nrow(tt_dmd) - length(m29)))
cat(sprintf("      = %+.4f log2 units  (= %.1f%% change on the linear scale)\n",
            RULER, 100 * (2^RULER - 1)))

## ---- 2. UTR-length-matched planting pool ----------------------------------
u <- utr[gene_symbol %in% rownames(mat), .(gene_symbol, utr3_length_nt)]
u[, decile := cut(utr3_length_nt, quantile(utr3_length_nt, seq(0, 1, .1)),
                  include.lowest = TRUE, labels = FALSE)]
by_dec <- split(u$gene_symbol, u$decile)
m29u <- intersect(m29, u$gene_symbol)
prof <- table(factor(u[gene_symbol %in% m29u, decile], levels = 1:10))
cat(sprintf("planted sets copy miR-29's 3'UTR decile profile; set size = %d genes\n", length(m29u)))

one <- function(f, i) {
  set.seed(400000L + round(f * 1000) * 1000L + i)
  planted <- unlist(lapply(names(prof), function(d) {
    k <- prof[[d]]; if (k == 0) character(0) else sample(by_dec[[d]], k) }), use.names = FALSE)
  m <- mat
  m[planted, is_bmd] <- m[planted, is_bmd] + f * RULER     # the plant, BMD samples only
  d <- model.matrix(~ 0 + grp); colnames(d) <- levels(grp)
  fit <- eBayes(contrasts.fit(lmFit(m, d),
         makeContrasts(BMDvsCtrl = BMD - control, levels = d)))
  tt <- topTable(fit, coef = "BMDvsCtrl", number = Inf, sort.by = "none")
  rk <- sort(setNames(tt$t, rownames(tt)), decreasing = TRUE)
  sets <- c(fs, list(PLANTED = planted))        # tested among all families, as the spec says
  r <- as.data.table(fgsea(sets, rk, minSize = 15, maxSize = Inf, nPermSimple = NPERM, eps = 0, nproc = 1))
  p <- r[pathway == "PLANTED"]
  data.table(f = f, rep = i, NES = p$NES, pval = p$pval, padj = p$padj,
             caught = isTRUE(p$NES > 0 && p$padj < 0.05))
}

t0 <- Sys.time(); invisible(lapply(1:3, function(i) one(0.5, i)))
per <- as.numeric(difftime(Sys.time(), t0, units = "secs")) / 3
cat(sprintf("\nTiming: %.1f s per run; projected %.1f min for %d x %d on %d cores\n",
            per, length(GRID) * N_REP * per / N_CORE / 60, length(GRID), N_REP, N_CORE))
if (TIME_ONLY == 1L) quit(status = 0)

jobs <- CJ(f = GRID, rep = seq_len(N_REP))
bp <- MulticoreParam(workers = N_CORE, RNGseed = 20261104)
cat(sprintf("Running %d simulations ...\n", nrow(jobs)))
t1 <- Sys.time()
res <- rbindlist(bplapply(seq_len(nrow(jobs)),
        function(k) one(jobs$f[k], jobs$rep[k]), BPPARAM = bp))
cat(sprintf("Done in %.1f min\n", as.numeric(difftime(Sys.time(), t1, units = "mins"))))

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
fwrite(res, "results/tables/10b_detection_runs.tsv", sep = "\t")
curve <- res[, .(n = .N, catch_rate = mean(caught),
                 mean_NES = mean(NES, na.rm = TRUE),
                 planted_log2 = f[1] * RULER,
                 planted_pct_linear = 100 * (2^(f[1] * RULER) - 1)), by = f]
setorder(curve, f)
fwrite(curve, "results/tables/10b_detection_curve.tsv", sep = "\t")

cat("\n=============== DETECTION LIMIT AT n = 11 BMD ===============\n")
print(curve[, .(f_pct = 100 * f, planted_log2 = round(planted_log2, 4),
                planted_pct_linear = round(planted_pct_linear, 1),
                mean_NES = round(mean_NES, 2), catch_rate = catch_rate, n)])
cat(sprintf("\nFALSE-ALARM RATE at f = 0: %.3f (we want this near 0.05 or below)\n",
            curve[f == 0, catch_rate]))
hit <- curve[catch_rate >= 0.8]
cat(if (nrow(hit)) sprintf("80%% detection reached at f = %.0f%% of the DMD miR-29 footprint\n",
                           100 * min(hit$f))
    else "80% detection NOT reached anywhere on the grid -- the limit is above 100% of the DMD footprint\n")

p <- ggplot(curve, aes(100 * f, catch_rate)) +
  geom_hline(yintercept = 0.8, linetype = 2, colour = "#B2182B") +
  geom_hline(yintercept = 0.05, linetype = 3, colour = "grey40") +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(x = "planted footprint, as a % of the DMD miR-29 footprint",
       y = "share of simulations where we caught it",
       title = sprintf("Detection limit for a microRNA footprint at n = 11 BMD (%d reps per point)", N_REP),
       subtitle = sprintf("Catch = NES > 0 and padj < 0.05. Ruler = %+.3f log2 units.\nRed dashed = 80%% power. Dotted = 5%%, the false-alarm level.", RULER)) +
  theme_bw(base_size = 9)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("results/figures/10b_detection_curve.png", p, width = 6.5, height = 4.3, dpi = 200)
cat("Wrote results/figures/10b_detection_curve.png\n")
saveRDS(list(curve = curve, runs = res, ruler = RULER), "results/rds/detection_limit.rds")
