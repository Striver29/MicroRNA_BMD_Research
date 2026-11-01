#!/usr/bin/env Rscript
# 08_shuffle_baseline.R -- the shared-control trap.  CLAUDE.md s6 step 8.
#
# Why this exists.  BMD-vs-control and DMD-vs-control are both measured against the SAME
# 6 controls. Any fluke in those 6 controls pushes BOTH comparisons the same way, so the
# two sets of family scores correlate even when nothing real is shared. Judging r against 0
# is therefore meaningless. We rebuild the whole pipeline on shuffled labels (keeping
# 6/11/17) and ask: how big is r when there is nothing real left to find?
suppressPackageStartupMessages({library(data.table); library(limma); library(fgsea)
                                library(BiocParallel); library(ggplot2)})
set.seed(20261101)

N_SHUFFLE  <- as.integer(Sys.getenv("N_SHUFFLE", "1000"))
NPERM_NULL <- as.integer(Sys.getenv("NPERM_NULL", "1000"))  # fgsea perms inside a shuffle
N_CORE     <- as.integer(Sys.getenv("N_CORE", "7"))
TIME_ONLY  <- as.integer(Sys.getenv("TIME_ONLY", "0"))

fs  <- readRDS("results/rds/family_sets.rds")
gm  <- readRDS("data/processed/gene_matrix.rds")
mat <- gm$mat; sheet <- as.data.table(gm$sheet)
setkey(sheet, gsm); sheet <- sheet[colnames(mat)]
grp0 <- sheet$group
stopifnot(identical(as.integer(table(grp0)[c("control","BMD","DMD")]), c(6L,11L,17L)))

agree <- readRDS("results/rds/bmd_dmd_agreement.rds")
r_real <- agree$r_real; sl_real <- agree$slope_real
cat(sprintf("REAL r = %.4f   REAL slope = %.4f  (across %d families)\n",
            r_real, sl_real, nrow(agree$xy)))

## one draw: relabel, refit limma, rerun fgsea for both control contrasts, return r + slope
one <- function(i) {
  set.seed(100000L + i)
  g <- factor(sample(grp0), levels = c("control", "BMD", "DMD"))
  d <- model.matrix(~ 0 + g); colnames(d) <- levels(g)
  f <- eBayes(contrasts.fit(lmFit(mat, d),
        makeContrasts(BMDvsCtrl = BMD - control, DMDvsCtrl = DMD - control, levels = d)))
  nes <- lapply(c("BMDvsCtrl", "DMDvsCtrl"), function(cn) {
    tt <- topTable(f, coef = cn, number = Inf, sort.by = "none")
    rk <- sort(setNames(tt$t, rownames(tt)), decreasing = TRUE)
    r <- fgsea(fs, rk, minSize = 15, maxSize = Inf, nPermSimple = NPERM_NULL, eps = 0)
    setNames(r$NES, r$pathway)
  })
  fams <- intersect(names(nes[[1]]), names(nes[[2]]))
  x <- nes[[2]][fams]; y <- nes[[1]][fams]          # x = DMD, y = BMD, as in step 7
  ok <- is.finite(x) & is.finite(y)
  list(r = cor(x[ok], y[ok]),
       slope = as.numeric(coef(lm(y[ok] ~ x[ok]))[2]),
       nes_bmd = nes[[1]], nes_dmd = nes[[2]])
}

## ---- timing first (CLAUDE.md: time 5 iterations and extrapolate) ----------
t0 <- Sys.time(); invisible(lapply(1:5, one)); el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
per <- el / 5
cat(sprintf("\nTiming: %.1f s for 5 iterations = %.2f s each (serial).\n", el, per))
cat(sprintf("Projected for %d shuffles on %d cores: %.1f min\n",
            N_SHUFFLE, N_CORE, N_SHUFFLE * per / N_CORE / 60))
if (TIME_ONLY == 1L) { cat("TIME_ONLY=1, stopping here.\n"); quit(status = 0) }

bp <- MulticoreParam(workers = N_CORE, RNGseed = 20261101, progressbar = FALSE)
cat(sprintf("\nRunning %d shuffles on %d cores ...\n", N_SHUFFLE, N_CORE))
t1 <- Sys.time()
out <- bplapply(seq_len(N_SHUFFLE), one, BPPARAM = bp)
cat(sprintf("Done in %.1f min\n", as.numeric(difftime(Sys.time(), t1, units = "mins"))))

null_r  <- vapply(out, `[[`, numeric(1), "r")
null_sl <- vapply(out, `[[`, numeric(1), "slope")
null_r  <- null_r[is.finite(null_r)]; null_sl <- null_sl[is.finite(null_sl)]

p_r  <- (sum(null_r  >= r_real)  + 1) / (length(null_r)  + 1)
p_sl <- (sum(null_sl >= sl_real) + 1) / (length(null_sl) + 1)

cat("\n=================== SHUFFLED-LABEL BASELINE ===================\n")
cat(sprintf("shuffles that finished : %d\n", length(null_r)))
cat(sprintf("null r   : mean %.3f  sd %.3f  median %.3f  95th pct %.3f  max %.3f\n",
            mean(null_r), sd(null_r), median(null_r), quantile(null_r, .95), max(null_r)))
cat(sprintf("REAL r   : %.4f    empirical p = %.4g\n", r_real, p_r))
cat(sprintf("null slope: mean %.3f  sd %.3f  median %.3f  95th pct %.3f  max %.3f\n",
            mean(null_sl), sd(null_sl), median(null_sl), quantile(null_sl, .95), max(null_sl)))
cat(sprintf("REAL slope: %.4f   empirical p = %.4g\n", sl_real, p_sl))

## ---- free by-product: per-family label-shuffle p --------------------------
nb <- do.call(rbind, lapply(out, function(z) z$nes_bmd[names(fs)]))
nd <- do.call(rbind, lapply(out, function(z) z$nes_dmd[names(fs)]))
real <- readRDS("results/rds/fgsea_all.rds")
fam_p <- rbindlist(lapply(names(fs), function(f) {
  rb <- real[contrast == "BMDvsCtrl" & pathway == f, NES]
  rd <- real[contrast == "DMDvsCtrl" & pathway == f, NES]
  cb <- nb[, f]; cd <- nd[, f]; cb <- cb[is.finite(cb)]; cd <- cd[is.finite(cd)]
  data.table(family_id = f,
             NES_BMDvsCtrl = rb, shuffle_p_BMDvsCtrl = (sum(abs(cb) >= abs(rb)) + 1) / (length(cb) + 1),
             NES_DMDvsCtrl = rd, shuffle_p_DMDvsCtrl = (sum(abs(cd) >= abs(rd)) + 1) / (length(cd) + 1))
}))
setorder(fam_p, shuffle_p_DMDvsCtrl, -NES_DMDvsCtrl)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
fwrite(fam_p, "results/tables/08_family_shuffle_p.tsv", sep = "\t")
cat("\nWrote results/tables/08_family_shuffle_p.tsv\n")
cat("\n== label-shuffle p for miR-29 ==\n"); print(fam_p[family_id == "miR-29-3p"])
cat("\nfamilies with label-shuffle p < 0.05 (DMD vs ctrl):",
    fam_p[shuffle_p_DMDvsCtrl < 0.05, .N], "of", nrow(fam_p), "\n")
cat("families with label-shuffle p < 0.05 (BMD vs ctrl):",
    fam_p[shuffle_p_BMDvsCtrl < 0.05, .N], "of", nrow(fam_p), "\n")

saveRDS(list(null_r = null_r, null_slope = null_sl, r_real = r_real, slope_real = sl_real,
             p_r = p_r, p_slope = p_sl, n_shuffle = length(null_r),
             nperm_null = NPERM_NULL, fam_p = fam_p),
        "results/rds/shuffle_baseline.rds")

## ---- the two histograms ---------------------------------------------------
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
mk <- function(v, real, lab, p) {
  ggplot(data.table(v = v), aes(v)) +
    geom_histogram(bins = 40, fill = "grey75", colour = "white") +
    geom_vline(xintercept = real, colour = "#B2182B", linewidth = 1) +
    annotate("label", x = real, y = Inf, vjust = 1.3, size = 2.8, colour = "#B2182B",
             label = sprintf("real = %.3f\nempirical p = %.3g", real, p)) +
    labs(x = lab, y = sprintf("number of shuffles (of %d)", length(v)),
         title = sprintf("Shuffled-label baseline for %s", lab),
         subtitle = "Grey = labels shuffled at random keeping 6/11/17, so nothing real is left.\nRed line = the value we actually measured.") +
    theme_bw(base_size = 9)
}
ggsave("results/figures/08_shuffle_r.png",
       mk(null_r, r_real, "r between BMD and DMD family NES", p_r), width = 6.5, height = 4.2, dpi = 200)
ggsave("results/figures/08_shuffle_slope.png",
       mk(null_sl, sl_real, "slope of BMD NES on DMD NES", p_sl), width = 6.5, height = 4.2, dpi = 200)
cat("Wrote results/figures/08_shuffle_r.png and 08_shuffle_slope.png\n")
