#!/usr/bin/env Rscript
# 09b_classify.R -- completes CLAUDE.md s6 step 9: the shared / BMD-only / DMD-only table.
#
# Why this script exists. 09_bmd_vs_dmd.R ran on 3 Nov at 01:49, before step 8 had finished
# (08_family_shuffle_p.tsv was not written until 05:46). Its final block needs that file, so
# the run halted there with a fread() error and no classification table was ever produced.
# The expensive part of step 9 -- 500 patient-label shuffles -- completed fine and its results
# are in 09_bmd_vs_dmd.tsv, so nothing needs recomputing. This reads the saved tables only.
#
# What "classify" means. Each microRNA family gets put in one box:
#   shared   = raised against controls in BOTH groups, and the direct BMD-vs-DMD test
#              finds no difference between them.
#   BMD-only = the direct test says the two groups differ, and only BMD is raised.
#   DMD-only = the direct test says they differ, and only DMD is raised.
#   neither  = not raised in either group (or differs but neither is individually raised).
#
# Worked example with made-up numbers. Family X: padj 0.01 vs controls in BMD, padj 0.02 in
# DMD, and the direct BMD-vs-DMD shuffle gives q = 0.60. Both raised, no difference -> "shared".
# Family Y: padj 0.80 in BMD, padj 0.01 in DMD, direct q = 0.01 -> differs, only DMD -> "DMD-only".
#
# TWO classifications are produced, because CLAUDE.md s10 requires both yardsticks:
#   (a) the spec one, using fgsea's padj to decide "raised vs controls";
#   (b) an honest one, using step 8's patient-label shuffle p instead.
# (b) matters because this project's central finding is that fgsea's padj is overconfident.
suppressPackageStartupMessages(library(data.table))

need <- c("results/tables/09_bmd_vs_dmd.tsv", "results/tables/07_fgsea_wide.tsv",
          "results/tables/08_family_shuffle_p.tsv")
miss <- need[!file.exists(need)]
if (length(miss)) stop("missing input(s): ", paste(miss, collapse = ", "))

res <- fread("results/tables/09_bmd_vs_dmd.tsv")
w   <- fread("results/tables/07_fgsea_wide.tsv")
s8  <- fread("results/tables/08_family_shuffle_p.tsv")

cls <- merge(w[, .(family_id, size, NES_BMDvsCtrl, padj_BMDvsCtrl, NES_DMDvsCtrl, padj_DMDvsCtrl)],
             res[, .(family_id, NES_BMDvsDMD, diff_p = shuffle_p, diff_q = shuffle_q)],
             by = "family_id")
cls <- merge(cls, s8[, .(family_id, shuffle_p_BMDvsCtrl, shuffle_p_DMDvsCtrl)], by = "family_id")

# "differs" is identical in both schemes: it comes from the direct comparison (step 9), where
# the 6 shared controls play no part at all.
cls[, differs := diff_q < 0.05]

box <- function(sig_B, sig_D, differs)
  fifelse(!differs & sig_B & sig_D, "shared",
  fifelse( differs & sig_D & !sig_B, "DMD-only",
  fifelse( differs & sig_B & !sig_D, "BMD-only",
  fifelse( differs & sig_B &  sig_D, "both, but different size", "neither"))))

## (a) spec classification -- fgsea padj decides "raised vs controls"
cls[, sig_BMD_fgsea := padj_BMDvsCtrl < 0.05]
cls[, sig_DMD_fgsea := padj_DMDvsCtrl < 0.05]
cls[, class_fgsea   := box(sig_BMD_fgsea, sig_DMD_fgsea, differs)]

## (b) honest classification -- step 8's patient-label shuffle decides it
cls[, sig_BMD_shuffle := shuffle_p_BMDvsCtrl < 0.05]
cls[, sig_DMD_shuffle := shuffle_p_DMDvsCtrl < 0.05]
cls[, class_shuffle   := box(sig_BMD_shuffle, sig_DMD_shuffle, differs)]

setcolorder(cls, c("family_id", "size", "NES_BMDvsCtrl", "NES_DMDvsCtrl", "NES_BMDvsDMD",
                   "padj_BMDvsCtrl", "padj_DMDvsCtrl", "shuffle_p_BMDvsCtrl",
                   "shuffle_p_DMDvsCtrl", "diff_p", "diff_q", "differs",
                   "class_fgsea", "class_shuffle"))
fwrite(cls, "results/tables/09_family_classification.tsv", sep = "\t")

cat("families classified:", nrow(cls), "\n")
cat("\n== (a) SPEC classification: 'raised vs controls' = fgsea padj < 0.05 ==\n")
print(cls[, .N, by = class_fgsea][order(-N)])
cat("\n== (b) HONEST classification: 'raised vs controls' = step 8 patient-shuffle p < 0.05 ==\n")
print(cls[, .N, by = class_shuffle][order(-N)])

cat("\n== families the DIRECT comparison says differ between BMD and DMD (q < 0.05) ==\n")
if (cls[differs == TRUE, .N] == 0)
  cat("NONE. 0 of", nrow(cls), "families. So no family can be BMD-only or DMD-only:\n",
      "every family is either 'shared' or 'neither'. This is the null result.\n", sep = " ") else
  print(cls[differs == TRUE][order(diff_q)])

cat("\n== the families called 'shared' under each rule ==\n")
cat("(a) fgsea:  ", cls[class_fgsea   == "shared", .N], "families\n")
cat("(b) shuffle:", cls[class_shuffle == "shared", .N], "families ->",
    paste(cls[class_shuffle == "shared", family_id], collapse = ", "), "\n")

cat("\n== miR-29 ==\n")
print(t(cls[family_id == "miR-29-3p"]))
cat("\nWrote results/tables/09_family_classification.tsv\n")
