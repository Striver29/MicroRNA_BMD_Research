#!/usr/bin/env Rscript
# 06_smoke_test.R -- the positive control: does miR-29 come out in DMD vs control?
# CLAUDE.md s6 step 6.  GO = NES > 0 and padj < 0.05.
suppressPackageStartupMessages({library(data.table); library(fgsea); library(ggplot2)})
set.seed(42)

fs <- readRDS("results/rds/family_sets.rds")
lim <- readRDS("results/rds/limma.rds")
tt  <- lim$tables$DMDvsCtrl

# the ranking vector: one number per gene, the limma moderated t for DMD - control.
# Large positive t = higher in DMD. fgsea walks this list from the top down.
ranks <- setNames(tt$t, tt$gene_symbol)
ranks <- ranks[!is.na(ranks)]
ranks <- sort(ranks, decreasing = TRUE)
cat("ranked genes:", length(ranks), "\n")
cat("t range:", round(range(ranks), 2), "\n")
cat("top 5 genes by t:", paste(names(head(ranks, 5)), collapse = ", "), "\n")
cat("bottom 5:", paste(names(tail(ranks, 5)), collapse = ", "), "\n")

cat("\nRunning fgsea over", length(fs), "muscle-expressed families ...\n")
res <- as.data.table(fgsea(pathways = fs, stats = ranks, minSize = 15, maxSize = Inf,
                           nPermSimple = 10000, eps = 0))
res[, absNES := abs(NES)]; setorder(res, padj, -absNES)
res[, rank_by_padj := .I]
res[, leadingEdge_chr := vapply(leadingEdge, function(x) paste(head(x, 50), collapse = ";"), "")]

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
fwrite(res[, .(family_id = pathway, size, NES, pval, padj, rank_by_padj,
               log2err, leading_edge_n = lengths(leadingEdge), leadingEdge_chr)],
       "results/tables/06_smoke_fgsea_DMDvsCtrl.tsv", sep = "\t")

m29_name <- grep("miR-29", names(fs), value = TRUE)
cat("\nmiR-29 family in the set list:", paste(m29_name, collapse = ", "), "\n")
r <- res[pathway %in% m29_name]
cat("\n================ SMOKE TEST: miR-29 in DMD vs control ================\n")
print(r[, .(family_id = pathway, size, NES = round(NES, 3), pval = signif(pval, 3),
            padj = signif(padj, 3), rank_by_padj)])

verdict <- if (nrow(r) && r$NES[1] > 0 && r$padj[1] < 0.05) "GO" else
           if (nrow(r) && r$NES[1] > 0 && r$pval[1] < 0.05) "BORDERLINE" else "NO-GO"
cat("\nVERDICT:", verdict, "\n")
cat("Direction rule (CLAUDE.md s1): NES > 0 means miR-29 targets bunch at the UP end of the\n")
cat("DMD-vs-control list, i.e. miR-29 is predicted LESS active in DMD. That is the published\n")
cat("direction (miR-29 is reduced in DMD), so NES > 0 is the result we needed to see.\n")

cat("\n== top 15 families by padj ==\n")
print(head(res[, .(family_id = pathway, size, NES = round(NES, 2),
                   padj = signif(padj, 3))], 15))
cat("\nfamilies at padj < 0.05:", res[padj < 0.05, .N], "of", nrow(res), "\n")
cat("  of those, NES > 0 (predicted less active):", res[padj < 0.05 & NES > 0, .N], "\n")
cat("  of those, NES < 0 (predicted more active):", res[padj < 0.05 & NES < 0, .N], "\n")

## ---- the enrichment plot (the "ranked walk") ------------------------------
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
if (nrow(r)) {
  p <- plotEnrichment(fs[[m29_name[1]]], ranks) +
    labs(title = sprintf("miR-29 target genes in DMD vs control  (NES = %.2f, padj = %.2g)",
                         r$NES[1], r$padj[1]),
         subtitle = paste0("x = all ", length(ranks),
                           " genes ranked by limma t (left = most up in DMD).\n",
                           "Each tick = one of the ", r$size[1],
                           " miR-29 targets. The curve climbs when ticks are dense."),
         x = "gene rank in the DMD-vs-control list", y = "running enrichment score")
  ggsave("results/figures/06_mir29_enrichment_DMDvsCtrl.png", p, width = 7.5, height = 5, dpi = 200)
  cat("Wrote results/figures/06_mir29_enrichment_DMDvsCtrl.png\n")
}

## ---- plain-numbers version of the same claim ------------------------------
# fgsea is a rank test, which is hard to feel. So also show the blunt version:
# the average logFC of miR-29 targets vs everything else.
tg29 <- fs[[m29_name[1]]]
a <- tt[gene_symbol %in% tg29, mean(logFC)]
b <- tt[!gene_symbol %in% tg29, mean(logFC)]
cat(sprintf("\nPlain version: mean DMD-vs-control log2FC of the %d miR-29 targets = %+0.3f\n",
            length(tg29), a))
cat(sprintf("               mean log2FC of the other %d genes              = %+0.3f\n",
            nrow(tt) - length(tg29), b))
cat(sprintf("               difference (the 'footprint size' ruler, s10b)  = %+0.3f log2 units\n", a - b))
wt <- wilcox.test(tt[gene_symbol %in% tg29, logFC], tt[!gene_symbol %in% tg29, logFC])
cat(sprintf("               Wilcoxon p for that shift = %.3g\n", wt$p.value))
saveRDS(list(res = res, ranks = ranks, mir29 = r, verdict = verdict,
             ruler_log2 = a - b), "results/rds/smoke_test.rds")
cat("Wrote results/rds/smoke_test.rds\n")
