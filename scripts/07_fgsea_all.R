#!/usr/bin/env Rscript
# 07_fgsea_all.R -- fgsea for all three contrasts + the leading-edge network. CLAUDE.md s6 step 7.
suppressPackageStartupMessages({library(data.table); library(fgsea); library(ggplot2)})
set.seed(42)

fs  <- readRDS("results/rds/family_sets.rds")
lim <- readRDS("results/rds/limma.rds")
pick <- function(f) { a <- file.path("inputs", f); if (file.exists(a)) a else file.path("inputs/auto", f) }
fib <- fread(pick("fibrosis_sets.tsv"))
fib_genes <- unique(fib$gene_symbol)

contrasts <- c("BMDvsCtrl", "DMDvsCtrl", "BMDvsDMD")
run1 <- function(cn) {
  tt <- lim$tables[[cn]]
  ranks <- sort(setNames(tt$t, tt$gene_symbol), decreasing = TRUE)
  r <- as.data.table(fgsea(fs, ranks, minSize = 15, maxSize = Inf, nPermSimple = 10000, eps = 0))
  r[, contrast := cn][]
}
all_res <- rbindlist(lapply(contrasts, run1))
saveRDS(all_res, "results/rds/fgsea_all.rds")

## ---- wide table: families as rows, NES/padj per comparison ----------------
wide <- dcast(all_res, pathway + size ~ contrast, value.var = c("NES", "padj"))
setnames(wide, "pathway", "family_id")
setcolorder(wide, c("family_id", "size",
                    paste0("NES_", contrasts), paste0("padj_", contrasts)))
setorder(wide, padj_DMDvsCtrl)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
fwrite(wide, "results/tables/07_fgsea_wide.tsv", sep = "\t")
cat("Wrote results/tables/07_fgsea_wide.tsv  (", nrow(wide), "families )\n")

cat("\n== significant families per comparison (fgsea padj < 0.05) ==\n")
print(all_res[, .(n_sig = sum(padj < 0.05), sig_NESpos = sum(padj < 0.05 & NES > 0),
                  sig_NESneg = sum(padj < 0.05 & NES < 0), n_tested = .N), by = contrast])

cat("\n== miR-29 across the three comparisons ==\n")
print(all_res[pathway == "miR-29-3p",
              .(contrast, size, NES = round(NES, 3), pval = signif(pval, 3), padj = signif(padj, 3))])

cat("\n== top 20 families by DMD-vs-control padj, side by side ==\n")
print(head(wide[, .(family_id, size,
                    NES_BMD = round(NES_BMDvsCtrl, 2), padj_BMD = signif(padj_BMDvsCtrl, 2),
                    NES_DMD = round(NES_DMDvsCtrl, 2), padj_DMD = signif(padj_DMDvsCtrl, 2),
                    NES_BvD = round(NES_BMDvsDMD, 2),  padj_BvD = signif(padj_BMDvsDMD, 2))], 20))

## ---- leading-edge table = the microRNA -> fibrosis-gene network ------------
top_fams <- head(wide$family_id, 25)
le <- all_res[pathway %in% top_fams, .(contrast, family_id = pathway, NES, padj,
                                       gene = leadingEdge)][, .(gene = unlist(gene)),
                                       by = .(contrast, family_id, NES, padj)]
le[, in_fibrosis_set := gene %in% fib_genes]
le <- merge(le, fib[, .(gene = gene_symbol, set_name)], by = "gene", all.x = TRUE, allow.cartesian = TRUE)
setorder(le, contrast, padj, family_id, -in_fibrosis_set, gene)
fwrite(le, "results/tables/07_leading_edge.tsv", sep = "\t")
cat("\nWrote results/tables/07_leading_edge.tsv  (", nrow(le), "rows )\n")

cat("\n== leading-edge fibrosis-gene content, top families, DMD vs control ==\n")
summ <- le[contrast == "DMDvsCtrl", .(leading_edge_n = uniqueN(gene),
                                      n_fibrosis = uniqueN(gene[in_fibrosis_set]),
                                      pct_fibrosis = round(100 * uniqueN(gene[in_fibrosis_set]) / uniqueN(gene), 1)),
           by = .(family_id, NES = round(NES, 2), padj = signif(padj, 2))]
setorder(summ, padj)
print(head(summ, 15))
fwrite(summ, "results/tables/07_leading_edge_summary.tsv", sep = "\t")

cat("\n== miR-29 leading-edge fibrosis genes in DMD vs control ==\n")
cat(paste(sort(unique(le[contrast == "DMDvsCtrl" & family_id == "miR-29-3p" & in_fibrosis_set == TRUE, gene])),
          collapse = ", "), "\n")

## ---- heatmap of NES for the top families ----------------------------------
hm <- melt(wide[family_id %in% top_fams, .(family_id, NES_BMDvsCtrl, NES_DMDvsCtrl, NES_BMDvsDMD)],
           id.vars = "family_id", variable.name = "contrast", value.name = "NES")
hm[, contrast := factor(sub("^NES_", "", contrast), levels = contrasts)]
ord <- wide[family_id %in% top_fams][order(NES_DMDvsCtrl), family_id]
hm[, family_id := factor(family_id, levels = ord)]
ph <- ggplot(hm, aes(contrast, family_id, fill = NES)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", NES)), size = 2.5) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  labs(x = NULL, y = NULL, title = "Predicted microRNA footprints (NES) for the top 25 families",
       subtitle = "NES > 0 (red) = targets bunch at the UP end = family predicted LESS active.\nNES < 0 (blue) = predicted MORE active.") +
  theme_bw(base_size = 8)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("results/figures/07_nes_heatmap.png", ph, width = 7, height = 7.5, dpi = 200)
cat("Wrote results/figures/07_nes_heatmap.png\n")

## ---- the real r and slope for step 8 --------------------------------------
xy <- wide[, .(family_id, x = NES_DMDvsCtrl, y = NES_BMDvsCtrl)]
r_real <- cor(xy$x, xy$y)
sl_real <- coef(lm(y ~ x, data = xy))[2]
cat(sprintf("\n== BMD vs DMD footprint agreement across %d families ==\n", nrow(xy)))
cat(sprintf("Pearson r (BMD NES vs DMD NES) = %.4f\n", r_real))
cat(sprintf("slope (BMD NES ~ DMD NES)      = %.4f\n", sl_real))
cat("These two numbers MUST be judged against the shuffled-label baseline (step 8), not against 0.\n")
saveRDS(list(xy = xy, r_real = r_real, slope_real = as.numeric(sl_real)),
        "results/rds/bmd_dmd_agreement.rds")
