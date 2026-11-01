#!/usr/bin/env Rscript
# 04_limma.R -- limma on the gene matrix, three contrasts.  CLAUDE.md s6 step 4.
suppressPackageStartupMessages({library(data.table); library(limma)})
set.seed(1)

x <- readRDS("data/processed/gene_matrix.rds"); mat <- x$mat; sheet <- as.data.table(x$sheet)
setkey(sheet, gsm); sheet <- sheet[colnames(mat)]
stopifnot(identical(sheet$gsm, colnames(mat)))

grp <- factor(sheet$group, levels = c("control", "BMD", "DMD"))
cat("design groups:\n"); print(table(grp))
design <- model.matrix(~ 0 + grp); colnames(design) <- levels(grp)

fit <- lmFit(mat, design)
cm <- makeContrasts(BMDvsCtrl = BMD - control,
                    DMDvsCtrl = DMD - control,
                    BMDvsDMD  = BMD - DMD,
                    levels = design)
fit2 <- eBayes(contrasts.fit(fit, cm))
cat("\nprior df (eBayes moderation):", round(fit2$df.prior, 2),
    " residual df:", unique(fit2$df.residual), "\n")

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/rds", recursive = TRUE, showWarnings = FALSE)
tt <- list()
for (cn in colnames(cm)) {
  t_ <- as.data.table(topTable(fit2, coef = cn, number = Inf, sort.by = "none"),
                      keep.rownames = "gene_symbol")
  t_ <- t_[, .(gene_symbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)]
  setorder(t_, -t)
  fwrite(t_, sprintf("results/tables/limma_%s.tsv", cn), sep = "\t")
  tt[[cn]] <- t_
  cat(sprintf("%-10s FDR<0.05: up %5d  down %5d  (of %d genes)\n", cn,
              t_[adj.P.Val < 0.05 & logFC > 0, .N],
              t_[adj.P.Val < 0.05 & logFC < 0, .N], nrow(t_)))
}
saveRDS(list(fit = fit2, tables = tt, sheet = sheet), "results/rds/limma.rds")
cat("Wrote results/rds/limma.rds and results/tables/limma_*.tsv\n")

## ---- sanity check: the collagen/crosslink genes ---------------------------
san_genes <- c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "LOX")
san <- rbindlist(lapply(names(tt), function(cn)
  tt[[cn]][gene_symbol %in% san_genes,
           .(contrast = cn, gene_symbol, logFC = round(logFC, 3),
             t = round(t, 2), adj.P.Val = signif(adj.P.Val, 3))]))
san <- dcast(san, gene_symbol ~ contrast, value.var = c("logFC", "adj.P.Val"))
cat("\n== SANITY: fibrosis genes (logFC on the log2 scale; +1 = 2x up) ==\n")
print(san)
fwrite(san, "results/tables/04_sanity_fibrosis_genes.tsv", sep = "\t")

cat("\nRead it like this: logFC_DMDvsCtrl = +2 means that gene sits 2^2 = 4x higher in DMD than control.\n")
