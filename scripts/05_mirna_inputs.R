#!/usr/bin/env Rscript
# 05_mirna_inputs.R -- build the family -> target-gene sets used by fgsea.  CLAUDE.md s6 step 5.
suppressPackageStartupMessages(library(data.table))
set.seed(1)
MIN_TARGETS <- 15L            # CLAUDE.md s10 default: >= 15 targets after intersecting, no cap

pick <- function(f) {
  a <- file.path("inputs", f); b <- file.path("inputs/auto", f)
  p <- if (file.exists(a)) a else b
  if (!file.exists(p)) stop("missing input: ", f)
  cat("using", p, "\n"); p
}
fam <- fread(pick("mirna_families.tsv"))
tg  <- fread(pick("mirna_targets.tsv"))

gm <- readRDS("data/processed/gene_matrix.rds")$mat
genes_on_array <- rownames(gm)
cat("genes on the array:", length(genes_on_array), "\n")

## ---- keep muscle-expressed families only ----------------------------------
fam_keep <- fam[expressed_in_muscle == TRUE]
cat("families total:", nrow(fam), " muscle-expressed:", nrow(fam_keep), "\n")

tg <- tg[family_id %in% fam_keep$family_id]
cat("target rows for muscle families:", nrow(tg),
    " families with any target:", uniqueN(tg$family_id), "\n")

## ---- keep only targets that exist in our gene matrix -----------------------
tg <- tg[gene_symbol %in% genes_on_array]
cat("target rows after intersecting with the array:", nrow(tg), "\n")

sizes <- tg[, .(n_targets = uniqueN(gene_symbol)), by = family_id]
setorder(sizes, -n_targets)
cat("families with >=", MIN_TARGETS, "targets:", sizes[n_targets >= MIN_TARGETS, .N],
    " dropped as too small:", sizes[n_targets < MIN_TARGETS, .N], "\n")
keep <- sizes[n_targets >= MIN_TARGETS, family_id]

family_sets <- split(tg[family_id %in% keep, unique(gene_symbol)],
                     tg[family_id %in% keep, family_id][
                       !duplicated(paste(tg[family_id %in% keep, family_id],
                                         tg[family_id %in% keep, gene_symbol]))])
family_sets <- lapply(split(tg[family_id %in% keep], by = "family_id", keep.by = FALSE),
                      function(d) sort(unique(d$gene_symbol)))
cat("\nfamily_sets:", length(family_sets), "families\n")
cat("set sizes: min", min(lengths(family_sets)), " median", median(lengths(family_sets)),
    " max", max(lengths(family_sets)), "\n")

dir.create("results/rds", recursive = TRUE, showWarnings = FALSE)
saveRDS(family_sets, "results/rds/family_sets.rds")
cat("Wrote results/rds/family_sets.rds\n")

out <- merge(sizes[family_id %in% keep],
             fam_keep[, .(family_id, seed_7mer, members, conservation, max_muscle_rpmm)],
             by = "family_id")
setorder(out, -n_targets)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
fwrite(out, "results/tables/05_family_sizes.tsv", sep = "\t")
cat("Wrote results/tables/05_family_sizes.tsv\n")

cat("\n== 10 largest families ==\n"); print(head(out[, .(family_id, seed_7mer, n_targets, conservation, max_muscle_rpmm = round(max_muscle_rpmm))], 10))
cat("\n== miR-29 (the positive control) ==\n")
m29 <- out[seed_7mer == "AGCACCA"]
print(m29[, .(family_id, seed_7mer, n_targets, members, conservation, max_muscle_rpmm = round(max_muscle_rpmm))])
if (!nrow(m29)) stop("miR-29 (seed AGCACCA) did not survive the filters -- the smoke test cannot run.")
