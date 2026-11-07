#!/usr/bin/env Rscript
# build_inputs.R -- build the inputs Hamza owns (CLAUDE.md s4) into inputs/auto/.
# Never overwrites anything in inputs/. Every file built here is logged in docs/handoffs.md.
suppressPackageStartupMessages({library(data.table); library(msigdbr)})
set.seed(1)
dir.create("inputs/auto", recursive = TRUE, showWarnings = FALSE)
TS <- "data/raw/targetscan80"
w <- function(dt, f) { p <- file.path("inputs/auto", f); fwrite(dt, p, sep = "\t")
                       cat(sprintf("wrote %-28s %6d rows\n", p, nrow(dt))) }

## ====================== 1. muscle-expressed miRNAs =========================
# Source: miRNA Tissue Atlas (CCB Saarland, 2025 release), human harmonised RPMM table.
# Skeletal muscle = Biotype "tissue", Organ_system "limb_muscle" (tissues limb_muscle +
# fast_twitch). Those are limb skeletal muscle, the same kind of muscle as vastus lateralis.
# A miRNA counts as muscle-expressed if its mean RPMM reaches MUSCLE_RPMM_MIN in either tissue.
MUSCLE_RPMM_MIN <- 10
atlas_file <- "data/raw/mirna_atlas/all_human_miRNA_rpmm_harmonized_meta_data.parquet"
muscle_csv <- "data/raw/mirna_atlas/muscle_rpmm.csv"   # written by build_inputs.py helper
stopifnot(file.exists(muscle_csv))
mus <- fread(muscle_csv)                               # cols: mirna, rpmm_muscle_max
mus[, expressed := rpmm_muscle_max >= MUSCLE_RPMM_MIN]
cat("muscle atlas: ", nrow(mus), " miRNAs, ", sum(mus$expressed),
    " at RPMM >= ", MUSCLE_RPMM_MIN, "\n", sep = "")

## ====================== 2. mirna_families.tsv ==============================
fam <- fread(file.path(TS, "miR_Family_Info.txt"), sep = "\t", quote = "")
setnames(fam, c("miR family", "Seed+m8", "Species ID", "MiRBase ID",
                "Mature sequence", "Family Conservation?", "MiRBase Accession"),
              c("family_id", "seed_7mer", "species", "mirbase_id",
                "mature_seq", "conservation", "mirbase_acc"), skip_absent = TRUE)
fam <- fam[species == 9606L]
cat("TargetScan 8.0 human miRNAs:", nrow(fam), " families:", uniqueN(fam$family_id), "\n")

fam <- merge(fam, mus[, .(mirbase_id = mirna, rpmm_muscle_max, expressed)],
             by = "mirbase_id", all.x = TRUE)
fam[is.na(expressed), `:=`(expressed = FALSE, rpmm_muscle_max = NA_real_)]

fam_tbl <- fam[, .(
  seed_7mer            = seed_7mer[1],
  members              = paste(sort(unique(mirbase_id)), collapse = ";"),
  n_members            = uniqueN(mirbase_id),
  conservation         = max(conservation),
  max_muscle_rpmm      = suppressWarnings(max(rpmm_muscle_max, na.rm = TRUE)),
  n_members_in_atlas   = sum(!is.na(rpmm_muscle_max)),
  expressed_in_muscle  = any(expressed),
  expression_source    = "miRNA Tissue Atlas 2025 (CCB Saarland), human harmonised RPMM, Organ_system=limb_muscle, mean RPMM >= 10"
), by = family_id]
fam_tbl[!is.finite(max_muscle_rpmm), max_muscle_rpmm := NA_real_]
setorder(fam_tbl, -expressed_in_muscle, -max_muscle_rpmm)
cat("families total:", nrow(fam_tbl), " muscle-expressed:", sum(fam_tbl$expressed_in_muscle), "\n")
cat("\n-- miR-29 check --\n")
print(fam_tbl[seed_7mer == "AGCACCA",
              .(family_id, seed_7mer, members, conservation, max_muscle_rpmm, expressed_in_muscle)])
w(fam_tbl, "mirna_families.tsv")

## ====================== 3. mirna_targets.tsv ===============================
tg <- fread(file.path(TS, "Predicted_Targets_Info.default_predictions.txt"),
            sep = "\t", quote = "", showProgress = FALSE)
setnames(tg, c("miR Family", "Gene Symbol", "Species ID", "Seed match"),
             c("family_id", "gene_symbol", "species", "site_type"), skip_absent = TRUE)
tg <- tg[species == 9606L, .(family_id, gene_symbol, site_type)]
tg <- unique(tg[, .(family_id, gene_symbol = toupper(gene_symbol), site_type)])
tg[, source := "TargetScanHuman 8.0, Predicted_Targets_Info.default_predictions.txt (conserved sites)"]
cat("human target rows:", nrow(tg), " genes:", uniqueN(tg$gene_symbol),
    " families:", uniqueN(tg$family_id), "\n")
cat("miR-29-3p targets:", tg[family_id == "miR-29-3p", uniqueN(gene_symbol)], "\n")
w(tg, "mirna_targets.tsv")

## ====================== 4. utr_lengths.tsv =================================
# TargetScan's UTR sequences are multiple-alignment rows, so '-' gap characters must go
# before counting.  Example: "GU--UUC-" has 8 characters but a true length of 5.
utr <- fread(file.path(TS, "UTR_Sequences_human9606.txt"), sep = "\t", quote = "",
             showProgress = FALSE)
setnames(utr, c("Refseq ID", "Gene Symbol", "UTR sequence"),
              c("transcript_used", "gene_symbol", "seq"), skip_absent = TRUE)
utr[, utr3_length_nt := nchar(gsub("-", "", seq, fixed = TRUE))]
utr <- utr[, .(gene_symbol = toupper(gene_symbol), utr3_length_nt, transcript_used)]
setorder(utr, gene_symbol, -utr3_length_nt)
utr <- utr[, .SD[1], by = gene_symbol]          # one representative transcript per gene
cat("UTR lengths for", nrow(utr), "genes; median", median(utr$utr3_length_nt),
    "nt, range", min(utr$utr3_length_nt), "-", max(utr$utr3_length_nt), "\n")
w(utr, "utr_lengths.tsv")

## ====================== 5. fibrosis_sets.tsv ===============================
mdb_ver <- as.character(msigdbr::msigdbr_collections()$db_version[1])
cgp <- as.data.table(msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CGP"))
hal <- as.data.table(msigdbr(species = "Homo sapiens", collection = "H"))
core <- unique(cgp[gs_name == "NABA_CORE_MATRISOME", toupper(gene_symbol)])
tgfb <- unique(hal[gs_name == "HALLMARK_TGF_BETA_SIGNALING", toupper(gene_symbol)])
stopifnot(length(core) > 200, length(tgfb) > 40)
fib <- rbind(
  data.table(set_name = "matrisome_core", gene_symbol = sort(core),
             source_version = sprintf("MSigDB %s C2:CGP NABA_CORE_MATRISOME (Naba et al. 2012, PMID 22159717)", mdb_ver)),
  data.table(set_name = "hallmark_tgfb", gene_symbol = sort(tgfb),
             source_version = sprintf("MSigDB %s H HALLMARK_TGF_BETA_SIGNALING", mdb_ver)))
cat("matrisome_core:", length(core), " hallmark_tgfb:", length(tgfb),
    " overlap:", length(intersect(core, tgfb)), "\n")
w(fib, "fibrosis_sets.tsv")

## ====================== 6. fibroblast_markers.tsv ==========================
# Published skeletal-muscle source: Rubenstein et al. 2020 single-cell atlas of human
# skeletal muscle, as curated in MSigDB C8. FAP = fibro-adipogenic progenitor, the
# fibroblast population of skeletal muscle.
c8 <- as.data.table(msigdbr(species = "Homo sapiens", collection = "C8"))
fapsets <- c("RUBENSTEIN_SKELETAL_MUSCLE_FAP_CELLS", "RUBENSTEIN_SKELETAL_MUSCLE_FBN1_FAP_CELLS")
fap <- unique(c8[gs_name %in% fapsets, toupper(gene_symbol)])
stopifnot(length(fap) > 100)
fmk <- data.table(gene_symbol = sort(fap),
                  source = sprintf("MSigDB %s C8: %s (Rubenstein et al. 2020, PMID 31937892)",
                                   mdb_ver, paste(fapsets, collapse = " + ")))
fmk[, in_fibrosis_set := gene_symbol %in% c(core, tgfb)]
cat("fibroblast markers:", nrow(fmk), " of which also in a fibrosis set:",
    sum(fmk$in_fibrosis_set), " -> usable for the covariate:", sum(!fmk$in_fibrosis_set), "\n")
w(fmk, "fibroblast_markers.tsv")

## ====================== 7. pig_mirnas.tsv ==================================
# CLAUDE.md s4: build this ONLY from Reschke et al. 2026 Table 1. Not obtained -> step 15 blocked.
cat("\npig_mirnas.tsv NOT built: Reschke et al. 2026 Table 1 not obtained. Step 15 stays BLOCKED.\n")
cat("\nDone.\n")
