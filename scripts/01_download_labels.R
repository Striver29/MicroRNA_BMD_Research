#!/usr/bin/env Rscript
# 01_download_labels.R -- get GSE109178 metadata + CEL files, build the sample sheet.
# CLAUDE.md s6 step 1.  Confirms 6 control / 11 BMD / 17 DMD / 15 LGMD and records control ages.
#
# NOTE (deviation, logged in docs/methods_log.md): ftp.ncbi.nlm.nih.gov returns HTTP 403 for
# file objects from this machine, so GEOquery::getGEO()/getGEOSuppFiles() cannot reach the
# series matrix or the supplementary tar. We therefore use the two official GEO web endpoints:
#   metadata : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109178&targ=gsm&form=text&view=brief
#   CEL tar  : https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE109178&format=file
# These serve byte-identical content to the FTP mirror. The parsed fields are exactly the ones
# pData() would expose (title, source_name_ch1, characteristics_ch1, description, supplementary_file).

suppressPackageStartupMessages(library(data.table))

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("inputs/auto", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

GSE <- "GSE109178"
soft_txt <- file.path("data/raw", paste0(GSE, "_gsm_brief.txt"))
raw_tar  <- file.path("data/raw", paste0(GSE, "_RAW.tar"))
cel_dir  <- "data/raw/CEL"

## ---- 1. metadata -----------------------------------------------------------
if (!file.exists(soft_txt)) {
  url <- sprintf("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s&targ=gsm&form=text&view=brief", GSE)
  download.file(url, soft_txt, quiet = TRUE)
}
txt <- readLines(soft_txt, warn = FALSE)
txt <- gsub("\r$", "", txt)

# split into one block per sample
starts <- grep("^\\^SAMPLE = ", txt)
stopifnot(length(starts) == 49L)
ends <- c(starts[-1] - 1L, length(txt))

field <- function(blk, key) {
  hit <- grep(paste0("^!Sample_", key, " = "), blk, value = TRUE)
  if (!length(hit)) return(NA_character_)
  sub(paste0("^!Sample_", key, " = "), "", hit[1])
}
charfield <- function(blk, key) {
  hit <- grep("^!Sample_characteristics_ch1 = ", blk, value = TRUE)
  hit <- sub("^!Sample_characteristics_ch1 = ", "", hit)
  kv  <- hit[grepl(":", hit, fixed = TRUE)]
  k   <- tolower(trimws(sub(":.*$", "", kv)))
  v   <- trimws(sub("^[^:]*:", "", kv))
  out <- v[k == tolower(key)]
  if (!length(out)) NA_character_ else out[1]
}

meta <- rbindlist(lapply(seq_along(starts), function(i) {
  blk <- txt[starts[i]:ends[i]]
  data.table(
    gsm         = sub("^\\^SAMPLE = ", "", blk[1]),
    title       = field(blk, "title"),
    description = field(blk, "description"),
    tissue      = field(blk, "source_name_ch1"),
    platform    = field(blk, "platform_id"),
    age_raw     = charfield(blk, "age"),
    sex_raw     = charfield(blk, "Sex"),
    severity_raw= charfield(blk, "pathology"),
    mutation    = charfield(blk, "mutation"),
    path_note   = charfield(blk, "pathology note"),
    cel_file    = basename(field(blk, "supplementary_file"))
  )
}))

## ---- 2. groups -------------------------------------------------------------
# Group comes from !Sample_description, which is the only unambiguous label:
#   "BMD-<n>-[age]-<sev>"  "Normal"  "DMD-<n>-[age]-<sev>"
#   "DYS-..."  = dysferlinopathy = LGMD2B (LGMD R2)
#   "FKRP-..." = FKRP-related    = LGMD2I (LGMD R9)
# !Sample_title is series boilerplate ("Asynchronous remodeling in DMD (...)") for every
# sample, so it must NOT be used for grouping.
meta[, group := fifelse(grepl("^BMD",    description), "BMD",
                 fifelse(grepl("^Normal", description), "control",
                 fifelse(grepl("^DMD",    description), "DMD",
                 fifelse(grepl("^DYS|^FKRP", description), "LGMD", NA_character_))))]
meta[, lgmd_subtype := fifelse(grepl("^DYS",  description), "LGMD2B_dysferlin",
                        fifelse(grepl("^FKRP", description), "LGMD2I_FKRP", NA_character_))]

meta[, age_years := suppressWarnings(as.numeric(age_raw))]
meta[, sex := fifelse(grepl("^M", sex_raw, ignore.case = TRUE) & !grepl("^Mut", sex_raw), "M",
              fifelse(grepl("^F", sex_raw, ignore.case = TRUE), "F", NA_character_))]
meta[, severity := fifelse(severity_raw %in% c("mild"), "mild",
                    fifelse(severity_raw %in% c("moderate", "mod"), "moderate",
                    fifelse(severity_raw %in% c("severe"), "severe", NA_character_)))]
meta[, include := group %in% c("control", "BMD", "DMD")]
meta[, notes := ""]
meta[gsm == "GSM2934804", notes := "disputed severity (GEO says severe); robustness item, s6 step 13"]
meta[group == "control", notes := "age and sex not reported in GEO -> potential age confound"]

## ---- 3. the hard check -----------------------------------------------------
tab <- table(meta$group)
cat("\n== group counts (all 49) ==\n"); print(tab)
cat("\n== LGMD breakdown ==\n"); print(table(meta$lgmd_subtype))
expected <- c(BMD = 11L, control = 6L, DMD = 17L, LGMD = 15L)
got <- as.integer(tab[names(expected)])
if (!identical(got, as.integer(expected))) {
  stop("GROUP COUNTS DO NOT MATCH 6/11/17/15. Got: ", paste(names(expected), got, collapse = " "))
}
cat("\nCounts match the spec: 6 control / 11 BMD / 17 DMD / 15 LGMD.\n")
cat("Included samples:", sum(meta$include), "\n")
stopifnot(sum(meta$include) == 34L)
cat("All on platform(s):", paste(unique(meta$platform), collapse = ","), "\n")
cat("All tissue(s):", paste(unique(meta$tissue), collapse = ","), "\n")

cat("\n== ages by group ==\n")
print(meta[, .(n = .N,
               n_age_known = sum(!is.na(age_years)),
               min_age = suppressWarnings(min(age_years, na.rm = TRUE)),
               median_age = suppressWarnings(as.numeric(median(age_years, na.rm = TRUE))),
               max_age = suppressWarnings(max(age_years, na.rm = TRUE))), by = group])
cat("\nCONTROL AGES:", if (all(is.na(meta[group == "control", age_years])))
    "NOT REPORTED IN GEO -> age confound stays unresolved [UNVERIFIED]" else "see above", "\n")

cat("\n== BMD severity ==\n"); print(table(meta[group == "BMD", severity], useNA = "ifany"))
cat("\n== DMD severity ==\n"); print(table(meta[group == "DMD", severity], useNA = "ifany"))

## ---- 4. CEL files ----------------------------------------------------------
if (!file.exists(raw_tar)) {
  download.file(sprintf("https://www.ncbi.nlm.nih.gov/geo/download/?acc=%s&format=file", GSE),
                raw_tar, quiet = TRUE, mode = "wb")
}
stopifnot(file.exists(raw_tar))
inside <- untar(raw_tar, list = TRUE)
stopifnot(length(inside) == 49L)
want <- meta[include == TRUE, cel_file]
stopifnot(all(want %in% inside))
dir.create(cel_dir, recursive = TRUE, showWarnings = FALSE)
todo <- want[!file.exists(file.path(cel_dir, want))]
if (length(todo)) {
  cat("\nExtracting", length(todo), "CEL files (only the 34 included) ...\n")
  untar(raw_tar, files = todo, exdir = cel_dir)
}
for (f in want) {
  gz <- file.path(cel_dir, f); plain <- sub("\\.gz$", "", gz)
  if (file.exists(gz) && !file.exists(plain)) R.utils::gunzip(gz, destname = plain, remove = TRUE)
}
have <- sum(file.exists(file.path(cel_dir, sub("\\.gz$", "", want))))
cat("CEL files ready in", cel_dir, ":", have, "of 34\n")

## ---- 5. write the sample sheet --------------------------------------------
sheet <- meta[, .(gsm, title = description, group, age_years, sex, severity,
                  include, notes, cel_file = sub("\\.gz$", "", cel_file),
                  lgmd_subtype, mutation, geo_title = title, path_note)]
setorder(sheet, group, gsm)

out <- if (file.exists("inputs/sample_sheet.csv")) "inputs/sample_sheet.csv" else "inputs/auto/sample_sheet.csv"
if (file.exists("inputs/sample_sheet.csv")) {
  cat("\nHamza's inputs/sample_sheet.csv exists -- NOT overwriting. Writing ours to inputs/auto/ for diffing.\n")
  out <- "inputs/auto/sample_sheet.csv"
}
fwrite(sheet, out)
cat("Wrote", out, "\n")
fwrite(sheet, "results/tables/01_sample_sheet.tsv", sep = "\t")
cat("Wrote results/tables/01_sample_sheet.tsv\n")
