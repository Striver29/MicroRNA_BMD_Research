# Fact-check register

Rule from CLAUDE.md §9: never invent a number, citation or claim. Anything unverified is
marked `[UNVERIFIED]` with the reason.

## Cleared in this session

| Claim | Status | Evidence |
|---|---|---|
| GSE109178 has 49 samples on GPL570, vastus lateralis | **CONFIRMED** | GEO sample records; `results/tables/01_sample_sheet.tsv`. All 49 are `GPL570`, all `vastus lateralis`. |
| 6 control / 11 BMD / 17 DMD / 15 LGMD | **CONFIRMED** | `scripts/01_download_labels.R` hard-asserts this and passed. LGMD splits 8 dysferlin (DYS) + 7 FKRP, matching "8 LGMD2B and 7 LGMD2I". |
| BMD ages run 2.5 to 60.8 years | **CONFIRMED** | min 2.5, max 60.8, median 12.2, all 11 known. |
| DMD patients are children | **CONFIRMED** | ages 0.9–8.0, median 3.5, all 17 known. |
| GSM2934804 severity | **RESOLVED as far as GEO goes** | GEO records it as `pathology: severe`, description `BMD-12-[6.4]-severe`, mutation `Deletion Exon 49`, age 6.4. The dispute is about whether that label is right, which GEO cannot settle. Kept as the step-13 robustness item. |
| miR-29 seed is `AGCACCA` | **CONFIRMED** | TargetScan 8.0 `miR_Family_Info.txt`: family `miR-29-3p`, seed `AGCACCA`, human members hsa-miR-29a-3p / 29b-3p / 29c-3p, conservation code 2 (highly conserved). |
| Eisenberg 2007: DMD COL3A1 ~6-fold up, COL1A2 ~5-fold up | **CONSISTENT with our data** | ours: COL3A1 logFC 3.26 = 9.6×, COL1A2 2.94 = 7.7×. Same direction, same order of magnitude. Different platform/normalisation, so "about" is the right word. |

## Still unverified — do not put in the paper yet

| Item | Why it is unverified |
|---|---|
| Control ages and sexes in GSE109178 | `[UNVERIFIED]` — **not reported in GEO at all.** All 6 control records have `age: NA` and no sex. CLAUDE.md §2 asked us to record them; they do not exist in the source. The age confound therefore **cannot be resolved** and must be stated as a limitation. See §Age confound below. |
| miRNA Tissue Atlas 2025 exact citation | `[UNVERIFIED]` — the download page lists three DOIs (gkw116, gkab808, gkae1036); which one is the 2025 release has not been checked against the papers. |
| Rubenstein et al. 2020 PMID 32042044 | `[UNVERIFIED]` — PMID taken from the MSigDB C8 set metadata, not checked against the paper. |
| Reschke et al. 2026 pig microRNA list and directions | `[UNVERIFIED]` — Table 1 not obtained. Step 15 BLOCKED. |
| Eisenberg "5 BMD", "GSE3307 = same dataset" | still not found. Do not use. |
| Xie 2024 gene counts; the "13 genes" claim; Dadgar 2014 as the source of GSE109178; the pig "18 upregulated"; GSE157668 contents | unchanged from CLAUDE.md §11. Do not use. |

## Age confound (must appear in Limitations)

- Control ages: **unknown**, not in GEO.
- BMD median 12.2 y (range 2.5–60.8). DMD median 3.5 y (range 0.9–8.0).
- So the **direct BMD-vs-DMD comparison confounds disease severity with age**: the BMD group is
  roughly a decade older. Any BMD-vs-DMD difference may be partly a child-vs-older-patient
  difference. We cannot adjust for this, because adjusting would also remove the disease
  signal and because control ages are missing entirely.
- State it; do not try to fix it (CLAUDE.md §2).
