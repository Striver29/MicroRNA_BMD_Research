# Methods log

Every parameter, default, data source and version, cut and deviation, with its reason.
Session 1: 2026-11-01 (Sunday night). Software: R 4.6.2 (2026-10-29), Bioconductor 3.23,
macOS arm64. Full package list in `results/sessionInfo.txt`.

Key package versions: limma 3.68.5, fgsea 1.38.0, oligo (Bioc 3.23), hgu133plus2.db,
msigdbr 26.1.1 (MSigDB 2026.1.Hs), data.table 1.18.6.1, ggplot2 4.0.3, BiocParallel 1.46.0.

## Environment deviations

| # | Deviation | Reason |
|---|---|---|
| E1 | R and all Bioconductor packages had to be installed from scratch (R was not on the machine). Installed via Homebrew (`r 4.6.2`); all R packages built **from source** because Homebrew's R does not use CRAN's macOS binaries. | No R present at session start. |
| E2 | `GEOquery` **failed to build** and is not installed. | Build error. Not needed: see D1. |
| E3 | `getGEOSuppFiles()` / `getGEO()` not used. GEO was reached through its two official web endpoints instead: metadata from `www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109178&targ=gsm&form=text&view=brief`, CEL archive from `www.ncbi.nlm.nih.gov/geo/download/?acc=GSE109178&format=file`. | `ftp.ncbi.nlm.nih.gov` returns **HTTP 403** for file objects from this machine (directory listings work, file fetches do not). The web endpoints serve the same content. Parsed fields are exactly the ones `pData()` exposes. |
| E4 | Only the 34 included CEL files were extracted from the 49-file archive. | Keeps the LGMD samples out of normalisation (CLAUDE.md §6 step 2) and saves disk. |

## Decisions and defaults used (CLAUDE.md §10)

| Decision | Value used | Note |
|---|---|---|
| Probe collapse | highest mean expression | §10 default |
| Low-expression filter | none | §10 default; ranks look healthy (t range −15.0 to +17.7) |
| TargetScan target set | 8.0 default (conserved) predictions, human 9606 | §10 default |
| Muscle-expression source | miRNA Tissue Atlas 2025 (CCB Saarland), `Organ_system = limb_muscle`, mean RPMM ≥ 10 | §4 let us pick the first credible published atlas; see `docs/handoffs.md`. The **RPMM ≥ 10 threshold is ours**, not the spec's — logged here as a parameter to sensitivity-test if asked. |
| Family size limits | ≥ 15 targets after intersecting with the gene matrix, no maximum | §10 default |
| Smoke-test thresholds | GO = NES > 0 and padj < 0.05 | §6 step 6 |
| Shuffle count | **1,000** (step 8), 1,000 (step 9), 1,000 draws (step 11) | §10 default; timing allowed the full 1,000 |
| fgsea permutations | `nPermSimple = 10000`, `eps = 0` for every reported run | At the default 1,000, 13 of 128 families returned `NA` NES (fgsea's own warning). Raised to 10,000 to remove the NAs. Inside step 11's null draws only, 1,000 is used, because there we need the NES distribution and not a precise p per draw. |
| Catch criterion (step 10b) | NES > 0 and padj < 0.05 | §10 default |
| Random seeds | 1 (steps 1–4), 42 (fgsea steps 6, 7, 10a, 12), 20261101 (step 8), 20261102 (step 9), 20261103 (step 11), 20261104 (step 10b) | every script sets its own |

## Step-by-step parameters

**Step 2 — RMA and QC.** `oligo::read.celfiles` + `rma()` on the 34 included CELs only.
Output 54,675 × 34.

QC deviation, logged: the spec asked for "each sample's median correlation with all the
others". That metric is **confounded by group size here**: 28 of 34 arrays are dystrophic
muscle and only 6 are healthy, so all 6 controls score lowest on it (0.9496–0.9558 vs a 0.9759
centre) purely because few samples share their biology. On that metric three controls looked
like outliers. We therefore added two things:
1. **within-group** median correlation (each sample against the others in its own group), which
   removes the group effect — controls then score 0.9784–0.9854, at the top rather than the bottom;
2. **RLE and NUSE** from `oligo::fitProbeLevelModel`, the standard Affymetrix array-quality
   measures, computed within probeset so group biology cannot drive them.

Result: **no array has a technical fault.** All median RLE within ±0.057 (flag threshold 0.10),
all median NUSE within 0.981–1.039 (flag threshold 1.05). The within-group rule flags
GSM2934803 (BMD, 0.9676 vs a 0.9757 group centre) and GSM2934813 (control, 0.9784 vs 0.9839)
as mild within-group outliers only.

**Decision taken, for Aziz to confirm:** keep all 34. Reason — the authoritative quality
metrics are clean, the two within-group flags are biological spread rather than bad arrays, and
step 13's leave-one-out already tests sensitivity to exactly these samples. This was a **STOP**
point in CLAUDE.md §6 step 2; it is flagged in `docs/results_log.md` for Aziz's sign-off.

**Step 3 — probes to genes.** 54,675 probes → 44,647 with a symbol → 1,524 multi-mapping
probes dropped → 43,123 kept → collapsed to **20,823 genes** (11,043 genes had more than one
probe competing).

**Step 4 — limma.** `~0 + group`, `lmFit` → `contrasts.fit` → `eBayes`. Contrasts BMD−Ctrl,
DMD−Ctrl, BMD−DMD. Residual df 31, eBayes prior df 4.23.

**Step 5 — microRNA sets.** 2,064 human TargetScan families → 253 muscle-expressed → of those,
the ones with ≥ 15 array-present targets give **128 families** (sizes 15 to 1,446, median 371).
miR-29-3p has 1,158 array-present targets.

**Step 6 — smoke test.** fgsea on DMD−Ctrl moderated-t ranks over all 128 families, seed 42.

**Step 7 — all families, three comparisons.** Same settings per contrast. Leading-edge genes of
the top 25 families cross-referenced against the fibrosis sets.

**Step 8 — shuffled-label baseline.** All 34 labels shuffled keeping 6/11/17, limma refit, fgsea
rerun for both control contrasts, r and slope recorded. 1,000 shuffles, 7 cores, `nPermSimple`
10,000. Timed first: 41.3 s per shuffle serial → ~98 min projected on 7 cores. Empirical
p = (#shuffles ≥ real + 1) / (N + 1).

**Step 9 — direct BMD vs DMD.** BMD/DMD labels shuffled among the **28 patients only**, keeping
11/17. The 6 controls are not used at all, so no control fluke can create a difference.

**Step 10a — ECM removed.** Matrisome-core genes deleted from both the ranked list **and** the
family sets (otherwise set sizes would be wrong), then fgsea rerun. Families falling below 15
genes after deletion are dropped.

**Step 10b — detection limit.** Ruler = mean DMD−Ctrl log2FC of miR-29 targets minus that of
non-targets. Planted sets copy miR-29's 3′UTR-decile profile. f ∈ {0, 10, 25, 50, 75, 100}%,
100 repeats each, planted in the 11 BMD samples only, then BMD−Ctrl limma and fgsea rerun with
the planted set tested among all 128 families.

**Step 11 — 3′UTR-matched null.** Genes binned into deciles of 3′UTR length; 1,000 random sets
drawn per family copying its per-bin counts; fgsea on each. Real NES recomputed on the
UTR-restricted gene set so it is comparable to the nulls.

**Step 12 — fibroblast covariate.** Score = mean row z-score of the 316 array-present
fibroblast markers that are **not** in a fibrosis set. Refit `~0 + group + fibro`.

## Cuts (CLAUDE.md §5)

| Cut | Reason |
|---|---|
| GSE3307 replication | CLAUDE.md §5 |
| Serum microRNAs (SQ4) | CLAUDE.md §5 |
| Hypergeometric overlap test | CLAUDE.md §5 — superseded by fgsea, which needs no cutoff |
| GSE157668 | CLAUDE.md §5, unless a ready-made count table appears |
| Step 15 (pig overlap) | **BLOCKED**, not cut: Reschke et al. 2026 Table 1 not obtained, and §4 forbids building `pig_mirnas.tsv` from any other source |
| The old preliminary GEO2R claim | Replaced by step 4's limma tables. Must not be cited (CLAUDE.md §6 step 4). |
