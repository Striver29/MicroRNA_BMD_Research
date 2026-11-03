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
| E5 | **Step 10b ran on a different machine**: Ubuntu 26.04, AMD EPYC 9R14, 16 physical cores (1 thread/core), 30 GB RAM, 3 Nov 2026. R 4.6.2 from CRAN's Ubuntu apt repo; packages as Posit Package Manager binaries, Bioconductor 3.23. **Every relevant version is identical to the Mac runs:** limma 3.68.5, fgsea 1.38.0, BiocParallel 1.46.0, data.table 1.18.6.1, ggplot2 4.0.3, scales 1.4.0. Only these six packages were installed; `00_setup.R` was not run. | The Mac (4 cores) needed ~15.5 CPU-h for 10b and kept sleeping. Inputs were the committed `gene_matrix.rds`, `family_sets.rds`, `limma.rds`, `utr_lengths.tsv`, so no CEL download or RMA was needed. |

## Decisions and defaults used (CLAUDE.md §10)

| Decision | Value used | Note |
|---|---|---|
| Probe collapse | highest mean expression | §10 default |
| Low-expression filter | none | §10 default; ranks look healthy (t range −15.0 to +17.7) |
| TargetScan target set | 8.0 default (conserved) predictions, human 9606 | §10 default |
| Muscle-expression source | miRNA Tissue Atlas 2025 (CCB Saarland), `Organ_system = limb_muscle`, mean RPMM ≥ 10 | §4 let us pick the first credible published atlas; see `docs/handoffs.md`. The **RPMM ≥ 10 threshold is ours**, not the spec's — logged here as a parameter to sensitivity-test if asked. |
| Family size limits | ≥ 15 targets after intersecting with the gene matrix, no maximum | §10 default |
| Smoke-test thresholds | GO = NES > 0 and padj < 0.05 | §6 step 6 |
| Shuffle count | **500** (steps 8 and 9); 1,000 draws (step 11); 100 reps (steps 10b, 14) | §10 permits 500 "if too slow". Measured throughput after the parallelism fix is **20.4 s per shuffle** on 7 workers (21-shuffle parallel batch, not the script's serial extrapolation), so 1,000 shuffles = 5.7 h for step 8 alone and ~11 h for steps 8+9 together. Dropped to 500, giving a minimum empirical p of 1/501 = 0.002 — ample resolution. **Logged per §10.** |
| fgsea permutations | `nPermSimple = 10000`, `eps = 0` for every reported run | At the default 1,000, 13 of 128 families returned `NA` NES (fgsea's own warning). Raised to 10,000 to remove the NAs. Inside step 11's null draws only, 1,000 is used, because there we need the NES distribution and not a precise p per draw. |
| Catch criterion (step 10b) | NES > 0 and padj < 0.05 | §10 default. **Failed its own check: false-alarm rate 1.00 at f = 0.** A second, calibrated rule (NES above the 95th percentile of the 100 unplanted runs) is reported alongside by `10b2_calibrated_catch.R`. Which one is the headline is **Aziz's call** — not switched. |
| Random seeds | 1 (steps 1–4), 42 (fgsea steps 6, 7, 10a, 12), 20261101 (step 8), 20261102 (step 9), 20261103 (step 11), 20261104 (step 10b), 20261105 (step 10b2) | every script sets its own |

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

**DECISION (Aziz, 2026-11-02): keep all 34 samples.** STOP point of CLAUDE.md §6 step 2 is
closed.

Basis for the decision:
- **RLE and NUSE flag nothing.** These are the measures that actually test array quality, and
  they are computed within a gene's own ~11 probes, so patient biology cannot inflate them.
  All 34 arrays sit well inside the conventional thresholds (RLE ±0.057 vs ±0.10 allowed;
  NUSE 0.981–1.039 vs 1.05 allowed).
- **The spec's original metric was confounded and is not used for the decision.** "Median
  correlation with all other samples" measures how many samples share your biology, not array
  quality. With 28 dystrophic arrays and only 6 healthy ones, all 6 controls were guaranteed to
  score lowest, and they did (lowest 6 of 34). Worked example, GSM2934814 (a control): it
  correlates with the other 5 controls at median r = 0.9854 — i.e. its array is excellent — but
  28 of its 33 comparisons are against diseased muscle (median r 0.9437 vs DMD, 0.9620 vs BMD),
  so its median-over-all collapses to 0.9496. Following that metric would have deleted controls
  for being healthy, in a 6-control study.
- **The within-group metric does still flag two samples, mildly**: GSM2934803 (BMD, 0.9676 vs a
  0.9757 group centre) and GSM2934813 (control, 0.9784 vs 0.9839). Both are kept. With RLE and
  NUSE clean these are patient-to-patient biological spread, not broken arrays, and dropping
  them would delete real variation in a disease known for a variable course.
- **The residual risk is tested, not assumed.** If either sample were carrying the miR-29
  result, step 13's leave-one-out (17 re-runs, dropping each BMD and each control in turn) will
  expose it. That is better evidence than a decision made now on a metric already shown to be
  unreliable here.

Reopen this only if step 13 shows miR-29 collapsing on the removal of one specific sample.

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
the planted set tested among all 128 families. Run with `N_REP=100 NPERM=10000 N_CORE=16`
(full spec, no reduction).

**Step 10b2 — calibrated catch rule.** Re-scores 10b's 600 saved runs; no limma or fgsea rerun.
Catch = planted-set NES > 95th percentile of the NES of the 100 unplanted (f = 0) runs (bar =
1.812). The f = 0 false-alarm rate is therefore 5% **by construction** (the bar is set on the
same runs it is tested on), not an independent estimate. The 80% crossing is a straight-line
interpolation between grid points 10% and 25%. Also computes how much of the ruler any
UTR-matched random set already shows (1,000 draws per contrast, same decile profile as 10b).
**What this rule still assumes, all in our favour:** the plant is the same shift in every BMD
patient; planted genes are random, so they don't move together across patients; and the null
varies the gene choice, not the patient labels. It is a best-case limit.

**Step 11 — 3′UTR-matched null.** Genes binned into deciles of 3′UTR length; 1,000 random sets
drawn per family copying its per-bin counts; fgsea on each. Real NES recomputed on the
UTR-restricted gene set so it is comparable to the nulls.

**Step 12 — fibroblast covariate.** Score = mean row z-score of the 316 array-present
fibroblast markers that are **not** in a fibrosis set. Refit `~0 + group + fibro`.

## Performance bug found and fixed (session 1, must not recur)

**Symptom.** Step 8 (1,000 shuffles) was projected at 98 minutes. After 99 minutes of wall
clock it had completed roughly **19%** of the work. Extrapolated finish: **~6.8 hours**.

**Diagnosis.** `fgsea::fgseaMultilevel` defaults to `nproc = 0`, which means "use every core".
Our shuffle loops call `fgsea` *inside* a `BiocParallel::bplapply` with 7 workers. So each of
the 7 workers spawned its own full-width parallel fgsea: 7 x 8 = up to 56 threads competing
for 8 cores.

**Evidence it was thrashing, not computing:**
- 20 live processes for a job configured with 7 workers;
- `top` reported **48% user / 52% system** CPU — more than half the machine was spent on
  kernel context-switching rather than arithmetic;
- the 7 worker processes had each accumulated only ~18.5 minutes of CPU time across ~94
  minutes of wall clock, i.e. each worker was idle ~80% of the time.

**Fix.** Pass `nproc = 1` to every `fgsea()` call that executes inside a parallel worker, so
the outer `bplapply` owns the parallelism and the inner fgsea stays serial. Applied to
`08_shuffle_baseline.R`, `09_bmd_vs_dmd.R`, `10b_detection_limit.R`, `11_utr_matched.R`,
`13_leave_one_out.R`, `14_subsample_dmd.R`.

**Deliberately NOT applied** to `06_smoke_test.R`, `07_fgsea_all.R`, `10a_ecm_removed.R`,
`12_fibroblast_cov.R`: these call fgsea at the top level with nothing wrapped around it, so
there full-width parallelism is correct and faster.

**Timing note for whoever re-runs step 8.** The original "41.3 s per iteration" figure was
measured in the script's serial timing loop, where fgsea silently used all 8 cores. With
`nproc = 1` the honest single-core cost is **67.0 s per iteration**. The projection printed by
the script (`N x per / cores`) is therefore ~160 min for 1,000 shuffles, but that figure has
not been validated against real parallel throughput — **measure a small parallel batch (e.g.
21 shuffles on 7 workers) before committing to the full run.** Step 8 was stopped unfinished
at the end of session 1 and must be re-run from the start.

## Measured runtimes (session 2, after the nproc fix)

Measured with a 21-shuffle parallel batch on 7 workers, which is the only trustworthy way to
project these jobs — the scripts' own `N x serial_per_iter / cores` estimate ignores the ~53%
efficiency loss from memory contention and has misled us once already.

| Quantity | Value |
|---|---|
| step 8, serial per shuffle (`nproc = 1`) | 67.0 s |
| step 8, **wall per shuffle on 7 workers** | **20.4 s** |
| implied parallel efficiency | ~47% (ideal would be 9.6 s) |
| step 8, 1,000 shuffles | 5.7 h → **rejected** |
| step 8, 500 shuffles | **2.8 h → chosen** |
| step 13, 17 leave-one-out runs x 2 contrasts | 7.3 min (measured) |
| step 10b, serial per simulation (EPYC, `nproc = 1`) | 102.0 s |
| step 10b, **32-run parallel batch on 16 workers** | 3.4 min = 2 waves × ~102 s → ~100% efficiency; CPU 100% user, 0% system |
| step 10b, full 600 runs on 16 workers | projected 64 min, **measured 76.2 min** |

## Cuts (CLAUDE.md §5)

| Cut | Reason |
|---|---|
| GSE3307 replication | CLAUDE.md §5 |
| Serum microRNAs (SQ4) | CLAUDE.md §5 |
| Hypergeometric overlap test | CLAUDE.md §5 — superseded by fgsea, which needs no cutoff |
| GSE157668 | CLAUDE.md §5, unless a ready-made count table appears |
| Step 15 (pig overlap) | **BLOCKED**, not cut: Reschke et al. 2026 Table 1 not obtained, and §4 forbids building `pig_mirnas.tsv` from any other source |
| The old preliminary GEO2R claim | Replaced by step 4's limma tables. Must not be cited (CLAUDE.md §6 step 4). |
