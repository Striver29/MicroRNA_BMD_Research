# Results log

Every number Aziz might quote, with the script and the output file that produced it.
Session 1: 2026-11-01. **Nothing here is from a step that did not actually run.**

---

## Step 1 — labels (`scripts/01_download_labels.R` → `results/tables/01_sample_sheet.tsv`)

| Number | Value |
|---|---|
| Samples in GSE109178 | 49, all GPL570, all vastus lateralis |
| Group counts | **6 control, 11 BMD, 17 DMD, 15 LGMD** — matches the spec exactly |
| LGMD split | 8 dysferlinopathy (LGMD2B) + 7 FKRP (LGMD2I) |
| Included in the analysis | 34 |
| BMD ages | median **12.2 y**, range **2.5–60.8**, all 11 known |
| DMD ages | median **3.5 y**, range **0.9–8.0**, all 17 known |
| **Control ages** | **not reported in GEO for any of the 6** `[UNVERIFIED]` |
| BMD severity | 6 mild, 3 moderate, 2 severe |
| DMD severity | 4 moderate, 13 severe |
| GSM2934804 | GEO says `severe`, `BMD-12-[6.4]-severe`, Deletion Exon 49, age 6.4 |

**Carry into Limitations:** the BMD group is about a decade older than the DMD group, and
control ages are unknown. The direct BMD-vs-DMD comparison therefore confounds disease with
age. Stated, not fixed (CLAUDE.md §2).

## Step 2 — RMA and QC (`scripts/02_rma_qc.R` → `results/tables/02_qc_metrics.tsv`)

| Number | Value |
|---|---|
| RMA matrix | **54,675 probesets × 34 samples** |
| Median RLE, all 34 arrays | within **±0.057** (flag threshold 0.10) → none flagged |
| Median NUSE, all 34 arrays | **0.981–1.039** (flag threshold 1.05) → none flagged |
| Within-group correlation flags | GSM2934803 (BMD, 0.9676 vs 0.9757 group centre), GSM2934813 (control, 0.9784 vs 0.9839) |
| Overall-correlation metric by group | control 0.9529, BMD 0.9704, DMD 0.9803 |
| Within-group metric by group | control **0.9839**, BMD 0.9757, DMD 0.9855 |

**✅ STOP point CLOSED — Aziz, 2026-11-02: keep all 34 samples.** RLE and NUSE, the measures
that actually test array quality, flag nothing. The spec's "median correlation with all other
samples" is confounded by group size here (28 of 34 arrays are dystrophic, only 6 healthy), so
all 6 controls scored lowest by construction; on the within-group metric those same controls
score *highest*. Two samples remain mildly flagged within group (GSM2934803, GSM2934813) and
are kept; step 13's leave-one-out tests whether either is carrying the result.
Full reasoning in `docs/methods_log.md`.

## Step 3 — probes to genes (`scripts/03_annotate_collapse.R`)

54,675 probes → 44,647 with a symbol → 1,524 multi-mapping dropped → 43,123 kept →
**20,823 genes**. 11,043 genes had more than one probe; the highest-mean probe won.

## Step 4 — limma (`scripts/04_limma.R` → `results/tables/limma_*.tsv`)

Residual df 31, eBayes prior df 4.23.

| Contrast | FDR < 0.05 up | FDR < 0.05 down |
|---|---|---|
| BMD − control | 4,549 | 4,659 |
| DMD − control | 5,547 | 6,654 |
| BMD − DMD | 4,564 | 4,039 |

**Fibrosis sanity table** (`results/tables/04_sanity_fibrosis_genes.tsv`), log2 fold change:

| Gene | BMD−Ctrl | adj.P | DMD−Ctrl | adj.P | BMD−DMD | adj.P |
|---|---|---|---|---|---|---|
| COL1A1 | +2.517 | 9.4e−08 | **+4.058** | 5.6e−13 | −1.542 | 4.6e−06 |
| COL1A2 | +2.043 | 1.4e−07 | +2.940 | 7.5e−12 | −0.897 | 4.2e−04 |
| COL3A1 | +2.311 | 8.3e−09 | +3.260 | 4.6e−13 | −0.949 | 1.3e−04 |
| COL5A1 | +1.919 | 1.4e−06 | +2.828 | 7.5e−11 | −0.910 | 8.5e−04 |
| LOX | +2.607 | 2.4e−07 | +3.573 | 3.7e−11 | −0.966 | 3.1e−03 |

All five up in **both** diseases, and all five **less** up in BMD than in DMD. On the linear
scale COL1A1 is 2^4.058 ≈ **16×** in DMD and 2^2.517 ≈ **5.7×** in BMD.
Consistent with Eisenberg 2007 (COL3A1 ~6×, COL1A2 ~5× in DMD): same direction, same order of
magnitude. This replaces the old preliminary GEO2R claim, which must not be cited.

## Step 5 — microRNA sets (`scripts/05_mirna_inputs.R` → `results/tables/05_family_sizes.tsv`)

2,064 human TargetScan 8.0 families → 253 muscle-expressed → **128 families** with ≥ 15
array-present targets. Sizes 15–1,446, median 371.
miR-29-3p: seed `AGCACCA`, members hsa-miR-29a/b/c-3p, conservation 2, muscle 1,895 RPMM,
**1,158 array-present targets**.

## Step 6 — smoke test (`scripts/06_smoke_test.R` → `results/tables/06_smoke_fgsea_DMDvsCtrl.tsv`)

### ✅ VERDICT: **GO**

| | |
|---|---|
| miR-29-3p, DMD vs control | **NES = +2.224, p = 3.4e−36, padj = 5.4e−35** |
| rank among the 128 families | **3rd** |
| direction | NES > 0 = targets up = brake **weaker** = miR-29 predicted **less active** in DMD — the published direction |
| mean log2FC of the 1,158 targets | **+0.304** |
| mean log2FC of the other 19,665 genes | **+0.036** |
| **the ruler** (footprint size, used by step 10b) | **+0.269 log2 units** ≈ 20.5% on the linear scale |

**⚠ The specificity caveat, which must travel with this result.**
103 of 128 families reached padj < 0.05 in this comparison, and **99 of them point the same way
(NES > 0)**. A competitive rank test on overlapping target sets does that when the whole
transcriptome shifts. So the smoke test shows the method **can** recover known biology; it does
**not** on its own show that miR-29 is special. The specificity evidence is steps 9, 10a, 11, 12.

## Step 7 — all families, three comparisons (`scripts/07_fgsea_all.R`)

| Contrast | padj < 0.05 | of which NES > 0 | NES < 0 |
|---|---|---|---|
| BMD − control | 104 / 128 | 99 | 5 |
| DMD − control | 103 / 128 | 99 | 4 |
| **BMD − DMD** | **7 / 128** | **0** | **7** |

miR-29-3p across the three:

| Contrast | NES | padj | rank |
|---|---|---|---|
| BMD − control | +2.169 | 1.4e−31 | 6 / 128 |
| DMD − control | +2.224 | 1.5e−34 | 3 / 128 |
| **BMD − DMD** | **−1.590** | **4.7e−10** | **1 / 128** |

**This is the headline shape.** When the shared controls are removed, the global skew collapses
(103 significant → 7) and **miR-29 is rank 1, about five orders of magnitude ahead of the next
family** (miR-30-5p, 1.9e−04; miR-1-3p/206, 4.1e−04).
NES < 0 for BMD−DMD means miR-29's targets are *lower* in BMD than in DMD → the miR-29 brake is
**less weakened in BMD than in DMD** → milder, in the predicted direction.

**miR-29 leading-edge fibrosis genes, DMD vs control** (`results/tables/07_leading_edge.tsv`):
COL15A1, COL16A1, COL1A1, COL1A2, COL21A1, COL27A1, COL3A1, COL4A1, COL4A2, COL4A4, COL5A1,
COL5A2, COL5A3, COL6A2, COL6A3, COL6A6, CRISPLD1, ELN, FBN1, FRAS1, HMCN1, LAMA2, LAMC1, MFAP2,
MFAP3, NID1, PXDN, SKI, SPARC, SVEP1.

miR-29's leading edge is **6.5% fibrosis-set genes — the highest of any top family** (next
highest 4.1%, most are 1–3%). That is a specificity signal in the content, independent of the p-value.

| BMD/DMD footprint agreement across 128 families | value |
|---|---|
| Pearson r (BMD NES vs DMD NES) | **0.9123** |
| slope (BMD NES on DMD NES) | **0.9769** |

**Do not quote these two against 0.** They must be read against the step-8 shuffled baseline.

## Step 10a — ECM removed (`scripts/10a_ecm_removed.R` → `results/tables/10a_ecm_removed.tsv`)

Dropping all 275 matrisome-core genes from both the ranked list and the family sets:

| | NES before | NES after | padj after | survives |
|---|---|---|---|---|
| miR-29, DMD − control | +2.224 | **+2.151** (96.7% retained) | 1.4e−29 | ✅ |
| miR-29, BMD − control | +2.169 | **+2.109** | 1.0e−27 | ✅ |

Across all families the median NES was essentially unchanged (1.679 → 1.684).
**The footprint is not just the scar genes.**

## Step 12 — fibroblast covariate (`scripts/12_fibroblast_cov.R`)

Score from 316 array-present fibroblast markers that are **not** fibrosis-set genes.

| Group | mean score |
|---|---|
| control | **−0.791** |
| BMD | **+0.040** |
| DMD | **+0.253** |

ANOVA on group p = 1.2e−08. BMD sits between healthy and DMD, as expected.

| Contrast | sig before | sig after | survive |
|---|---|---|---|
| BMD − control | 104 | 80 | 74 |
| DMD − control | 103 | **51** | 45 |
| BMD − DMD | 7 | 3 | 3 |

miR-29 survives in all three:

| Contrast | NES raw → adjusted | padj adjusted |
|---|---|---|
| BMD − control | +2.169 → +1.596 | 1.1e−10 |
| DMD − control | +2.224 → **+1.762** (79% retained) | **1.7e−16** |
| BMD − DMD | −1.590 → −1.300 | 7.1e−03 |

**Second specificity signal:** adjusting for fibroblast content halves the number of significant
families (103 → 51) but miR-29 retains more of its NES than the other top families (79% vs
~55–60%), and afterwards it leads by eleven orders of magnitude (1.7e−16 vs 3.1e−05 next).
Caveat per CLAUDE.md §12: this covariate also absorbs real disease signal, so read it as
"does the footprint survive", not as the headline estimate.

---

## Step 11 — 3′UTR-length-matched null (`scripts/11_utr_matched.R`)

**This is the most sobering result so far, and it must be in the paper.**

Genes with long 3′UTRs appear on many target lists *and* behave alike in disease, so a plain
random comparison is biased in our favour. We drew 1,000 random gene sets per family that copy
that family's 3′UTR-length profile decile by decile, and scored them identically.

**A UTR-matched random set already scores NES ≈ +1.6.** Not zero.

| | real NES | UTR-matched null mean | null SD | z | empirical p |
|---|---|---|---|---|---|
| miR-29, DMD vs control | 2.222 | **1.592** | 0.096 | **6.56** | 0.000999 |
| miR-29, BMD vs control | 2.172 | **1.669** | 0.093 | **5.39** | 0.000999 |

*(p = 0.000999 is the floor, 1/1001 — no random set out of 1,000 ever reached the real value.)*

**How to read this.** About **72%** of miR-29's apparent enrichment (1.59 of 2.22) is explained
by 3′UTR length alone. The footprint is still real — the excess sits **6.5 standard deviations**
above a properly matched baseline, and not one of 1,000 matched sets came close — but the raw
NES of 2.22 **badly overstates** the effect. Never quote NES as if zero were the baseline.

**All 25 top families survived at q < 0.05**, so this step validates the method broadly rather
than singling out miR-29. By z-score miR-29 ranks 2nd in DMD vs control (behind miR-1/206).

## Step 13 — leave-one-out (`scripts/13_leave_one_out.R`) — 7.3 min

Dropped each of the 11 BMD and 6 control samples in turn, 17 full re-runs.

| | full data | range across 17 re-runs | significant in |
|---|---|---|---|
| miR-29, BMD vs control | +2.169 | +2.120 to +2.186 | **17/17** |
| miR-29, BMD vs DMD | −1.590 | −1.674 to −1.500 | **17/17** |

miR-29 stayed **rank 1 of 128 in every single BMD-vs-DMD run**. Top-20 family list overlapped
the full-data list by 18–20 of 20 in every run.

**This closes the step-2 QC question.** The two within-group-flagged samples are harmless:
- **GSM2934803** (BMD): largest effect of any sample, and still 18/20 overlap; miR-29 rank 1.
- **GSM2934813** (control): 19/20 overlap.

**GSM2934804** (disputed severity, CLAUDE.md §2 robustness item) is also immaterial: dropping it
gives NES +2.156 (vs +2.169) and −1.507 (vs −1.590), rank 1 either way. **Item closed.**

## Step 14 — DMD subsampled to n = 11 (`scripts/14_subsample_dmd.R`) — 31.1 min

**⚠ A negative result that changes how the headline must be phrased.**

| | miR-29 NES vs control | families significant |
|---|---|---|
| full DMD (n = 17) | 2.224 | 103 / 128 |
| **DMD cut to 11**, 100 draws | **median 2.167** (2.074–2.258) | median 102 (98–107) |
| **actual BMD (n = 11)** | **2.169** | 104 / 128 |

**53.0%** of the DMD-at-n=11 draws land at or below BMD's value — a coin flip. BMD sits in the
middle of the distribution DMD itself produces at matched sample size.

**So: measured against healthy controls, BMD's miR-29 footprint is NOT detectably weaker than
DMD's.** The apparent gap (2.169 vs 2.224) is entirely within sample-size noise.

**Why this does not contradict the direct comparison.** NES is a *rank* measure and it
saturates: it asks "are the targets bunched at the up end?", not "how far did they move?" In
both diseases miR-29's targets sit firmly at the up end, so both score ≈ 2.2. The magnitudes
differ greatly — COL1A1 is **+2.52 in BMD vs +4.06 in DMD** — and only the direct BMD-vs-DMD
contrast can see that.

**Required phrasing for the paper.** Do **not** write "BMD's miR-29 footprint is weaker than
DMD's" on the basis of the vs-control scores. Write:

> BMD shows the **same miR-29 footprint pattern** as DMD. The **direct** BMD-vs-DMD comparison
> shows miR-29 targets are significantly **less elevated** in BMD, consistent with milder
> disease; but the enrichment strength measured against healthy controls is indistinguishable
> between the two groups at matched sample size.

## Step 10b — detection limit (`scripts/10b_detection_limit.R` + `scripts/10b2_calibrated_catch.R`) — 76.2 min on 16 cores

Outputs: `results/tables/10b_detection_curve.tsv`, `10b_detection_runs.tsv` (all 600 runs),
`10b_log.txt`; `10b2_detection_curve_calibrated.tsv`, `10b2_ruler_utr_share.tsv`, `10b2_log.txt`;
figures `10b_detection_curve.png`, `10b2_detection_curve_calibrated.png`.

**⚠ The spec's catch rule failed its own sanity check. The false-alarm rate is 1.00, not ≤ 0.05.**

| | Value |
|---|---|
| Ruler (DMD miR-29 footprint) | **+0.2686 log2** ≈ 20.5% linear — matches step 6 |
| Planted set size | 1,158 genes, copying miR-29's 3′UTR-decile profile |
| Unplanted (f = 0) NES, 100 runs | mean **1.679**, SD 0.087, range 1.455–1.867 — matches step 11's UTR-matched null (1.669) |
| Spec rule at f = 0 | **100 / 100 caught**; even the *least* significant unplanted set had padj **2.1e−06** |
| Calibrated bar (95th pct of unplanted NES) | **1.812** |

| f (% of DMD footprint) | planted log2 | median NES | z vs unplanted | catch, spec rule | **catch, calibrated** |
|---|---|---|---|---|---|
| 0 | 0.000 | 1.67 | 0.00 | 1.00 | **0.05** (by construction) |
| 10 | 0.027 | 1.87 | 2.28 | 1.00 | **0.73** |
| 25 | 0.067 | 2.17 | 5.59 | 1.00 | **1.00** |
| 50 | 0.134 | 2.59 | 10.52 | 1.00 | 1.00 |
| 75 | 0.201 | 2.94 | 14.63 | 1.00 | 1.00 |
| 100 | 0.269 | 3.23 | 17.81 | 1.00 | 1.00 |

- **Spec rule:** the curve is flat at 100%, so it carries no information. The script's own
  printout "80% detection reached at f = 0%" is an artefact of this. **Do not quote it.**
- **Calibrated rule:** 80% is crossed between f = 10% and 25%, at **≈ 14% (+0.037 log2)** by
  straight-line interpolation. The grid is coarse there.

**How much of the ruler is just 3′UTR length** (`10b2_ruler_utr_share.tsv`, 1,000 UTR-matched draws):

| Contrast | miR-29 gap (targets − others) | UTR-matched random set gap | share explained by UTR length |
|---|---|---|---|
| DMD − control | +0.2686 (the ruler) | +0.1124 | **42%** |
| BMD − control | +0.1923 (= 72% of the ruler) | +0.0973 | **51%** |

So BMD's miR-29-specific excess is 0.1923 − 0.0973 = **+0.095 log2 ≈ 35% of the ruler**. On the
planting scale that sits at f ≈ 35%, and BMD's real miR-29 NES (2.169) does match the planted
median at f = 25–35% (2.17 at f = 25%). Cross-check: planted z at f = 25% is 5.59; step 11's
UTR-matched z for real BMD miR-29 was 5.39.

**How to read this (the important part).**
1. Even under the best case (the same shift in every patient, uncorrelated genes,
   gene-sampling null), 10b2 says BMD's real miR-29 footprint (f ≈ 35%) is **well above** the
   detectable size (f ≈ 14%). It would be caught 100% of the time.
2. Yet the honest patient-label shuffle (step 8) gives miR-29 BMD-vs-control **p = 0.084**.
3. So the binding limit at n = 11 is **not** the size of the gene-level signal. It is
   **patient-to-patient variability**, which a constant plant on random genes doesn't
   simulate. **10b as specified can't measure the limit that matters for our null result.**
   The 10b2 curve is a best-case floor: the real detection limit is at least this large,
   probably much larger.
4. The spec rule's 100% false-alarm rate is a third, concrete example of the fgsea
   overconfidence finding. Any random set of 1,158 long-UTR genes gets padj ≤ 2.1e−06.

**Do NOT write:** "our method detects footprints down to 14% of DMD's at n = 11." Write
instead: "Under idealised planting, a footprint ≈ 14% of the DMD miR-29 footprint was
detectable against a 3′UTR-matched null. The observed BMD miR-29 footprint exceeds this, yet
does not survive patient-label permutation (p = 0.084). At n = 11 the detection limit is set by
between-patient variability, which planting simulations of this kind do not capture."

## Still running / not yet run

| Step | Status |
|---|---|
| 8 — shuffled-label baseline, 500 shuffles | done (`08_family_shuffle_p.tsv`); headline numbers in the CLAUDE.md status block, full section not yet written here |
| 9 — direct BMD vs DMD, 500 patient-label shuffles | done (`09_bmd_vs_dmd.tsv`); as above |
| 10b — detection limit | **done** (above) |
| 15 — pig overlap | **BLOCKED**: Reschke et al. 2026 Table 1 not obtained |
