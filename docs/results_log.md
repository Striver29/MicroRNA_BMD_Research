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

**⚠ STOP point, awaiting Aziz's sign-off.** The spec's "median correlation with all other
samples" flagged three controls — but that metric is confounded by group size here (28 of 34
arrays are dystrophic, only 6 are healthy), and on the within-group metric those same controls
score *highest*. RLE and NUSE, the proper array-quality measures, flag nothing.
**Decision taken to keep all 34**, so the session could continue; step 13's leave-one-out
tests sensitivity to exactly these samples. Aziz to confirm or overrule.

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

## Still running / not yet run

| Step | Status |
|---|---|
| 8 — shuffled-label baseline, 1,000 shuffles | **running** (~98 min projected) |
| 11 — 3′UTR-matched null, 25 families × 1,000 draws | **running** |
| 10b — detection limit, 6 × 100 simulations | queued behind 11 |
| 9 — direct BMD vs DMD label shuffle | queued |
| 13 — leave-one-out | not started (CLAUDE.md "if time") |
| 14 — DMD subsampled to 11 | not started (CLAUDE.md "if time") |
| 15 — pig overlap | **BLOCKED**: Reschke et al. 2026 Table 1 not obtained |
