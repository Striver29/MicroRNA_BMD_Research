# Write-up plan: what to report, and where it goes

**Compiled 7 November 2026, after the computation freeze.** First draft due Sunday 8 November.

This file is the drafting map. Every number here traces to a script that ran; sources are named
in square brackets. If a number is not in this file, check `docs/results_log.md` before using it
— and if it is in neither, do not write it.

**Three documents, three jobs:**

| File | Job |
|---|---|
| `docs/FULL_RESEARCH_DOCUMENT.md` | The full reference. Why each step exists, worked examples, teaching. |
| **`docs/WRITEUP_PLAN.md` (this file)** | **What goes in the paper, in which section, with the numbers.** |
| `docs/results_log.md` | Provenance for every number. |

> ### ⚠ Four framing corrections this file applies to the reference document
> These are not stylistic. Get any of them wrong and the section they sit in becomes
> indefensible.
>
> 1. **The methods finding is a confirmation, not a discovery.** Cite the prior art (§6.3).
> 2. **Our positive evidence uses the same kind of null we criticise.** Say so first, ourselves (§6.4).
> 3. **The collagen result is age-confounded in the direction of the finding.** Caveat travels with it (§5.2, §7.3).
> 4. **The duplicate-accession claim is scoped to one verified array pair.** Do not generalise (§5.5).

---

# 1. Title and framing

**Title:** *Can microRNA target footprints recover known microRNA biology in dystrophic muscle?
Detection limits at n = 11.*

**One-line verdict:** We asked whether the miR-29–fibrosis axis known in Duchenne also operates
in Becker, built a method that provably recovers the known Duchenne biology, then tested it with
the correct null — and 11 Becker patients is not enough. We state quantitatively how not-enough.

**Narrative order (do not reverse it):** lead with the biological question, let the methods
finding emerge as *the reason the biological answer is null*. The methods result is stronger and
more certain, but leading with it makes the paper look like it is hiding a failure.

**Claim language throughout:** "predicted", "candidate", "may contribute", "consistent with".
Never "causes", "proves", "demonstrates", "validated biomarker". State explicitly that microRNA
activity is **inferred from target behaviour and never measured** — the array cannot see
microRNAs. With 11 BMD patients every result is hypothesis-generating; say so.

---

# 2. Abstract — draft to adapt

> Becker muscular dystrophy (BMD) is the milder, in-frame form of Duchenne (DMD), and both
> progress through fibrotic replacement of muscle. In DMD, miR-29 is reduced and de-represses
> collagen genes, but whether the same axis operates in human BMD is untested — a question made
> pressing by exon-skipping therapies that aim to convert DMD into a BMD-like state. Using 6
> control, 11 BMD and 17 DMD vastus lateralis biopsies (GSE109178, Affymetrix HG-U133 Plus 2.0),
> we inferred microRNA activity from the collective behaviour of predicted target genes across
> 128 muscle-expressed seed families. The approach recovered the expected biology directionally:
> miR-29 targets were coordinately elevated in DMD muscle (NES +2.22, 2nd of 128 families), its
> leading edge was the most fibrosis-dense of any family (6.5%), and the signal persisted after
> removing all 275 core-matrisome genes (96.7% retained), adjusting for fibroblast content (79%
> retained), matching on 3'UTR length (6.5 SD above a matched null) and leave-one-out resampling
> (17/17 runs). All of these, however, assess significance by permuting genes. When significance
> was instead assessed by permuting **patient labels** — the only null that tests whether disease
> status drives the footprint — no result survived: miR-29 reached p = 0.052 in DMD and p = 0.084
> in BMD, no family differed between the two diseases (0 of 128, all q > 0.95), and the apparent
> agreement between their footprints (r = 0.91) fell inside the shuffled-label distribution
> (p = 0.064). In a controlled test where the ground truth was known to be nothing, unmodified
> random gene sets were called significant in 100 of 100 trials (weakest padj 2.1e-06),
> confirming in this setting a known hazard of competitive gene-set tests on large, correlated
> sets. Simulation indicates the binding limit at n = 11 is between-patient variability rather
> than the magnitude of the gene-level signal. We report a null result for the Becker comparison
> with a measured detection floor, alongside one robust per-gene finding: collagen genes are
> elevated in both diseases and significantly less so in BMD.

**Numbers the abstract commits you to:** 6/11/17; 128 families; NES +2.22; 2nd of 128; 6.5%;
96.7%; 79%; 6.5 SD; 17/17; p = 0.052; p = 0.084; 0 of 128; q > 0.95; r = 0.91; p = 0.064;
100 of 100; padj 2.1e-06. All verified.

---

# 3. Introduction — what goes in, in order

1. **Dystrophin, and frame rule.** DMD = out-of-frame, essentially no dystrophin, onset 2–5 y.
   BMD = in-frame, partly functional protein, milder and hugely variable. Same gene, two
   severities.
2. **Fibrosis is a driver, not a side effect.** Dead fibres are replaced by collagen scar;
   muscle is literally replaced by non-contractile tissue.
3. **MicroRNAs as brakes.** ~22 nt, seed = positions 2–8, binds 3'UTR, transcript degraded or
   blocked. One microRNA tunes hundreds of genes. Families share a seed and are treated as a unit.
4. **miR-29 is the brake on collagen, and it is down in DMD.** Cite Eisenberg et al. 2007 PNAS:
   miR-29c about 6-fold down in DMD, COL3A1 about 6-fold up, COL1A2 about 5-fold up. *Write
   "about" — the units are not stated in the source.*
5. **Why BMD is open.** The same paper had only **6 BMD samples**; in BMD only miR-146b and
   miR-221 passed their fold-change cutoff, and their microRNA–gene linkage reached significance
   only in DMD and Miyoshi myopathy — **not BMD**. Nobody has shown the axis operates in human BMD.
6. **Why it matters now.** Exon-skipping therapy aims to convert DMD into a BMD-like state.
   Naturally occurring BMD is the closest available proxy for a successfully treated DMD patient
   years later.
7. **The inference problem.** Our platform measures mRNA, not microRNA. We infer activity from
   **footprints**: when a brake weakens, its targets drift up together. State the direction rule
   — NES > 0 means targets bunched at the up end, implying the microRNA is *less* active.
8. **The question, and the three possible answers:** weaker copy / different pattern /
   undetectable. All three were publishable. State that we landed on the third.

**Hamza owns the Intro and Related Work drafts.** He needs items 4, 5 and 6 sourced.

---

# 4. Methods — every subsection, with the numbers

Write in past tense. Each bullet below is roughly one sentence to one short paragraph.

## 4.1 Data and sample selection
- **GSE109178**, NCBI GEO. Affymetrix HG-U133 Plus 2.0 (**GPL570**), **54,675 probesets**.
- **Vastus lateralis** needle biopsies, all samples, same muscle. **49 samples total; 34 used.**
- Groups: **6 control, 11 BMD, 17 DMD**. **15 LGMD excluded** (8 LGMD2B/dysferlinopathy +
  7 LGMD2I/FKRP) as not dystrophin-related.
- **Ages:** BMD median **12.2 y** (range **2.5–60.8**), 6 mild / 3 moderate / 2 severe.
  DMD median **3.5 y** (range **0.9–8.0**), 4 moderate / 13 severe.
  **Control ages are not reported in GEO at all** — mark `[UNVERIFIED]` and carry it to Limitations.
- **Disclose the download route.** `GEOquery` would not build on the analysis machine and NCBI
  FTP returned HTTP 403, so metadata and the raw archive were retrieved from GEO's two official
  web endpoints (`acc.cgi?...form=text&view=brief` and `download/?acc=...&format=file`). Only
  the 34 included CEL files were extracted. Group counts matched the expected 6/11/17/15 exactly.
  [`01_download_labels.R`]

## 4.2 Normalisation and quality control
- **RMA** (background correction, quantile normalisation, probeset summarisation) on the
  **34 included arrays only**. State the reason: quantile normalisation lets every array
  influence every other, so including the 15 LGMD arrays would have let a different disease
  shape our values.
- **QC, and a metric we rejected.** The pre-specified metric (each array's median correlation
  with all others) is confounded by group size here: 28 of 34 arrays are dystrophic, so the 6
  controls are correlated against 28 biologically unlike arrays and only 5 like them, and are
  **guaranteed** to score lowest. On that metric three controls looked like outliers. Overall
  correlation by group: control **0.9529**, BMD 0.9704, DMD 0.9803. **Within-group** correlation
  reverses it: control **0.9839**, BMD 0.9757, DMD 0.9855 — the same controls now rank highest.
- **The metrics that actually test the chip** are computed within each probeset's own ~11 probes,
  so patient biology cannot inflate them: **median RLE within ±0.057** (flag threshold 0.10) and
  **median NUSE 0.981–1.039** (flag threshold 1.05). **Nothing flagged. All 34 arrays retained.**
- Two arrays remained mildly flagged on the within-group metric (**GSM2934803** BMD, 0.9676;
  **GSM2934813** control, 0.9784) and were kept, with leave-one-out testing whether either
  carried the result (§5.4 — neither did). [`02_rma_qc.R`]

## 4.3 Probe-to-gene mapping
- `hgu133plus2.db`. Funnel: **54,675** probesets → **44,647** with a gene symbol → **−1,524**
  dropped for mapping to more than one gene → **43,123** kept → **20,823 genes** after collapsing.
- **11,043 genes had multiple probesets. Collapse rule: highest mean expression** (project
  default). Rationale: among probesets for one gene the brightest generally hybridises properly.
- **No low-expression filter was applied** (project default). **Disclose that this was never
  sensitivity-tested.** [`03_annotate_collapse.R`]

## 4.4 Differential expression
- `limma`, design `~0 + group`, contrasts **BMD−Control, DMD−Control, BMD−DMD**;
  `lmFit` → `contrasts.fit` → `eBayes`.
- **Residual df 31, prior df 4.23** — state these; the prior df quantifies how much the empirical
  Bayes moderation is doing (roughly 4 extra samples' worth of information per gene).
  [`04_limma.R`]

## 4.5 microRNA families and target sets
- **TargetScanHuman 8.0**, human (species 9606):
  - `miR_Family_Info.txt` (dated 2018-04-06) for seeds — the `Seed+m8` column, positions 2–8.
    **2,064 human families.**
  - `Predicted_Targets_Info.default_predictions.txt` — **conserved-site default predictions**;
    115,308 rows, 219 families, 13,034 genes.
  - `UTR_Sequences.txt` for 3'UTR lengths; alignment gaps stripped, longest transcript per gene;
    **19,429 genes**.
- **Muscle filter: miRNA Tissue Atlas 2025** (Center for Bioinformatics, Saarland), rows with
  `Organ_system = limb_muscle`, **mean RPMM ≥ 10** → **353 of 2,656** microRNAs pass.
  **The threshold is ours, not from any publication, and was never sensitivity-tested — disclose
  it in Methods and again in Limitations.**
  - *Face-validity check worth one sentence:* the atlas's top limb-muscle microRNAs are
    miR-486-5p, miR-10b-5p, miR-133a-3p, miR-1-3p, miR-22-3p, miR-378a-3p, miR-206 — the
    textbook myomiRs, so the table behaves as muscle should.
- **Funnel: 2,064 families → 253 muscle-expressed → 128 with ≥ 15 array-present targets.**
  Sizes **15 to 1,446, median 371**. No upper limit (project default).
- **miR-29-3p:** seed **`AGCACCA`**; members hsa-miR-29a/b/c-3p; TargetScan conservation class 2;
  **1,158 array-present targets**; miR-29a-3p **1,895 RPMM** in limb muscle. [`05_mirna_inputs.R`]

## 4.6 Enrichment testing
- `fgsea` on genes ranked by the **limma moderated t**, no cutoff.
- **`nPermSimple = 10000`, `eps = 0`.** State why: at the default 1,000, **13 of 128 families
  returned NA NES**.
- Report **both** an fgsea adjusted p and a patient-label permutation p for every headline
  family. This is the paper's central methodological commitment — introduce it here, in Methods,
  not as a surprise in Results. [`06_smoke_test.R`, `07_fgsea_all.R`]

## 4.7 Patient-label permutation — define both nulls explicitly
This subsection carries the paper. Define the two nulls side by side:
- **Gene permutation (fgsea default):** patients and the ranked list are held fixed; random gene
  sets of matched size are drawn. Answers *"is this set special compared with other sets?"*
  Assumes genes are independent, which co-regulated collagens emphatically are not.
- **Patient permutation (ours):** the gene set is held fixed; the 34 disease labels are dealt out
  at random **preserving 6/11/17**, then limma, ranking and fgsea are rerun in full. Answers
  *"could disease status have produced this, or would random patient groups do as well?"*
- **500 shuffles** (plan said 1,000; measured throughput made that 5.7 h for this step alone, and
  the project rules permit 500 — log it as a deviation). Empirical p floor **1/501 = 0.002**.
- **p = (k + 1) / (N + 1)**. Family-level tests are **two-sided on |NES|**; the r and slope tests
  are **one-sided**. State both — they differ, and a referee checking reproducibility will notice.
- Step 9 shuffles **BMD/DMD among the 28 patients only, preserving 11/17**, with controls excluded
  from the test entirely. [`08_shuffle_baseline.R`, `09_bmd_vs_dmd.R`, `09b_classify.R`]

## 4.8 Bias controls
- **ECM removal:** all **275 core-matrisome genes** deleted from **both** the ranked list and the
  family sets (deleting from both keeps set sizes and scores comparable); families falling below
  15 genes excluded. [`10a_ecm_removed.R`]
- **3'UTR-matched null:** genes binned into **deciles of 3'UTR length**; **1,000** random sets
  drawn per family copying that family's per-bin counts exactly; empirical p floor **1/1001**.
  [`11_utr_matched.R`]
- **Fibroblast covariate:** per-sample score = mean row z-score of **316 array-present fibroblast
  markers that are not themselves fibrosis-set genes** (using fibrosis genes would be circular);
  limma refit as **`~0 + group + fibro`**. [`12_fibroblast_cov.R`]
- **Leave-one-out:** each BMD sample (11) and each control (6) dropped in turn — **17 complete
  pipeline re-runs**. [`13_leave_one_out.R`]
- **Sample-size control:** DMD randomly subsampled to 11, **100 draws**. [`14_subsample_dmd.R`]
- **Detection limit:** ruler = the DMD miR-29 footprint (**+0.2686 log2**, targets minus
  non-targets); a UTR-matched random set of **1,158** genes was shifted by *f* × ruler **in the 11
  BMD samples only**; BMD−Control limma and fgsea rerun; *f* ∈ {0, 10, 25, 50, 75, 100}%,
  **100 repeats each = 600 simulations**. A second script re-scored the **already-saved** 600 runs
  against a bar calibrated to the 95th percentile of the unplanted runs (**NES > 1.812**) — no
  recomputation. [`10b_detection_limit.R`, `10b2_calibrated_catch.R`]

## 4.9 Gene sets used
- **Fibrosis:** MSigDB 2026.1.Hs via `msigdbr` 26.1.1 — **`NABA_CORE_MATRISOME`** (Naba et al.
  2012, PMID 22159717) **275 genes**; **`HALLMARK_TGF_BETA_SIGNALING`** **54 genes**; 2 overlap.
- **Fibroblast markers:** MSigDB C8 `RUBENSTEIN_SKELETAL_MUSCLE_FAP_CELLS` +
  `..._FBN1_FAP_CELLS` (Rubenstein et al. 2020), union **372**; **326** not in a fibrosis set;
  **316** present on our array.

## 4.10 Software, seeds, reproducibility
- **R 4.6.2** (2026-10-29), **Bioconductor 3.23**, limma **3.68.5**, fgsea **1.38.0**, msigdbr
  **26.1.1**, BiocParallel **1.46.0**, data.table **1.18.6.1**, ggplot2 **4.0.3**, plus oligo,
  pd.hg.u133.plus.2, hgu133plus2.db, AnnotationDbi. All compiled from source. Full `sessionInfo()`
  in `results/sessionInfo.txt`.
- **Seeds:** 1 (steps 1–4); 42 (6, 7, 10a, 12); 20261101 (8); 20261102 (9); 20261103 (11);
  20261104 (10b); 20261105 (10b2).
- **Disclose the hardware change.** Step 10b ran on Ubuntu 26.04 / AMD EPYC 9R14 / 16 cores /
  30 GB rather than the 4-core Mac used for everything else. **Every relevant package version was
  identical.**
- **Runtimes** (optional, for a supplement): step 8 **206.7 min** (7 cores, 500 shuffles);
  step 10b **76.2 min** (16 cores, 600 simulations); step 13 **7.3 min**; step 14 **31.1 min**.

---

# 5. Results — in this order

Five results. R1 builds the method, R2 is the hard per-gene finding, R3 is the collapse, R4 is
the floor, R5 is the technical note.

## 5.1 R1 — the method recovers the known DMD biology (directionally)
**Where:** first Results subsection. **Figure 1** (ranked walk), **Figure 2** (NES heatmap).

- miR-29-3p, DMD vs control: **NES +2.224**, **2nd of 128 by NES**, 3rd by padj.
  Direction NES > 0 → targets up → brake weaker → miR-29 **less active** in DMD = **the published
  direction.**
- Mean log2FC of the **1,158** targets **+0.304** vs **+0.036** for the other 19,665 genes.
  Footprint size (the "ruler") **+0.269 log2 ≈ 20.5%** on the linear scale.
- **Leading edge is a roll-call of collagens** — COL1A1, COL1A2, COL3A1, COL4A1/2/4, COL5A1/2/3,
  COL6A2/3/6, COL15A1, COL16A1, COL21A1, COL27A1, plus ELN, FBN1, LAMA2, LAMC1, NID1, SPARC.
  **miR-29's leading edge is 6.5% fibrosis-set genes, the highest of any family** (next 4.1%;
  most 1–3%). *Lean on this: it is a statement about content, not p-values, and no permutation
  problem can touch it.*
- Survives the bias controls: ECM removal **+2.224 → +2.151 (96.7% retained)**; fibroblast
  adjustment **+2.224 → +1.762 (79% retained**, vs ~55–60% for other top families, and afterwards
  leading the field by eleven orders of magnitude: 1.7e−16 vs 3.1e−05 for the next family);
  UTR-matched null **z = 6.56**, 0 of 1,000 matched sets reached the real value; leave-one-out
  **17/17**, range +2.120 to +2.186.
- **Three qualifications that must appear in this same subsection, not later:**
  1. **103 of 128 families reached padj < 0.05 in this comparison and 99 pointed the same way.**
     The positive control was therefore **not** evidence that miR-29 is special — essentially
     every muscle family looked de-repressed. What is specific is the direction, the rank, the
     leading-edge content and the matched-null z-scores.
  2. **The cause is a globally tilted ranked list:** roughly 9,000–12,000 of 20,823 genes shift
     at FDR < 0.05 in dystrophic muscle (BMD−ctrl 4,549 up / 4,659 down; DMD−ctrl 5,547 / 6,654;
     BMD−DMD 4,564 / 4,039). A **competitive** rank test finds almost any set enriched when half
     the list is moving.
  3. **NES has no zero baseline.** A 3'UTR-matched random set already scores **NES ≈ +1.6**
     (1.592 in DMD, 1.669 in BMD; independently reproduced at 1.679 by step 10b's unplanted runs).
     **About 72% of miR-29's apparent enrichment is 3'UTR length alone.** Never quote an NES
     against zero.

## 5.2 R2 — the one robust per-gene finding
**Where:** second Results subsection. **Table 1.**

| Gene | BMD−Ctrl | adj.P | DMD−Ctrl | adj.P | BMD−DMD | adj.P |
|---|---|---|---|---|---|---|
| COL1A1 | +2.517 | 9.4e−08 | **+4.058** | 5.6e−13 | −1.542 | 4.6e−06 |
| COL1A2 | +2.043 | 1.4e−07 | +2.940 | 7.5e−12 | −0.897 | 4.2e−04 |
| COL3A1 | +2.311 | 8.3e−09 | +3.260 | 4.6e−13 | −0.949 | 1.3e−04 |
| COL5A1 | +1.919 | 1.4e−06 | +2.828 | 7.5e−11 | −0.910 | 8.5e−04 |
| LOX | +2.607 | 2.4e−07 | +3.573 | 3.7e−11 | −0.966 | 3.1e−03 |

- All five up in **both** diseases; all five up **significantly less in BMD**. COL1A1 is
  **2^4.058 ≈ 16.7×** in DMD versus **2^2.517 ≈ 5.7×** in BMD — Becker scars, roughly a third
  as hard.
- **Consistent with Eisenberg 2007** in direction and order of magnitude: they reported COL3A1
  about 6× and COL1A2 about 5× up in DMD; we get 9.6× and 7.7× on a different cohort and platform.
- **This is per-gene limma. No permutation is involved, so none of §5.3's collapse applies.**
- **⚠ The caveat that must travel in the same paragraph.** Nothing was adjusted for age, DMD
  median is 3.5 y and BMD median 12.2 y, and muscle ECM/collagen expression is highest during
  early postnatal growth. **The developmental expectation points the same way as the disease
  finding.** Report the difference as real and significant, and state that age and disease
  severity cannot be separated in this cohort. *(Hamza: this needs one sourced citation on
  age-dependence of muscle collagen expression. If none can be found by Saturday, write the
  caveat without the citation — do not drop it.)*
- **Delete the sentence in the reference document (Part 3.3) saying no BMD-vs-DMD differences
  were found so the age confound never got the chance to mislead us.** This table is a BMD-vs-DMD
  difference. The two statements contradict each other.

## 5.3 R3 — under patient permutation, nothing survives
**Where:** third Results subsection. **The core of the paper. Figure 3.**

| Claim | fgsea (gene shuffle) | **Patient shuffle** |
|---|---|---|
| miR-29, DMD vs control | padj 1.5e−34 | **p = 0.052** |
| miR-29, BMD vs control | padj 1.4e−31 | **p = 0.084** |
| miR-29, BMD vs DMD | padj 4.7e−10 | **p = 0.239** (q = 0.951) |
| Families differing BMD vs DMD | 7 of 128 | **0 of 128** |
| Families altered vs controls in both groups | 101 of 128 | **1 of 128** (miR-196-5p) |
| r between BMD and DMD footprints | 0.9123 | **p = 0.064** |

- **The correlation.** Real **r = 0.9123**, slope **0.9769**. Shuffled-label null: r mean 0.510,
  median **0.553**, SD 0.345, **95th percentile 0.933**, max 0.984 → **p = 0.0639**. Slope null
  median 0.555, 95th percentile 1.094, max 1.311 → **p = 0.1138**. **The real r does not clear
  the 95th percentile of pure noise.** Explain why with the shared-control mechanism: BMD−Ctrl
  and DMD−Ctrl use the same 6 controls, so any fluke in those 6 enters both comparisons in the
  same direction and manufactures a straight line through the origin. r = 0.91 is simply what
  this design produces when nothing is happening.
- **The direct comparison** (step 9) excludes controls entirely, so no control fluke can create a
  BMD–DMD difference. **0 of 128 families differ** under patient permutation (all q > 0.95);
  miR-29 NES −1.590, fgsea padj 4.7e−10, **shuffle p = 0.239**, null mean 0.053 / SD 1.386 from
  496 usable shuffles. Note the direction (brake predicted *more* active in Becker, the
  milder-disease direction) **only as a non-significant observation, if at all.**
- **Because 0 of 128 differ, no family can be BMD-only or DMD-only.** The classification lands
  entirely on which yardstick defines "altered versus controls": by fgsea padj, **101 shared /
  27 neither**; by patient shuffle, **1 shared / 127 neither**. miR-29 is "shared" under the
  first and "neither" under the second. **This one table is the paper in six numbers.**
- **The controlled test, which is the cleanest demonstration.** At *f* = 0 nothing was planted —
  1,158 random UTR-matched genes, unmodified, ground truth known to be nothing. fgsea called them
  significant in **100 of 100 runs**, and the **least** significant had **padj = 2.1e−06**. The
  false-alarm rate was **1.00 where 0.05 was expected.** Unplanted NES mean 1.679, SD 0.087,
  range 1.455–1.867.
- **Honest reporting of the borderline values.** With 500 shuffles the Monte Carlo standard error
  on an empirical p near 0.05 is √(0.05 × 0.95 / 500) ≈ **±0.01**. So **p = 0.052 is not
  resolvable from 0.05** — write that we cannot resolve it, and do not quote these p-values to
  more than two significant figures. **Never write "trending", "marginal", or "approaching
  significance".**
- **Do not double-count.** "miR-29 ranks 3rd of 128 on the honest test" is **p = 0.052 restated as
  a rank**, not independent evidence. Use one or the other. And it is **true for DMD only** — in
  **BMD vs control miR-29 ranks 45th**. Never let "3rd of 128" sit next to a Becker sentence.
- **Drop the "6.4 expected by chance" benchmark**, or caveat it heavily. It is 128 × 0.05, which
  assumes the 128 families are independent tests; target lists overlap heavily and miR-29 alone is
  5.6% of the measured transcriptome, so the variance of that count far exceeds binomial. The null
  conclusion does not depend on it.
- **Do not list the 10 families with BMD shuffle p < 0.05.** None survives correction for 128
  tests. **And the 14q32 / DLK1-DIO3 observation stays out entirely — it is `[UNVERIFIED]` and the
  counts are at chance level.**

## 5.4 R4 — the detection floor, and what actually limits us
**Where:** fourth Results subsection. **Figure 4** (detection curve), **Figure 5** (subsampling).

- **BMD's footprint is not detectably weaker than DMD's.** Full DMD NES 2.224; **DMD subsampled
  to n = 11**: median **2.167** (range 2.074–2.258); **actual BMD 2.169**. **53.0% of the
  DMD-at-n=11 draws fall at or below BMD's value** — a coin flip. BMD sits almost exactly at the
  median of what Duchenne produces at Becker's sample size.
- **Explain the apparent paradox with R2 explicitly, or a reader will think you contradicted
  yourself.** NES is a **rank** measure and rank measures saturate: in both diseases miR-29's
  targets sit firmly at the up end, so both score ≈ 2.2, even though BMD's collagens rose about a
  third as far. **The magnitude difference is real and visible only in the per-gene limma results
  or the direct contrast — never in the vs-control NES values.**
- **Robustness.** 17 leave-one-out runs: BMD vs control +2.120 to +2.186 (17/17 significant);
  BMD vs DMD −1.674 to −1.500, **rank 1 of 128 in every run**; top-20 family overlap 18–20 of 20.
  **GSM2934804** (disputed severity) dropped → NES +2.156 vs +2.169. **Immaterial; item closed.**
  State plainly that this proves the result is a property of the group, not one patient, and that
  **it does not rescue anything from R3** — every one of those runs used a gene-permutation p-value.
- **The detection curve, calibrated.** Against a bar set at the 95th percentile of the unplanted
  runs (NES > 1.812): catch rate **0.05 / 0.73 / 1.00 / 1.00 / 1.00 / 1.00** at
  *f* = 0 / 10 / 25 / 50 / 75 / 100%. 80% detection at **≈ 14%** of the DMD footprint
  (+0.037 log2). **Note that the 5% at f = 0 is true by construction** — the bar was set from
  those same 100 runs — so it is not an independent false-alarm estimate.
- **Keep the 14% out of the abstract and out of the main-text headline.** The whole of step 10b
  uses gene-sampling nulls, while the paper's conclusion rests on patient permutation, so the
  number is in the wrong currency for the claim it would sit beside. **Report the flat 100% line
  of the uncalibrated rule as the headline result** — it is the cleanest demonstration of R3 —
  and put the calibrated curve and its floor in the supplement.
- **The real conclusion of this subsection, and it is a finding in its own right.** BMD's actual
  miR-29-specific excess is **+0.095 log2 ≈ 35% of the ruler** (BMD gap 0.1923 minus the
  UTR-matched 0.0973; 3'UTR length alone accounts for **42%** of the DMD ruler and **51%** of
  BMD's gap). That is comfortably above the 14% floor, and at *f* = 35% the simulation catches it
  **100% of the time** — **yet the honest test does not detect it (p = 0.084)**. The reason is
  that planting adds **the same constant shift to every BMD patient**, whereas our Becker patients
  span **2.5 to 60.8 years** and mild to severe. **The binding limit at n = 11 is
  patient-to-patient variability, which a planting simulation of this kind cannot model.** The
  calibrated curve is therefore a **best-case floor**: the true detection limit is at least this
  large and probably much larger.

## 5.5 R5 — a technical note on duplicate accessions
**Where:** a short Results subsection or a Methods note. Two or three sentences.

- While searching for additional Becker samples we found that **GSE38417 and GSE109178 share raw
  arrays**. **GSM941830** (GSE38417 control) and **GSM2934813** (GSE109178 control) are
  **bit-identical**: **1,354,896 probe cells compared, r = 1.00000000, 0 cells differing.** Both
  are HG-U133 Plus 2.0, both carry the same lab-internal array identifiers
  (`MIDO-CTRL_Pos-1UP` …), and only the accession and the CEL encoding (text vs binary) differ.
- **Scope the claim to the evidence: one verified identical pair plus matching internal
  identifiers.** Do **not** write that all arrays across both accessions were compared — they
  were not. Either weaken the sentence accordingly or run `scripts/check_duplicate_arrays.R`
  across all pairs before submission and report what it finds.
- **Frame it as a technical caution, not an allegation.** Both deposits are from the same group
  (Hoffman, CNMC); same-group re-deposit is routine. The point is that **metadata will not reveal
  it and merging would double-count samples**, halving standard errors with no new information —
  a real temptation when n = 11 is the field's binding constraint. Ship
  `scripts/check_duplicate_arrays.R` as the practical output.
- **Do not write that merging "would have crossed 0.05".** That is speculation and is recorded as
  such in `docs/fact_check.md`.

---

# 6. Discussion — the argument in order

## 6.1 What the method achieved
It recovers the known Duchenne biology directionally and specifically: right direction, 2nd of
128 by NES, the most fibrosis-dense leading edge of any family, and survival of matrisome
removal, fibroblast adjustment, a UTR-matched null and leave-one-out. **That is a real validation
of the footprint approach as a descriptive tool.**

## 6.2 What the biology says
We cannot confirm that the miR-29–fibrosis axis operates detectably in human BMD, and we cannot
distinguish BMD's predicted microRNA activity from DMD's. The field's default assumption is that
Becker is Duchenne-dialled-down; **we tested that directly, with the correct null, and could not
confirm it at n = 11** — while being able to state quantitatively what we could and could not
have seen. **A null result with a measured detection floor is a contribution; a null result
without one is a failed experiment.** Note that the per-gene collagen magnitudes *do* differ
(R2) — so "Becker scars less" is supported, while "Becker's microRNA footprint differs" is not.
Keep those two statements rigorously apart.

## 6.3 What the methods result is — and is not
**The framing here decides whether this section is a strength or an embarrassment.**

- **✅ Claim:** on microRNA target sets this large and this correlated, in a tissue whose
  transcriptome is globally shifted, gene-permutation adjusted p-values were **irreconcilable
  with patient-label permutation** — padj 1.5e−34 against p = 0.052 for the same gene set and
  contrast — with a **measured false-alarm rate of 1.00 where 0.05 was expected**.
- **🚫 Do not claim this is a new phenomenon, and do not claim fgsea is broken.** The
  gene-sampling versus subject-sampling distinction is established: **Goeman & Bühlmann 2007**
  (the canonical treatment of competitive vs self-contained tests), **Wu & Smyth 2012** (`camera`,
  which exists specifically to correct enrichment tests for inter-gene correlation and ships in
  limma), **Subramanian et al. 2005** (GSEA's own recommendation of phenotype permutation), and
  **Efron & Tibshirani 2007** (restandardisation). **Cite them.** Our contribution is the
  **magnitude in a specific, severe regime** — not the discovery.
- Frame the failure as a **mismatch between the null and the hypothesis**: fgsea correctly answers
  *"is this set special compared with other sets?"*, and we were using that answer for *"does
  disease status drive this?"* That framing is accurate and much harder to attack than "the tool
  is wrong".
- **Anticipate the obvious question — "why not just use `camera` or `fry`?"** Answer it honestly:
  patient-label permutation was the pre-specified test, computation was frozen on 6 November, and
  a correlation-adjusted parametric test is the natural next step. Put it in Future Work rather
  than letting a referee put it in their report.
- State the two features that make our case especially severe, so the finding is not
  overgeneralised: **set size** (miR-29's 1,158 targets are 5.6% of the measured transcriptome,
  and target lists overlap heavily) and the **globally tilted list** (~half the transcriptome
  moves; 103 of 128 families "significant", 99 in the same direction).

## 6.4 Reconcile our own evidence — before a referee does it for us
**This paragraph is not optional.** Our bias controls — ECM removal, the UTR-matched null, the
fibroblast covariate, leave-one-out — all assess significance by **sampling genes**, which is the
very null we criticise. Say so ourselves, and draw the distinction:

> Those controls answer a different and legitimate question — *is miR-29's footprint an artifact
> of scar content, fibroblast infiltration or 3'UTR length, and is it special among comparable
> gene sets?* — and the answer is no, it is not an artifact. **None of them can establish that
> disease status drives it.** Only patient permutation addresses that, and it does not.

Without this paragraph the paper reads as arguing gene permutation is catastrophically wrong and
then relying on it for every positive claim.

## 6.5 What would settle it
A direct miR-29 measurement in BMD muscle — one experiment would settle what fourteen analysis
steps could not. Failing that, a cohort large enough to absorb Becker's clinical heterogeneity.
**No public dataset fixes n = 11:** GSE38417 duplicates our arrays (R5); GSE13608 buys 6 → 12
controls for ~18% precision (0.508 → 0.417) at the cost of unknown biopsy sites, ages 26–85 and
one pooled sample; GSE291383 has 6 BMD but **no controls** and is RNA-seq.

---

# 7. Limitations — these seven must appear

Written out, not hedged into invisibility. A referee will find them anyway.

1. **Sample size.** 11 BMD, 6 controls. The binding constraint on everything, and no public
   dataset fixes it.
2. **Six controls shared between two comparisons.** Quantified rather than hand-waved: shuffled
   labels produce r ≥ 0.933 on 5% of runs, against our real r = 0.9123.
3. **Age confound.** BMD median 12.2 y (2.5–60.8); DMD median 3.5 y (0.9–8.0); **control ages are
   not reported in GEO at all** `[UNVERIFIED]`. Age and disease severity are entangled, nothing
   was adjusted, and **the confound runs in the same direction as the collagen finding**. It is
   also plausibly part of *why* the result is null: Becker's 58-year age span is itself a source
   of the between-patient variability that R4 identifies as the binding limit.
4. **No microRNA measurements.** Everything is inferred from predicted targets; the array cannot
   see microRNAs. Targets are **TargetScan conserved-site predictions**, which contain false
   positives and miss real sites.
5. **3'UTR-length bias.** About **72%** of miR-29's raw enrichment is attributable to UTR length
   alone. Any study quoting NES against a baseline of zero on target sets is subject to this.
6. **Tissue composition.** Dystrophic biopsies contain fewer muscle cells and more fibroblasts
   (fibroblast score: control **−0.791**, BMD **+0.040**, DMD **+0.253**, ANOVA
   **p = 1.2e−08** — Becker sits neatly in between, which is itself a clean independent result
   worth one sentence). Step 12 adjusts for it, **but the covariate also absorbs real disease
   signal**, so neither the adjusted nor the unadjusted number is a clean effect estimate.
7. **Analytical choices that were never sensitivity-tested:** the **mean RPMM ≥ 10** muscle
   threshold is ours and appears in no publication; **no low-expression filter** was applied, so
   all 20,823 genes including many weakly expressed probesets entered the ranking; **500** label
   shuffles rather than 1,000; and the **detection simulation is a best case** that plants a
   constant shift on uncorrelated genes and cannot model between-patient variability. Also:
   **no cross-species validation** — the intended pig comparison was dropped because the source
   table could not be obtained (mention as an intended analysis, never as a pending result).

---

# 8. Figures and tables

## Main
| # | File | Shows | Caption must include |
|---|---|---|---|
| **1** | `06_mir29_enrichment_DMDvsCtrl.png` | The ranked walk | NES +2.224; direction means miR-29 predicted **less** active |
| **2** | `07_nes_heatmap.png` | NES across three comparisons, top families | The vs-control tilt and its collapse in BMD−DMD |
| **3** | `08_shuffle_r.png` | **The pivotal figure** | Real r = 0.9123; 500 shuffles; **95th pct 0.933 marked**; p = 0.064 |
| **4** | `10b2_detection_curve_calibrated.png` | Detection curve, both rules | **The flat line at 100% for the spec rule IS the finding** |
| **5** | `14_subsample_dmd.png` | BMD at the median of DMD-at-n=11 | 53.0% of draws at or below 2.169 |
| **T1** | `04_sanity_fibrosis_genes.tsv` | The collagen table (§5.2) | Age caveat in the footnote |

**Figure 3 is the one a referee will stare at longest.** It must show the histogram of 500
shuffled r values, the real r marked, **and the 95th percentile at 0.933 marked so the reader can
see the real value sits below it.**

**Caption to adapt:** *"Correlation between BMD and DMD microRNA footprints (r = 0.9123, red)
against 500 permutations of patient labels preserving group sizes (grey). Shuffled labels produce
r ≥ 0.933 on 5% of permutations. Empirical p = 0.064. The apparent agreement between the two
footprints is not distinguishable from what the shared-control design produces by chance."*

## Supplementary
`02_qc_boxplot.png`, `02_qc_nuse.png`, `02_qc_within_group_cor.png` (including the group-size
artefact); `09_mir29_bmd_vs_dmd_shuffle.png`; `12_fibroblast_score.png`;
`13_mir29_leave_one_out.png`; the calibrated detection curve and its floor; the classification
table (101 vs 1).

**🚫 No volcano or scatter plots.** If a BMD-vs-DMD NES scatter is unavoidable it must be printed
**next to** Figure 3 — on its own it looks like overwhelming evidence and is misleading for
exactly the reason in §5.3.

Per-figure axis notes and worked example points are in `docs/figure_notes.md`.

---

# 9. Sentences that are forbidden

| Never write | Why |
|---|---|
| miR-29 is significantly altered in BMD or DMD | Honest p = 0.084 / 0.052. Neither clears 0.05. |
| BMD's footprint is a weaker copy of DMD's | 53% of DMD-at-n=11 draws fall at or below BMD; r p = 0.064 |
| BMD and DMD differ in microRNA family X | **0 of 128**, all q > 0.95 |
| The direct comparison shows miR-29 significantly less elevated in BMD | fgsea padj 4.7e−10, **honest p = 0.239** |
| We detect footprints down to 14% of DMD's | Best case only; the real limit is patient variability |
| r = 0.91 shows the footprints agree | Shuffled labels reach 0.933 on 5% of runs |
| Any NES quoted against a baseline of zero | The matched baseline is **≈ +1.6** |
| 80% detection reached at f = 0% | A script artefact of the broken catch rule |
| "miR-29 ranks 3rd of 128 on the honest test" near a BMD sentence | True for **DMD only**; BMD is **45th** |
| fgsea is broken / this invalidates GSEA | Null–hypothesis mismatch, not a bug. See §6.3 |
| "overconfident by 32 orders of magnitude" / "by a factor of 2.3e8" | The permutation p floor is 1/501 ≈ 0.002; you have no resolution to quantify the gap, and a ratio of p-values is not a meaningful quantity. Write "padj 1.5e−34 versus p = 0.052 — not significant" |
| p = 0.052 is "trending" / "marginal" / "approaching significance" | ±0.01 Monte Carlo error at 500 shuffles; unresolvable |
| Anything citing the old preliminary GEO2R result | Superseded by the limma tables |
| The 14q32 / DLK1-DIO3 cluster observation | `[UNVERIFIED]`, counts at chance level |
| Merging the duplicate "would have crossed 0.05" | Speculation, recorded as such |

---

# 10. Clear before submission

**Do not cite any of these until checked against the source.** From `docs/fact_check.md`.

| Item | Status | Owner |
|---|---|---|
| Control ages | **Absent from GEO entirely.** Write `[UNVERIFIED]`; it cannot be resolved | — |
| **Spot-check 2–3 GSM titles by hand on the GEO page** | Labels were parsed from GEO web endpoints and validated only against the expected 6/11/17/15 counts — a swap *inside* a group would be invisible. Ten minutes; it is the foundation under everything | Aziz |
| miRNA Tissue Atlas 2025 exact citation | Lineage recorded (gkw116 → gkab808 → gkae1036); which is the 2025 release is unchecked | Hamza |
| Rubenstein et al. 2020 PMID 32042044 | Taken from MSigDB metadata, not checked against the paper | Hamza |
| Age-dependence of muscle collagen expression | Needed for the §5.2 caveat; if unsourced by Saturday, write the caveat without it | Hamza |
| **Goeman & Bühlmann 2007; Wu & Smyth 2012; Subramanian 2005; Efron & Tibshirani 2007** | **Required for §6.3. Highest-priority citations in the paper** | Hamza |
| Eisenberg "5 BMD", "GSE3307 = same dataset" | Not found. **The checked figure is 6 BMD** | — |
| Xie 2024 gene counts; the "13 genes" claim; Dadgar 2014 as the source of GSE109178; GSE157668 contents; the pig "18 upregulated" | Unverified. **Do not use** | — |
| GSM2934804 severity | Disputed; immaterial either way (NES +2.156 vs +2.169) | closed |

**Cleared and safe to cite:** GSE109178 = 49 samples on GPL570, all vastus lateralis; the
6/11/17/15 split with LGMD as 8 dysferlin + 7 FKRP; BMD ages 2.5–60.8; DMD ages 0.9–8.0; miR-29
seed `AGCACCA` with members 29a/b/c-3p; and Eisenberg 2007's 6 BMD samples, miR-146b/miR-221 as
the only BMD passes, linkage significant only in DMD and Miyoshi, and DMD miR-29c about 6-fold
down with COL3A1 about 6× and COL1A2 about 5× up.

---

# 11. Division of labour for the remaining days

| | Aziz | Hamza |
|---|---|---|
| **Fri 7 Nov** | Methods (§4) and Results (§5) from this file | Intro/Related Work; the §6.3 citations |
| **Sat 8 Nov** | Discussion (§6) and Limitations (§7); figure captions | Limitations review; Data section; the age-collagen citation |
| **Sun 8 Nov** | Both: check every citation against its source, proofread, submit | |

**Computation is frozen.** The only permitted runs are the two verification tasks in §10 (the GSM
label spot-check, and optionally `check_duplicate_arrays.R` across all pairs if you want the
stronger version of R5). Nothing else is to be rerun.

---

*Compiled 7 November 2026 from `docs/FULL_RESEARCH_DOCUMENT.md`, `docs/results_log.md`,
`docs/methods_log.md`, `docs/fact_check.md` and `docs/figure_notes.md`, with the four framing
corrections listed at the top. Every number traces to a script in this repository that ran to
completion.*
