# CLAUDE.md: BMD microRNA–fibrosis project (full research)

> ## ⚠ STATUS BLOCK — updated 3 Nov 2026. Read before anything else.
>
> **Steps 0–14 are COMPLETE** (10b ran 3 Nov on a 16-core Linux box, see below). **Step 15 is
> BLOCKED** (Reschke Table 1 unavailable). Computation freezes **Friday 6 Nov**.
>
> **Step 10b: the spec's catch rule failed (false alarm 1.00 at f = 0).** Any UTR-matched random
> set scores NES ≈ 1.68 and padj ≤ 2.1e−06. Calibrated against the unplanted runs
> (`10b2_calibrated_catch.R`), 80% detection comes at ≈ 14% of the DMD footprint, but that is a
> best case. BMD's real miR-29 footprint (≈ 35% on that scale) would pass it, yet fails the
> patient shuffle (p = 0.084). **The binding limit is patient-to-patient variability, which 10b
> can't simulate.** Which catch rule to headline is Aziz's call. ❌ Don't write "we detect
> footprints down to 14%".
>
> ### The headline changed. The project is a NULL RESULT.
>
> The smoke test passed on fgsea's own terms (miR-29 DMD-vs-ctrl NES +2.22, padj 5.4e-35),
> but **fgsea permutes genes, not people**. Under the patient-label shuffle that §5 and §10
> demanded, everything dissolves:
>
> | Claim | fgsea said | Honest (patient-shuffle) |
> |---|---|---|
> | miR-29, DMD vs ctrl | padj 1.5e-34 | **p = 0.052** |
> | miR-29, BMD vs ctrl | padj 1.4e-31 | **p = 0.084** |
> | miR-29, BMD vs DMD | padj 4.7e-10 | **p = 0.239** |
> | families differing BMD vs DMD | 7 of 128 | **0 of 128** |
> | r between BMD and DMD footprints | 0.9123 | **p = 0.064** (null 95th pct = 0.933) |
>
> **fgsea is overconfident by up to 34 orders of magnitude on this data.** §10's warning was
> correct and the magnitude is a publishable finding in itself.
>
> ### Do NOT write these sentences
> - ❌ miR-29 is significantly altered in BMD or DMD
> - ❌ BMD's footprint is a weaker copy of DMD's (step 14: 53% of DMD-at-n=11 draws fall at or
>   below BMD; step 8: r p = 0.064)
> - ❌ BMD and DMD differ in any microRNA family
>
> ### What survives
> - The method recovers the right biology directionally: miR-29 targets up, collagen genes
>   leading the edge, surviving ECM removal (96.7% of NES), fibroblast adjustment (79%), a
>   3′UTR-matched null (6.5 SD), and 17/17 leave-one-out runs.
> - miR-29 ranks **3rd of 128 even on the honest test** — the signature of *underpowered*,
>   not *absent*.
> - **The paper is now the §6-step-6 fallback:** *"Can target footprints recover known
>   microRNA biology in dystrophic muscle? Detection limits at n = 11."* Step 10b carries it.
>
> ### Two findings worth publishing on their own
> 1. **fgsea's gene-permutation p-values are catastrophically overconfident** on correlated
>    target sets (above).
> 2. **GSE38417 is a duplicate re-deposit of GSE109178** — identical raw intensities across
>    all 1,354,896 probe cells (r = 1.000000, 0 differing). Merging it would have faked
>    significance. **Run `scripts/check_duplicate_arrays.R` before merging ANY accession.**
>
> ### Closed decisions (do not reopen)
> - **QC STOP (§6 step 2): keep all 34.** RLE/NUSE clean on every array; the spec's
>   correlation metric was confounded by group size (28 diseased vs 6 healthy guarantees
>   controls score lowest). Step 13 confirmed no sample matters.
> - **GSM2934804 severity robustness item: CLOSED.** Dropping it moves miR-29 from +2.169 to
>   +2.156. Immaterial.
> - **Shuffle count: 500**, not 1,000 (§10's permitted fallback). Measured throughput was
>   20.4 s/shuffle; 1,000 would have been 5.7 h for step 8 alone.
> - **No dataset expansion.** GSE38417 duplicate; GSE13608 adds muscle-site and age confounds
>   for ~18% precision; GSE291383 (6 BMD, no controls) is RNA-seq. **11 BMD is the binding
>   constraint and no public data fixes it.**
>
> ### Engineering rules learned the hard way
> - `fgsea` inside `bplapply` **must** take `nproc = 1`, else it grabs every core per worker
>   (52% of CPU lost to context-switching; 98 min job projected to 6.8 h).
> - `setorder(dt, -abs(x))` is **invalid** in data.table — crashed a step after 2 h of work.
> - Scripts must `saveRDS` expensive objects **before** any optional block.
> - Never run two heavy R jobs concurrently on this machine (8 GB, 4 real cores).
> - Project runtimes from a **parallel** batch, never a serial timing loop.
> - R lives at `/opt/homebrew/bin` — not on the default PATH.
> - NCBI FTP returns 403 here; use GEO's web endpoints (already wired into script 01).
>
> **Plain-language write-up for Aziz: `docs/NOTION_session_01.md`.**

---


Read this whole file before doing anything in this repo. It is the single source of truth for the analysis. Hamza’s narrower brief lives in `docs/CLAUDE-hamza.md`.

**First draft due: Sunday 8 November 2026.** We are behind schedule. Favour momentum.

---

## 0. How to work with Aziz

- **Aziz Al-Shujaa** leads the analysis. He is a CS student with limited biology and no strong stats background, and he must be able to explain and defend every step to the professor.
- **You write and run the code.**
    - **Before each step:** explain in 2–4 lines what it does and why. Include a tiny worked example: a few numbers, a small table, or an ASCII sketch.
    - **After each step:** show the key output and how to read it.
    - Never explain a stats point in prose alone.
- **Scatter plots, correlation (r) and slopes don’t land with him on their own.** Whenever you show one, add a worked example with real numbers. For the BMD–DMD comparison, always pair it with the shuffled-baseline histogram that marks the real value.
- **Keep it direct and concise.**
    - Where this file gives a default, use it, log it (§9) and carry on.
    - Stop only at the checkpoints marked **STOP**.
    - If you got something wrong, say so plainly and fix it.
- **Biology for CS readers:** explain at a basic level whenever it comes up.

---

## 1. The project

1. Becker muscular dystrophy (BMD) is the milder form of Duchenne (DMD). In both, muscle is slowly replaced by collagen scar (**fibrosis**).
2. MicroRNAs are ~22-letter RNAs that act as **brakes** on many genes at once. Their 7-letter **seed** (letters 2–8) pairs with sites in a message’s tail (the **3′UTR**), and that message gets destroyed or blocked. In DMD, **miR-29** (seed `AGCACCA`) is reduced, and it is a brake on collagen genes.
3. Nobody has tested whether the same microRNA–fibrosis link runs in human BMD. This matters because exon-skipping therapy turns DMD into a BMD-like state.
4. Our data contain **no microRNA measurements**. We infer microRNA activity from **footprints**: when a brake weakens, its target genes drift up together.
5. Core question: is BMD’s footprint a weaker copy of DMD’s, a different pattern, or nothing we can detect?

**Direction rule (use everywhere):**

| What fgsea shows for a family’s targets | Prediction |
| --- | --- |
| NES > 0: targets bunch at the UP end | Family predicted **less active** (brake weaker) |
| NES < 0: targets bunch at the DOWN end | Family predicted **more active** |

A measured microRNA that goes **up** should therefore show **NES < 0**.

---

## 2. Data

- **GSE109178.** Affymetrix HG-U133 Plus 2.0 (GPL570), vastus lateralis (thigh) muscle biopsies, 49 samples.
- **Keep 34 samples: 6 control, 11 BMD, 17 DMD.** Drop the 15 limb-girdle samples (7 LGMD2I and 8 LGMD2B).
- **Don’t assume label column names.** Inspect `pData()` and match on the actual titles and characteristics. Confirm 6 / 11 / 17 before going on.
- **Known context:**
    - BMD ages run from 2.5 to 60.8 years, and most cases are mild.
    - DMD patients are children.
    - **Control ages are unknown so far.** Record them from GEO. If the controls are adults, that is an age confound: state it, don’t fix it.
- GSM2934804 is a BMD sample with a disputed severity. It is a robustness item only (§6, step 13).

---

## 3. Repo layout

```
CLAUDE.md
scripts/        00_setup.R … 15_pig_overlap.R  (one script per step in §6)
inputs/         Hamza's files (exact formats in §4)
inputs/auto/    files you built because Hamza's weren't there yet
data/raw/       CEL files (gitignored)
data/processed/ RMA matrix, gene-level matrix (.rds)
results/tables/ results/figures/ results/rds/
docs/           methods_log.md, results_log.md, handoffs.md, figure_notes.md, fact_check.md, CLAUDE-hamza.md
```

Every script must run on its own: `Rscript scripts/NN_*.R`. Each script reads from disk and writes to disk, and sets `set.seed()` wherever randomness enters.

---

## 4. Inputs and Hamza’s part: don’t wait

Hamza owns the input files below. **If a step needs a file that isn’t in `inputs/` yet, build it yourself to the same format, write it to `inputs/auto/`, and carry on.**

- **Loading rule:** read `inputs/<file>` if it exists, else `inputs/auto/<file>`. Print which one was used.
- **Logging:** record every auto-built file in `docs/handoffs.md` with its source, version, date and method, so Hamza can check it.
- **When Hamza’s version arrives:**
    1. Diff it against the auto version.
    2. Report the differences to Aziz.
    3. Switch to Hamza’s version.
    4. Rerun the affected steps, unless results are already frozen (§7).

Gene symbols are HGNC, uppercase. Record the source **and version** of every list.

| File | Columns | Needed by |
| --- | --- | --- |
| `sample_sheet.csv` | `gsm, title, group` (control/BMD/DMD/LGMD), `age_years, sex, severity` (mild/moderate/severe/NA), `include` (TRUE/FALSE), `notes` | Step 1 |
| `mirna_families.tsv` | `family_id, seed_7mer` (positions 2–8), `members` (`;`-separated), `conservation, expressed_in_muscle` (TRUE/FALSE), `expression_source` | Step 5 |
| `mirna_targets.tsv` | `family_id, gene_symbol, source` (e.g. TargetScan 8.0), `site_type` | Step 5 |
| `utr_lengths.tsv` | `gene_symbol, utr3_length_nt, transcript_used` | Step 11 |
| `fibrosis_sets.tsv` | `set_name` (matrisome_core / hallmark_tgfb), `gene_symbol, source_version` | Step 10 |
| `fibroblast_markers.tsv` | `gene_symbol, source, in_fibrosis_set` (TRUE/FALSE) | Step 12 |
| `pig_mirnas.tsv` | `mirna_as_in_paper, human_name, direction` (up/down, BMD-like vs WT), `log2fc, padj, table_ref` | Step 15 |

**How to build them if needed:**

- **Seed families and targets:** use TargetScanHuman 8.0 (data download page). The family info file gives seeds and conservation. The default predictions file gives targets per family; filter it to human, species 9606. Confirm the exact file names on the page and log them.
    - **Unblocker:** the smoke test (step 6) needs only the miR-29 family (seed `AGCACCA`). Build that first if everything else is late.
- **Muscle filter:** propose a published source with a citation, such as a tissue microRNA atlas, to Aziz. Use the first credible one if he doesn’t answer within the session, and log it.
- **3′UTR lengths:** from TargetScan’s human UTR sequences (strip gap characters), one representative transcript per gene.
- **Fibrosis sets:** the human core matrisome list (Naba et al., MatrisomeDB) and `HALLMARK_TGF_BETA_SIGNALING` via `msigdbr`.
- **Fibroblast markers:** a published skeletal-muscle source. Set `in_fibrosis_set` to TRUE for any marker that is also in the fibrosis sets.
- **Pig list:** only from Reschke et al. 2026 (J Cachexia Sarcopenia Muscle), Table 1, with the directions. If you can’t read the table, don’t build this file. Mark step 15 as blocked.

---

## 5. Decisions already made (don’t reopen without Aziz)

- **Main method:** `fgsea`, ranked by the limma moderated **t**, on **muscle-expressed seed families**. There is no cutoff. The hypergeometric overlap test is **cut**.
- **Positive control first (the smoke-alarm test):** the miR-29 family must come out positive in DMD vs control. If it doesn’t, the paper becomes a method / detection-limit paper. The go/no-go is decided **Monday 2 Nov evening**.
- **The shared-control trap.** BMD-vs-control and DMD-vs-control use the **same 6 controls**, so any fluke in those 6 enters both comparisons:
    
    
    | Gene (no true difference) | Ctrl avg | BMD avg | DMD avg | BMD−Ctrl | DMD−Ctrl |
    | --- | --- | --- | --- | --- | --- |
    | 1 | 9.70 | 10.05 | 9.98 | +0.35 | +0.28 |
    | 2 | 10.30 | 9.96 | 10.02 | −0.34 | −0.28 |
    - In simulation with pure noise and these group sizes, the BMD and DMD footprint scores still correlate at **r ≈ 0.7 with a slope below 1**. That is the “weaker copy” picture, produced by noise.
    - **Never judge r or the slope against 0. Judge them against the shuffled-label baseline** (step 8).
    - “Is BMD different from DMD?” is answered by the **direct BMD vs DMD comparison** (step 9), where the controls cancel out.
- **Detection limit:** use the **DMD miR-29 footprint** as the ruler, not the collagen change (step 10b).
- **Cut:** GSE3307 replication, serum microRNAs (SQ4), the hypergeometric test, and GSE157668 unless it has a ready-made count table. Log each cut with its reason in `docs/methods_log.md`.

---

## 6. Pipeline

Each step lists what “done” looks like.

**0. Setup** (`00_setup.R`)

- Install via BiocManager: `GEOquery`, `oligo`, `pd.hg.u133.plus.2`, `hgu133plus2.db`, `AnnotationDbi`, `limma`, `fgsea`, `msigdbr`, `BiocParallel`, `data.table`, `ggplot2`.
- Write `sessionInfo()` to `results/sessionInfo.txt`.

**1. Download and labels** (`01_download_labels.R`)

- Start `getGEOSuppFiles("GSE109178")` first; it’s the slow part. Untar into `data/raw/`.
- Take the labels from `getGEO("GSE109178")` → `pData()`.
- Build or verify the sample sheet. Confirm 6 / 11 / 17 / 15 and record control ages.
- **Done:** the sample sheet exists, and the counts are printed and match.

**2. RMA and QC** (`02_rma_qc.R`)

- Read only the 34 included CEL files with `oligo::read.celfiles`, then run `rma()`. Restricting to 34 keeps the LGMD samples from shaping the normalization.
- QC without scatter plots:
    - a boxplot of the 34 samples (the boxes should line up);
    - each sample’s median correlation with all the others. Flag any sample clearly below the rest.
- **Done:** a 54,675 × 34 matrix in `data/processed/rma.rds` and a QC figure. **STOP** if any sample is flagged: Aziz decides keep or drop.

**3. Probes to genes** (`03_annotate_collapse.R`)

- Map probes to symbols with `hgu133plus2.db`. Drop unmapped probes and probes that map to several genes.
- Collapse to one probe per gene, keeping the probe with the **highest mean expression**.
- **Done:** a gene × 34 matrix, with the counts (probes kept, genes kept) logged.

**4. limma** (`04_limma.R`)

- Design `~0 + group` with three contrasts: **BMD−Ctrl, DMD−Ctrl, BMD−DMD**. Run `lmFit`, `contrasts.fit`, `eBayes`.
- Save the full tables (logFC, t, P, adj.P) to `results/tables/limma_<contrast>.tsv`.
- Sanity check: are COL1A1, COL3A1, COL5A1 and LOX up in DMD and in BMD? Show a small table of logFC and adj.P.
- Also count the genes at FDR < 0.05, up and down.
- **Done:** three tables plus the sanity table. This replaces the old preliminary GEO2R claim, which must not be cited.

**5. MicroRNA inputs** (`05_mirna_inputs.R`)

- Load the families and targets (§4). Keep muscle-expressed families only, and keep only targets present in the gene matrix.
- Size filter: at least **15** targets after intersecting with the gene matrix. Don’t cap large families.
- **Done:** `results/rds/family_sets.rds` (a named list) and a table of the families with their sizes.

**6. Smoke-alarm test** (`06_smoke_test.R`)

- Run `fgsea` on the DMD−Ctrl t-ranks over **all** retained families, with a fixed seed.
- Report the miR-29 family’s NES, p, padj and rank, and save its enrichment plot. This is the “ranked walk” Aziz already knows.
- Verdict:
    
    
    | Verdict | Condition |
    | --- | --- |
    | **GO** | NES > 0 and padj < 0.05 |
    | **BORDERLINE** | NES > 0 and nominal p < 0.05, but padj ≥ 0.05 |
    | **NO-GO** | Anything else |
- **STOP.** Show the result. The go/no-go is a joint call with Hamza on Monday evening.
- **If NO-GO:** the paper becomes “Can target footprints recover known microRNA biology in dystrophic muscle? Detection limits at n = 11.” Steps 10b, 11 and 12 then carry it, and steps 8–9 shrink.

**7. All families, three comparisons** (`07_fgsea_all.R`)

- Run fgsea for BMD−Ctrl, DMD−Ctrl and BMD−DMD.
- Output:
    - one table with families as rows and NES / padj for each comparison as columns;
    - the leading-edge genes of each top family, marked if they are in the fibrosis sets. This edge table is the microRNA → fibrosis-gene “network”.

**8. Shuffled-label baseline** (`08_shuffle_baseline.R`)

- Real values: r and slope between the BMD−Ctrl NES and the DMD−Ctrl NES across families.
- Null: shuffle all 34 labels **keeping 6/11/17**, rerun limma (both contrasts) and fgsea, and record r and slope. Default **1,000 shuffles**.
- Time 5 iterations first and extrapolate. If 1,000 won’t finish overnight, drop to 500 and log it. Use `BiocParallel`.
- Report the empirical p as (number of shuffles with r ≥ real r, plus 1) ÷ (N + 1), and the same for the slope.
- Free by-product: each family’s null NES distribution gives a **label-shuffle p** for the headline families (§10).
- **Figure:** a histogram of the shuffled r values with the real r marked, plus the same for the slope. Show it next to any scatter plot.

**9. Direct BMD vs DMD** (`09_bmd_vs_dmd.R`)

- Use the BMD−DMD fgsea from step 7, plus a label-shuffle test that shuffles BMD/DMD labels among the 28 patients (11/17) and reruns.
- Report which families differ beyond the shuffles. Shared / BMD-only / DMD-only families come from here, not from step 8.

**10. Two bias checks**

- **(a) ECM removed** (`10a_ecm_removed.R`): drop the matrisome-core genes from the ranked lists and rerun fgsea for DMD−Ctrl and BMD−Ctrl. Does miR-29, or each top family, survive? This asks whether the footprint is just scar genes.
- **(b) Detection limit** (`10b_detection_limit.R`):
    1. **Ruler:** the size of the DMD miR-29 footprint, measured as the mean DMD−Ctrl logFC of miR-29 targets minus that of non-targets.
    2. **Plant:** take a random gene set the size of the miR-29 family, matched on 3′UTR length. Add *f* × ruler to those genes in the 11 BMD samples only.
    3. **Rerun:** BMD−Ctrl limma and fgsea, with the planted set tested among all families.
    4. **Count a catch:** NES > 0 and padj < 0.05.
    5. **Grid:** *f* ∈ {0, 10, 25, 50, 75, 100}%, about 100 repeats each (time it first). At *f* = 0 the catch rate is the false-alarm rate.
    - **Figure:** catch rate against *f*. The detection limit is where the curve crosses 80%.

**11. 3′UTR-length-matched null** (`11_utr_matched.R`)

- Genes with long tails sit on many target lists, so a plain random baseline is biased.
- For each top family:
    1. Bin all genes into 3′UTR-length deciles.
    2. Draw 1,000 random sets with the family’s per-bin counts.
    3. Run fgsea on them.
    4. Empirical p = the share of matched sets with NES at least as extreme.

**12. Fibroblast covariate** (`12_fibroblast_cov.R`)

- Fibroblast score per sample: the mean row-z-score of the markers where `in_fibrosis_set` is FALSE.
- Refit limma as `~0 + group + fibro`, rerun fgsea and compare with step 7.
- Caveat to log: scar load differs by group, so the covariate also absorbs some real disease signal. Read this as “does the footprint survive”, not as the main result.

**13–15. If time** (the robustness checks need no new inputs)

- **13:** leave-one-out. Drop each BMD sample (11 runs) and each control (6 runs). Does the top-family list hold? Also rerun without GSM2934804.
- **14:** subsample DMD to 11, 100 draws, to show whether BMD’s weaker signal is only a sample-size effect.
- **15:** pig overlap.
    - Map Reschke’s microRNAs to human names, then to families.
    - Expected direction: pig microRNA **up** → our BMD−Ctrl NES **< 0**.
    - Report an agreement table, descriptive only. Caveats: pig, 3.5 months old, an engineered deletion.

---

## 7. Schedule (compressed, as of Sun 1 Nov night)

| When | Analysis (Aziz plus you) | Hamza |
| --- | --- | --- |
| **Sun night** | Steps 0–2 | Sample sheet; microRNA inputs; fibrosis sets |
| **Mon** | AM: steps 3–4. PM: steps 5–6. **Go/no-go in the evening.** | Fact-checks; pig list; sentence fixes |
| **Tue–Wed** | Steps 7–12 | Fibroblast markers; Intro and Related work drafts |
| **Thu** | Buffer for anything that slipped; 13–15 only if Wed finished clean. **Results freeze tonight. No new analyses after this.** | Reviews the figures; Data section |
| **Fri–Sat** | Aziz writes Methods and Results from the logs. You only fix bugs, polish figures and check numbers. | Limitations; finishes his sections |
| **Sun 8** | Both: check every citation against its source, proofread, submit |  |

---

## 8. Figures

- Planned figures:
    - QC boxplot;
    - collagen sanity table or bar chart;
    - miR-29 enrichment plot in DMD;
    - NES across the three comparisons for the top families, as a table or heatmap;
    - the step 8 histograms;
    - the detection curve.
- Avoid volcano and other scatter plots unless Aziz asks. If one is needed, pair it with the shuffle histogram.
- For **every** figure, add 3–5 lines to `docs/figure_notes.md`: what each axis is, one worked example point, and which pattern would mean what.

---

## 9. Logs, reproducibility, claims

- **`docs/methods_log.md`:** every parameter, default, package version, data source and version, cut, and deviation from the proposal, each with its reason.
- **`docs/results_log.md`:** every number Aziz might quote, each with the script and output file that produced it.
- **`docs/handoffs.md`:** the auto-built inputs (§4) and their status.
- **Results come only from scripts in this repo that actually ran.** Never describe an analysis as done unless it ran. Never invent a number, citation or claim. If something can’t be verified, write `[UNVERIFIED]` and say why.
- **Claim language:** “predicted”, “candidate”, “may contribute”. We infer microRNA activity; we don’t measure it, and nothing here shows causation or a validated biomarker. With 11 BMD samples, every result is hypothesis-generating.
- Keep raw data out of git, and keep intermediate `.rds` files so steps can be rerun separately.

---

## 10. Decisions to confirm with Aziz (use the default until he says otherwise)

| Decision | Default |
| --- | --- |
| Probe collapse rule | Highest mean expression |
| Low-expression gene filter | None beyond collapsing; flag if the ranks look noisy |
| TargetScan target set | Default (conserved) predictions, human |
| Muscle-expression source | The first credible published atlas you find, cited |
| Family size limits | At least 15 targets, no maximum |
| Smoke-test thresholds | As in step 6 |
| Shuffle count | 1,000 (500 if too slow) |
| Headline significance | Report fgsea padj **and** the label-shuffle p from step 8. fgsea shuffles genes, not people, and correlated genes can make its p-values look stronger than they are. |
| Catch criterion (detection limit) | NES > 0 and padj < 0.05 |

---

## 11. Literature facts you may rely on (checked) and ones you may not

- **Checked: Eisenberg et al. 2007 (PNAS).**
    - 6 BMD samples; only miR-146b and miR-221 passed their fold-change cutoff in BMD.
    - microRNA–gene links were significant only in DMD and Miyoshi myopathy.
    - In DMD, miR-29c was about 6-fold down, with COL3A1 about 6-fold up and COL1A2 about 5-fold up. The units aren’t stated, so say “about”.
- **Unchecked or wrong (don’t use until `docs/fact_check.md` clears them):**
    - Eisenberg “5 BMD” and “GSE3307 = same dataset” (not found);
    - Xie 2024 gene counts (not in the abstract);
    - the “13 genes” claim (undefined);
    - Dadgar 2014 as the source of GSE109178;
    - the pig “18 upregulated” (the directions are unverified);
    - GSM2934804’s severity;
    - GSE157668’s contents.