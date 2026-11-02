# Figure notes

Per CLAUDE.md §8: for every figure, what each axis is, one worked example point, and which
pattern would mean what.

---

## `results/figures/02_qc_boxplot.png` — QC boxplot

- **x:** the 34 samples, grouped BMD / control / DMD. **y:** log2 RMA expression; one box
  summarises all 54,675 probesets on that array.
- **Worked example:** a box whose middle line sits at 6.2 means half that array's probesets are
  below 2^6.2 ≈ 73 units and half above.
- **What to look for:** after RMA the boxes should line up, because quantile normalisation
  forces every array to the same distribution. **They do.** A box shifted up or down, or much
  taller or shorter than the rest, would mean normalisation failed on that array.

## `results/figures/02_qc_within_group_cor.png` — array agreement within group

- **x:** median Pearson correlation of that sample with the other samples **in its own group**.
  **y:** sample.
- **Worked example:** GSM2934803 scores 0.9676. Its 10 fellow BMD arrays correlate with it
  between roughly 0.95 and 0.98, and 0.9676 is the middle of those 10 numbers.
- **Why within group and not against all 34:** 28 of 34 arrays are dystrophic muscle and only 6
  are healthy. Measured against all 34, every control looks like an outlier simply because few
  arrays share its biology — all 6 controls scored the lowest 6 values (0.9496–0.9558 against a
  0.9759 centre). Within group they score 0.9784–0.9854, i.e. the **highest**. The "outlier
  controls" were an artefact of the metric.
- **What would mean trouble:** one sample far below its own group, **together with** a bad RLE
  or NUSE. Here no array has either.

## `results/figures/02_qc_nuse.png` — per-array noise (NUSE)

- **x:** median NUSE (Normalised Unscaled Standard Error), scaled so the typical array here is
  1.00. **y:** sample. Red dashed line = 1.05, the usual flag threshold.
- **Worked example:** GSM2934807 has NUSE 1.039, so its probe-level fits are about 4% noisier
  than the typical array — well inside tolerance.
- **What to look for:** every bar left of 1.05. **All 34 are** (range 0.981–1.039). An array
  above 1.05 would be a genuine candidate for dropping.

## `results/figures/06_mir29_enrichment_DMDvsCtrl.png` — the ranked walk

This is the central picture of the smoke test. Read it slowly.

- **x:** all 20,823 genes, lined up left to right by their limma moderated *t* for DMD vs
  control. Far left = most up in DMD; far right = most down in DMD.
- **y:** the running enrichment score. Start at 0 at the far left. Walk right. Every time you
  step on one of the 1,158 miR-29 target genes the line goes **up**; every time you step on a
  non-target it drifts slightly **down**. The step sizes are set so that if targets were
  sprinkled evenly the line would wander around 0 and finish at 0.
- **Worked example with 10 genes and 3 targets (T):** order `T T T n n n n n n n`. Each T adds
  +1/3, each n subtracts 1/7. The line goes 0.33 → 0.67 → **1.00** → 0.86 → 0.71 → … → 0.
  The peak is 1.00 and it comes early. Now order them `n n n n n n n T T T`: the line goes
  −0.14 → −0.29 → … → −1.00 at gene 7, then climbs back to 0. Peak height is the same
  magnitude but it is **negative** and comes late. NES is that peak, normalised so sets of
  different sizes are comparable.
- **Our result:** NES = **+2.22**, padj = 5.4 × 10⁻³⁵. The curve climbs steeply near the left,
  i.e. miR-29's targets are crowded into the up-in-DMD end.
- **Which pattern means what (CLAUDE.md §1 direction rule):**
  - peak early and positive → targets up → the brake is **weaker** → miR-29 predicted **less
    active** in DMD. ← **this is what we see, and it is the published direction**
  - peak late and negative → targets down → brake **stronger** → predicted **more active**
  - line wandering near 0 → no footprint.
- **Honesty note to carry into the text:** this plot alone does not show that miR-29 is
  *special*. 99 of 128 families also came out NES > 0 in this comparison. The specificity
  evidence is elsewhere: step 9 (direct BMD vs DMD), step 10a, step 11, step 12.

## `results/figures/07_nes_heatmap.png` — NES across the three comparisons

- **rows:** the top 25 families. **columns:** BMD−Ctrl, DMD−Ctrl, BMD−DMD. **colour:** NES,
  red positive, blue negative, white 0. The number is printed in each cell.
- **Worked example:** the miR-29 row reads +2.17 (BMD−Ctrl), +2.22 (DMD−Ctrl), −1.59 (BMD−DMD).
  First two: targets up versus healthy muscle in **both** diseases, so the brake is weaker in
  both. Third: targets **lower in BMD than in DMD**, so the brake is less weakened in BMD than
  in DMD — milder, in the predicted direction.
- **What to look for:** the first two columns are almost uniformly red for every family (the
  global skew). The third column is where families separate: most cells there are not
  significant, and miR-29's is the strongest.

## `results/figures/08_shuffle_r.png` and `08_shuffle_slope.png` — the shared-control baseline

**Always show these next to any claim that BMD looks like a weaker copy of DMD.**

- **x:** the value of r (or of the slope) obtained after shuffling all 34 labels at random while
  keeping 6 control / 11 BMD / 17 DMD. **y:** how many of the 1,000 shuffles landed there.
  Red line = the value we actually measured.
- **Why the grey histogram is not centred on 0:** both comparisons are measured against the
  **same 6 controls**. CLAUDE.md §5 gives the toy case: two genes with no true difference, where
  a fluke in the controls pushes BMD−Ctrl to +0.35 and DMD−Ctrl to +0.28 for one gene and to
  −0.34 / −0.28 for the other. The two columns move together because the controls are shared,
  not because BMD and DMD share biology. With these group sizes that alone produces r ≈ 0.7
  with a slope below 1 — exactly the "weaker copy" picture, from pure noise.
- **Worked example of reading the p:** if 3 of 1,000 shuffles reached our r, the empirical
  p is (3 + 1) / (1,000 + 1) = 0.004.
- **Which pattern means what:**
  - red line far to the right of the grey bulk → the agreement is more than the shared controls
    can manufacture → BMD and DMD really do share footprint structure;
  - red line inside the grey bulk → **our r means nothing**, and the "weaker copy" reading must
    be dropped. In that case the only trustworthy answer to "is BMD different from DMD?" is
    step 9.

## `results/figures/09_mir29_bmd_vs_dmd_shuffle.png` — miR-29, controls removed

- **x:** miR-29's NES for BMD vs DMD after shuffling the BMD/DMD labels among the **28 patients
  only** (keeping 11/17). **y:** count over 1,000 shuffles. Red = measured value.
- **Why this test is the clean one:** the 6 controls never enter it, so the shared-control trap
  cannot operate at all. This is what answers "is BMD different from DMD?"
- **What to look for:** red line outside the grey bulk → BMD and DMD genuinely differ for
  miR-29. Red inside → the apparent difference is reshuffling noise.
- **Confound to state regardless of the result:** the BMD group is about a decade older than the
  DMD group (medians 12.2 y vs 3.5 y). A BMD-vs-DMD difference may be partly age.

## `results/figures/10b_detection_curve.png` — the detection limit

- **x:** the size of a footprint we deliberately planted, as a percentage of the real DMD miR-29
  footprint. **y:** the share of 100 simulated runs in which our own pipeline caught it
  (NES > 0 and padj < 0.05).
- **Worked example:** the ruler is +0.269 log2 units. At x = 50 we add 0.5 × 0.269 = +0.135 to
  each of ~1,158 planted genes, **in the 11 BMD samples only** — about a 10% rise on the linear
  scale. If 62 of 100 runs then flag the planted set, the point is at y = 0.62.
- **x = 0 is the control:** nothing is planted, so anything caught there is a false alarm. That
  point should sit near or below 0.05.
- **Which pattern means what:** the detection limit is where the curve crosses **80%**. If it
  crosses at 50%, then with 11 BMD samples we can only reliably find footprints at least half
  the size of DMD's miR-29 footprint — and a null result for a smaller footprint means
  "too small to see here", not "absent".

## `results/figures/12_fibroblast_score.png` — how much fibroblast is in each biopsy

- **x:** group. **y:** the fibroblast score = mean row z-score of 316 fibroblast markers that
  are **not** themselves fibrosis-set genes.
- **How a row z-score works, with 2 genes × 3 samples:** gene A reads 8, 9, 10 → mean 9, sd 1 →
  z = −1, 0, +1. Gene B reads 5, 6, 7 → z = −1, 0, +1. Sample scores: −1, 0, +1. So the score
  says "high or low on fibroblast genes **relative to the other samples here**"; it is not an
  absolute cell fraction, and every marker gets equal weight.
- **Our values:** control −0.79, BMD +0.04, DMD +0.25 (ANOVA p = 1.2 × 10⁻⁸). BMD sits between
  healthy and DMD, as the disease picture predicts.
- **Caveat that must travel with this figure (CLAUDE.md §12):** because scar load really does
  differ by group, putting this score in the model also removes some genuine disease signal.
  Read the adjusted result as "does the footprint survive", never as the headline estimate.

## `results/figures/13_mir29_leave_one_out.png` — does one sample carry the result?

- **y:** which sample was removed (red = a BMD patient, teal = a control). **x:** miR-29's NES
  when the analysis is re-run without that sample. Two panels: BMD-vs-control and BMD-vs-DMD.
  Red dashed line = the full-data value.
- **Worked example.** Removing GSM2934812 gives NES +2.120 in BMD-vs-control. The full-data
  value is +2.169. So that sample was worth 0.049 of NES — about 2% of the signal.
- **What to look for:** every bar reaching the dashed line. **They do** — the 17 re-runs span
  only +2.120 to +2.186 (BMD-vs-control) and −1.674 to −1.500 (BMD-vs-DMD).
- **What would mean trouble:** one noticeably short bar. That would say "this result exists
  because of one patient", which at n = 11 would be fatal. Nothing like that appears.
- **Why this figure settles the QC argument from step 2:** GSM2934803 and GSM2934813 were the
  two samples flagged as mild within-group outliers. Their bars are indistinguishable from the
  rest, so keeping them changed nothing.

## `results/figures/14_subsample_dmd.png` — is BMD only weaker because n is smaller?

- **x:** miR-29's NES for DMD vs control, after randomly cutting DMD from 17 patients down to
  11 — the same number BMD has. **y:** how many of 100 draws landed there.
  **Red solid** = the real BMD value (2.169). **Blue dashed** = full 17-sample DMD (2.224).
- **Worked example of the comparison being made.** Draw 1 picks 11 random DMD patients and
  scores 2.19; draw 2 picks a different 11 and scores 2.13; and so on 100 times. The grey
  histogram is the answer to "what does DMD look like when it only has 11 patients?"
- **What we see:** the red line sits **in the middle of the grey bulk** — 53% of draws fall at
  or below it. The blue dashed line sits at the upper edge, which is expected: more samples
  give cleaner statistics and a slightly higher score.
- **Which pattern means what:**
  - red line *inside* the grey bulk → BMD is within what DMD itself produces at n = 11, so the
    vs-control difference is **sample size, not biology**. ← **this is what we see**
  - red line *left of* the grey bulk → BMD really is weaker than DMD even at matched n.
- **The consequence, which must reach the paper:** we may **not** claim "BMD's miR-29 footprint
  is weaker than DMD's" from the vs-control scores. NES is a rank measure that saturates — it
  asks whether targets are bunched at the up end, not how far they moved. The magnitudes do
  differ (COL1A1 +2.52 in BMD vs +4.06 in DMD), but only the **direct** BMD-vs-DMD contrast
  can see that.
