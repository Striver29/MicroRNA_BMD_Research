# Step 10b — detection limit. Handoff to Hamza

> **✅ DONE, 3 Nov 2026. Do not rerun.** Ran in full (6 × 100, NPERM 10000) on a 16-core
> Linux box in 76.2 min, with package versions identical to the Mac runs. All five deliverables
> in §4 are in the repo. **The f = 0 false-alarm check failed (1.00)**, so a calibrated
> re-scoring was added (`scripts/10b2_calibrated_catch.R`). See `docs/results_log.md`, step 10b.

**What's needed:** one R script run to completion, and four output files sent back.
**Why it moved to you:** it needs ~15.5 CPU-hours. Aziz's machine is a fanless MacBook Air M1
(4 real performance cores, 8 GB RAM) and it kept sleeping; after 10 hours it was 64% done.
On any desktop with 8+ real cores this is a **2–4 hour** job.

**Deadline: all computation frozen Friday 6 November.** This is the last outstanding analysis —
everything else in the pipeline is finished.

---

## 1. What this analysis is for (the short version)

Our main result is a **null**: after testing properly, no microRNA family shows a detectable
difference between BMD and DMD, and miR-29's footprint versus healthy controls sits at
p = 0.052, just outside significance.

A null result is only meaningful if you can say **how big an effect you would have caught.**
Otherwise "we found nothing" is indistinguishable from "our method is blind."

Step 10b answers exactly that. It **plants a fake microRNA footprint of known size** into the 11
BMD samples, then checks how often our own pipeline detects it. Repeat across a range of effect
sizes and you get a detection curve. Where that curve crosses 80% is our detection limit.

This is now the **central result of the paper**, not a side check.

### How the planting works

1. **The ruler.** Measure the real DMD miR-29 footprint: mean log2 fold-change of the 1,158
   miR-29 target genes minus that of all other genes. It comes out at **+0.2686 log2 units**
   (a 20.5% rise). The script recomputes this; it should print that number.
2. **The plant.** Pick a random gene set the same size as miR-29's, **matched on 3′UTR length
   decile by decile** (long-UTR genes behave alike, so an unmatched random set would be unfair).
   Add `f × 0.2686` to those genes **in the 11 BMD samples only**.
3. **The test.** Re-run BMD-vs-control limma and fgsea, with the planted set thrown in among all
   128 real families. A "catch" is NES > 0 and padj < 0.05.
4. **The grid.** f ∈ {0, 10, 25, 50, 75, 100}%, 100 repeats each = 600 simulations.
   **At f = 0 nothing is planted, so the catch rate there is our false-alarm rate** — it should
   land near 0.05 or below. If it doesn't, something is wrong and we need to know.

---

## 2. Setup

```bash
git clone https://github.com/Striver29/MicroRNA_BMD_Research.git
cd MicroRNA_BMD_Research
Rscript scripts/00_setup.R        # installs everything; slow on first run
```

**You do NOT need to download any CEL files.** Every input this script needs is already in the
repo:

| File | What it is |
|---|---|
| `data/processed/gene_matrix.rds` | the 20,823 × 34 gene expression matrix (5 MB) |
| `results/rds/family_sets.rds` | the 128 microRNA family → target-gene sets |
| `results/rds/limma.rds` | fitted limma object, used to compute the ruler |
| `inputs/auto/utr_lengths.tsv` | 3′UTR lengths, for the matched planting |

Sanity check before you start — this should print `TRUE` four times:

```bash
Rscript -e 'for (f in c("data/processed/gene_matrix.rds","results/rds/family_sets.rds",
  "results/rds/limma.rds","inputs/auto/utr_lengths.tsv")) cat(f, file.exists(f), "\n")'
```

---

## 3. Run it

```bash
# First: time it without committing. Prints a projection, then exits.
TIME_ONLY=1 N_CORE=8 Rscript scripts/10b_detection_limit.R

# Then the real run. Set N_CORE to your PHYSICAL core count, not hyperthreads.
N_REP=100 NPERM=10000 N_CORE=8 Rscript scripts/10b_detection_limit.R 2>&1 | tee results/tables/10b_log.txt
```

**Scale:** 600 runs × ~93 s each = **15.5 CPU-hours**. On 8 real cores with ~80% efficiency,
roughly **2.5 hours**. Memory is modest (~1 GB per worker).

### Three things that will bite you

1. **Do not raise `N_CORE` above your physical cores.** fgsea is memory-bandwidth-bound; we
   measured only 65% efficiency at 7 workers on 4 real cores. Hyperthreads make it worse.
2. **Do not run anything else heavy at the same time.** We lost ~50 minutes to two R jobs
   competing for the same 8 cores.
3. **`nproc = 1` is already set inside the script and must stay.** `fgseaMultilevel` defaults to
   grabbing every core, so nesting it inside `bplapply` spawns workers × cores threads. That
   cost us a 98-minute job turning into a projected 6.8 hours, with 52% of CPU lost to kernel
   context-switching. If you see system CPU anywhere near user CPU, that's the symptom.

### If it is still too slow

Drop to `N_REP=50`. Catch-rate precision goes from about ±5% to ±7%, which does not move where
the curve crosses 80%. **Tell us if you do this** so it goes in the methods log.

---

## 4. What to send back

| File | Why |
|---|---|
| `results/tables/10b_detection_curve.tsv` | the curve itself — the headline numbers |
| `results/tables/10b_detection_runs.tsv` | all 600 individual simulations |
| `results/figures/10b_detection_curve.png` | the figure for the paper |
| `results/rds/detection_limit.rds` | everything, for re-plotting |
| `results/tables/10b_log.txt` | the console log, including the printed ruler |

A pull request or a zip both work.

### What the result should look like

The script prints a table like this:

```
 f_pct  planted_log2  planted_pct_linear  mean_NES  catch_rate   n
     0        0.0000                 0.0      ~0.0        0.0x  100
    10        0.0269                 1.9       ...         ...  100
    25        0.0672                 4.8       ...         ...  100
    50        0.1343                 9.7       ...         ...  100
    75        0.2015                15.0       ...         ...  100
   100        0.2686                20.5       ...         ...  100
```

**Two things we care about most:**
- **the f = 0 row** — this is the false-alarm rate and should be ≈ 0.05 or lower;
- **where catch_rate crosses 0.80** — that is the detection limit we will quote.

If the curve never reaches 80% even at f = 100%, that is a legitimate and important finding:
it would mean that with 11 BMD samples we could not reliably detect even an effect the size of
the real DMD miR-29 footprint. Don't treat that as a failure — report it.

---

## 5. Context you may want

Full background is in `docs/NOTION_session_01.md` (written for a non-specialist) and
`docs/results_log.md` (every number with the script that produced it). Two findings worth
knowing because they shape the paper:

- **fgsea's p-values are overconfident by up to 34 orders of magnitude** on this data. It
  permutes genes, not people. miR-29 DMD-vs-control: fgsea padj 1.5e-34, patient-label-shuffle
  p = 0.052. This is why step 10b matters.
- **GSE38417 is a duplicate re-deposit of GSE109178** — identical raw arrays, verified across
  all 1,354,896 probe cells. Before merging any GEO accession, run
  `scripts/check_duplicate_arrays.R`. Matching metadata is a red flag; identical intensities
  are proof.
