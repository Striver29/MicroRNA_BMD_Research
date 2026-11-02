# Handoff — where we are, what to do next

> **Updated end of session 2 (Mon 2 Nov, ~01:25 Tue 3 Nov).** See `## SESSION 2 STATUS` at the
> top of this file for the current position; the session-1 notes below it are still accurate
> for environment and setup.

## SESSION 2 STATUS (read this first)

### Compute left running overnight

When session 2 paused, a chained job was **still running unattended**. It runs three steps back
to back and writes its progress to `results/tables/queue_log.txt`:

| Step | Status at pause | ETA |
|---|---|---|
| 9 — direct BMD vs DMD, 500 patient shuffles | **58% done** | ~02:00 |
| 10b — detection limit, 600 simulations | queued | ~03:50 |
| 8 — shared-control baseline, 500 shuffles | queued | ~05:50 |

**First thing to do tomorrow:**

```bash
cd /Users/abdulazizalshujaa/Documents/MicroRNA_BMD_Research
cat results/tables/queue_log.txt          # did it reach "ALL DONE"?
ls -la results/rds/bmd_vs_dmd.rds results/rds/detection_limit.rds results/rds/shuffle_baseline.rds
```

- If `queue_log.txt` ends with **ALL DONE** and all three `.rds` files exist → everything
  finished. Go straight to writing up.
- If it stopped partway (laptop slept, etc.) → each step writes its outputs only at the very
  end, so there are **no corrupt partial files**. Just re-run whichever steps are missing:

```bash
export PATH="/opt/homebrew/bin:$PATH"
N_SHUFFLE=500 NPERM_NULL=10000 N_CORE=7 Rscript scripts/09_bmd_vs_dmd.R      2>&1 | tee results/tables/09_log.txt
N_REP=100     NPERM=10000      N_CORE=7 Rscript scripts/10b_detection_limit.R 2>&1 | tee results/tables/10b_log.txt
N_SHUFFLE=500 NPERM_NULL=10000 N_CORE=7 Rscript scripts/08_shuffle_baseline.R 2>&1 | tee results/tables/08_shuffle_log.txt
```

**To stop it instead:** `pkill -f "Rscript scripts/"`

### Still owed to Aziz (not yet done)

1. The full Notion page refresh with steps 8–14 results (`docs/NOTION_session_01.md` currently
   covers only steps 0–7, 10a, 12).
2. **Update `CLAUDE.md`** with the decisions taken, the environment realities and the
   performance lessons.
3. **Push to GitHub.** Push access is verified (dry run clean). 3 commits are currently unpushed:
   `git push origin main`.

### What session 2 settled

| Step | Result |
|---|---|
| **11** — 3′UTR-matched null | Survives at **6.5 SD** (p = 0.000999, the 1/1001 floor). **But a UTR-matched random set already scores NES ≈ 1.6**, so ~72% of miR-29's apparent 2.22 is UTR-length bias. **Never quote NES against a zero baseline.** All 25 top families survived. |
| **13** — leave-one-out | miR-29 significant in **17/17** runs; **rank 1 of 128 in every BMD-vs-DMD run**. Closes the step-2 QC question (GSM2934803 and GSM2934813 are harmless) and the GSM2934804 disputed-severity item. |
| **14** — DMD cut to n = 11 | **Negative result.** 53% of DMD-at-11 draws fall at or below BMD's miR-29 NES, so **measured against controls, BMD is NOT detectably weaker than DMD** at matched n. NES is a rank measure and saturates. The "milder" claim may rest only on the **direct** contrast and the raw fold-changes (COL1A1 +2.52 BMD vs +4.06 DMD). |

### Hardware ceiling (don't try to tune past it)

This Mac has **8 GB of RAM**. Measured parallel efficiency is **~65%** (4.5 effective cores of
7) for these jobs — memory bandwidth, not a bug. Earlier readings suggesting 20% efficiency
were lifetime averages distorted by slow starts; a 60-second delta is the honest measure:

```bash
PIDS=$(pgrep -f "08_shuffle_baseline" | tr '\n' ',' | sed 's/,$//')
cpu() { ps -o time= -p $PIDS | tr -d ' ' | awk -F: '{s+=$1*60+$2} END {print s}'; }
A=$(cpu); sleep 60; B=$(cpu); echo "effective cores: $(echo "($B-$A)/60" | bc -l)"
```

---

**Session:** Sun 2026-11-01 evening → Mon 2026-11-02 early hours.
**Schedule position (CLAUDE.md §7):** Sunday night was meant to be steps 0–2. We finished
**steps 0–7, 10a and 12**, i.e. most of Monday and some of Tuesday–Wednesday as well.
Step 8 was started, hit a performance bug, and was stopped unfinished. It must be re-run.
**Draft still due Sunday 8 November.**

---

## 1. How to restart the machine

R is installed at `/opt/homebrew/bin` and is **not on the default PATH**. Every session must
start with:

```bash
export PATH="/opt/homebrew/bin:$PATH"
cd /Users/abdulazizalshujaa/Documents/MicroRNA_BMD_Research
```

Check it works: `Rscript --version` → should print `Rscript (R) version 4.6.2`.

Everything is cached on disk, so **no step below needs re-running to continue**:

| Cache | Holds |
|---|---|
| `data/raw/CEL/` | the 34 CEL files, already un-gzipped |
| `data/processed/rma.rds` | the 54,675 × 34 RMA matrix |
| `data/processed/plm_qc.rds` | the RLE/NUSE probe-level model stats |
| `data/processed/gene_matrix.rds` | the 20,823 × 34 gene matrix |
| `results/rds/limma.rds` | the fitted limma object + all three result tables |
| `results/rds/family_sets.rds` | the 128 microRNA family → target-gene sets |
| `results/rds/fgsea_all.rds` | fgsea for all three contrasts |
| `inputs/auto/*.tsv` | the six rebuilt input files |

---

## 2. Run this next (in this order)

**Step 8 produced no output and must be re-run from the start.** It was stopped partway after a
performance bug was found (details in §5 and in `docs/methods_log.md`). The bug is **already
fixed in the scripts**; nothing needs editing before you resume.

```bash
export PATH="/opt/homebrew/bin:$PATH"
cd /Users/abdulazizalshujaa/Documents/MicroRNA_BMD_Research
```

### 2a. FIRST: measure real throughput before committing to a long run

The runtime figure the scripts print is **not trustworthy** — it is derived from a serial
timing loop and does not reflect parallel throughput. Measure the real thing with a short
parallel batch and extrapolate:

```bash
N_SHUFFLE=21 NPERM_NULL=10000 N_CORE=7 Rscript scripts/08_shuffle_baseline.R 2>&1 | tail -20
```

That runs 21 shuffles (3 per worker) and prints `Done in X min`. Multiply by 1000/21 to get the
true cost of the full run. **Only then** decide between 1,000 shuffles and the 500 that
CLAUDE.md §10 permits as a fallback.

> Note: that command writes a 21-shuffle `results/rds/shuffle_baseline.rds` and overwrites the
> step-8 figures with meaningless versions. That is fine — the real run overwrites them again.
> Just don't quote numbers from the timing run.

### 2b. Then run the remaining steps, ONE AT A TIME

```bash
N_SHUFFLE=1000 NPERM_NULL=10000 N_CORE=7 Rscript scripts/08_shuffle_baseline.R 2>&1 | tee results/tables/08_shuffle_log.txt
for job in 11_utr_matched 10b_detection_limit 09_bmd_vs_dmd 13_leave_one_out 14_subsample_dmd; do
  N_CORE=7 Rscript "scripts/${job}.R" 2>&1 | tee "results/tables/${job}_log.txt"
done
```

| Step | What it answers | Runtime |
|---|---|---|
| **8** | Does the BMD/DMD agreement beat the shared-control baseline? | **measure first** (§2a) |
| **11** | Is the footprint just a long-3′UTR artefact? | unmeasured since the fix |
| **10b** | How big must a footprint be to be visible at n = 11? | unmeasured since the fix |
| **9** | Does BMD genuinely differ from DMD? (the clean test) | unmeasured since the fix |
| **13** | Does any single sample carry the result? **Settles the two QC-flagged samples.** | fast (17 re-runs) |
| **14** | Is BMD's weaker signal just a sample-size effect? | fast (100 draws) |

**Never run two of these at once.** Session 1 lost ~50 minutes to exactly that: step 10b was
given 3 cores while step 8 held 7 on an 8-core machine, and its projected runtime blew up from
~30 minutes to 11.5 hours.

**Step 15 (pig overlap) stays BLOCKED.** It needs Reschke et al. 2026, *J Cachexia Sarcopenia
Muscle*, Table 1. CLAUDE.md §4 forbids building `inputs/pig_mirnas.tsv` from any other source.
If Hamza supplies the table, build the file and run step 15; otherwise the paper states the
comparison was not possible.

## 3. The scientific state of play

### The result, in one paragraph

The method works: in DMD vs control, miR-29's predicted footprint comes out at **NES = +2.22,
padj = 5.4 × 10⁻³⁵**, rank 3 of 128 families, in the direction the literature predicts (targets
up ⇒ brake weaker ⇒ miR-29 less active). But that alone proves little, because **99 of 128
families point the same way** — a competitive rank test does that when the whole transcriptome
shifts. The finding that actually carries the paper is the **direct BMD-vs-DMD comparison**,
where the shared controls cancel out: there the global skew collapses from 103 significant
families to **7**, and **miR-29 is rank 1 at padj = 4.7 × 10⁻¹⁰**, about five orders of
magnitude clear of the next family. Its sign (NES = −1.59) says the miR-29 brake is **less
weakened in BMD than in DMD** — milder, in the predicted direction.

### Three independent specificity signals, none of which is a p-value

1. **Leading-edge content.** miR-29's leading edge is **6.5% fibrosis-set genes**, the highest
   of any top family (next 4.1%, most 1–3%), and the gene list is a collagen roll-call:
   COL1A1, COL1A2, COL3A1, COL4A1/2/4, COL5A1/2/3, COL6A2/3/6, COL15A1, COL16A1, COL21A1,
   COL27A1, ELN, FBN1, LAMA2, LAMC1, NID1, SPARC.
2. **Survives deleting the scar genes** (step 10a). Removing all 275 matrisome-core genes from
   both the ranked list and the family sets leaves miR-29 at **96.7% of its NES**
   (+2.224 → +2.151, padj 1.4 × 10⁻²⁹). The footprint is not just scar.
3. **Survives adjusting for fibroblast content** (step 12). Adding a fibroblast-abundance
   covariate halves the significant families (103 → 51), but miR-29 keeps **79%** of its NES
   where other top families keep ~55–60%, and afterwards it leads by eleven orders of magnitude
   (padj 1.7 × 10⁻¹⁶ vs 3.1 × 10⁻⁵ for the next family).

### What is still unknown

- **Step 8 (started, stopped unfinished, NO result).** Whether the BMD/DMD agreement
  (**r = 0.9123**, slope **0.9769**) beats the shared-control baseline. Until this is re-run
  and lands, **do not claim BMD is a weaker copy of DMD**; the spec's own simulation says noise
  alone yields r ≈ 0.7 here.
- **Step 9.** Whether miR-29's BMD-vs-DMD difference survives shuffling patient labels.
- **Step 11.** Whether the footprint survives a 3′UTR-length-matched null.
- **Step 10b.** The detection limit at n = 11 — needed to interpret every null result.
- **Steps 13–14.** Sample-level robustness.

---

## 4. Decisions already taken (do not reopen without Aziz)

| Decision | Outcome |
|---|---|
| **QC STOP point (§6 step 2)** | **CLOSED — keep all 34.** RLE/NUSE clean on every array. The spec's correlation metric was confounded by group size and is not used for the decision. Two samples (GSM2934803, GSM2934813) remain mildly flagged within group and are kept; step 13 tests them. Full reasoning in `docs/methods_log.md`. |
| **Smoke-test go/no-go (§6 step 6)** | **GO.** Still needs Hamza's joint sign-off per the spec, but the criterion (NES > 0, padj < 0.05) is met with enormous margin. |
| Muscle-expression source | miRNA Tissue Atlas 2025 (CCB Saarland), limb skeletal muscle, **mean RPMM ≥ 10**. The threshold is ours, not the spec's — logged as a sensitivity-test candidate. |
| fgsea permutations | `nPermSimple = 10000` everywhere (the 1,000 default left 13 of 128 families with NA NES). |

---

## 5. Known traps for whoever picks this up

1. **Never nest `fgsea` inside `bplapply` without `nproc = 1`.** `fgseaMultilevel` defaults to
   `nproc = 0` = "use every core", so 7 workers each spawned a full-width fgsea: up to 56
   threads on 8 cores. Step 8 spent **52% of CPU in kernel context-switching** and completed
   ~19% of its work in 99 minutes (projected finish ~6.8 h). Fixed in scripts 08, 09, 10b, 11,
   13, 14. Left alone in 06, 07, 10a, 12, where fgsea is top-level and full width is correct.
   **Diagnostic recipe:** if `top` shows system CPU anywhere near user CPU, or `ps -o time=`
   on the workers shows far less CPU time than wall time, you are thrashing, not computing.
2. **`setorder(dt, -abs(x))` does not work in data.table.** It crashed step 11 after it had
   done the work, and the same line was sitting in step 9 where it would have destroyed a
   45-minute run. Fixed in both. Always `dt[, absX := abs(x)]; setorder(dt, p, -absX)`.
3. **Never judge r or the slope against 0.** Only against the step-8 shuffled baseline.
   Any scatter plot must be shown beside the shuffle histogram (CLAUDE.md §0, §8).
4. **NCBI FTP returns 403 from this machine.** Use GEO's web endpoints (already wired into
   `scripts/01_download_labels.R`). `GEOquery` is not installed and is not needed.
5. **Disk was down to 9 GB at session start.** The 795 MB TargetScan UTR archive is downloaded,
   filtered to human, and deleted by the build script. Don't re-add it.
6. **Claim language (CLAUDE.md §9):** "predicted", "candidate", "may contribute". We infer
   microRNA activity from target behaviour; we never measured a microRNA. With 11 BMD samples
   everything is hypothesis-generating.

---

## 6. Two confounds that must appear in Limitations

1. **Control ages are absent from GEO entirely.** All 6 control records carry `age: NA` and no
   sex. CLAUDE.md §2 asked us to record them; they do not exist in the source. The age confound
   therefore **cannot be resolved** — state it, don't fix it.
2. **BMD patients are about a decade older than DMD patients** (median 12.2 y, range 2.5–60.8
   vs median 3.5 y, range 0.9–8.0). The direct BMD-vs-DMD comparison — the one carrying the
   headline — therefore confounds disease with age. It cannot be adjusted away, because
   adjusting for age would also remove the disease signal.

---

## 7. Where everything lives

| Path | Contents |
|---|---|
| `scripts/00–14` | one script per pipeline step, each runnable standalone |
| `scripts/build_inputs.R` + `build_inputs_muscle.py` | rebuild the six auto inputs |
| `docs/methods_log.md` | every parameter, default, version, deviation, cut — with reasons |
| `docs/results_log.md` | every number Aziz might quote, with its script and output file |
| `docs/figure_notes.md` | per figure: axes, a worked example point, what each pattern means |
| `docs/fact_check.md` | what is confirmed, what is `[UNVERIFIED]`, and why |
| `docs/handoffs.md` | the auto-built inputs and their sources/versions, for Hamza to diff |
| `results/tables/`, `results/figures/`, `results/rds/` | outputs |
