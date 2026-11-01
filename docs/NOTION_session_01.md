# BMD microRNA–Fibrosis Project — Session 1

**Date:** Sunday 1 November 2026, evening → early Monday 2 November
**Covered:** CLAUDE.md steps 0–7, 10a, 12 (step 8 started, not finished)
**Status:** smoke test **PASSED**. Headline result identified. Four robustness steps still to run.

> **How to read this page.** Nothing below assumes biology or statistics. Every technical
> term is defined where it first appears, and every statistical claim comes with a worked
> example using real numbers from our data.

---

## Part 1 — The question, in plain words

### What the disease is

**Duchenne muscular dystrophy (DMD)** is caused by a broken *dystrophin* gene. Dystrophin is a
shock absorber inside muscle fibres. Without it, fibres tear every time the muscle contracts,
die, and are gradually replaced by **collagen scar tissue**. That replacement process is called
**fibrosis**. Boys with DMD usually lose the ability to walk in their early teens.

**Becker muscular dystrophy (BMD)** is the *same broken gene*, but the break is less severe —
the body still makes a shortened, partly working dystrophin. Same disease process, slower and
milder. Some people with BMD stay mobile into their sixties. Our dataset contains one BMD
patient aged 60.8.

### Why anyone should care about the comparison

The leading DMD therapy, **exon skipping**, works by converting a DMD-type mutation into a
BMD-type one. It deliberately turns DMD into something closer to BMD. So the question *"what
does the molecular biology of BMD actually look like?"* is the question *"what is exon-skipping
therapy aiming at?"*

### What a microRNA is

Think of a gene as an instruction that gets copied into a working message called mRNA, which
the cell then reads to build a protein.

A **microRNA** is a very short piece of RNA — about 22 letters — that acts as a **brake**. It
finds messages and shuts them down, so less protein gets built.

The targeting works through a 7-letter stretch called the **seed** (letters 2–8 of the
microRNA). The seed pairs with a matching 7-letter site in the tail end of a message — a region
called the **3′UTR**. Match found, message destroyed or blocked.

One microRNA has one seed, so it brakes **every** message carrying that matching site —
often hundreds or thousands of genes at once.

> **Worked example.** miR-29 has the seed `AGCACCA`. In our data, **1,158** genes present on the
> array carry a matching site. miR-29 is one brake pressing on all 1,158 at the same time.
> Many of them are collagen genes — the scar proteins.

### The known fact we are building on

In **DMD**, miR-29 levels go **down**. The brake on collagen is released, collagen piles up,
scarring accelerates. This is established.

**Nobody has tested whether the same thing happens in human BMD.** That gap is our paper.

### The problem we had to work around

**Our dataset contains no microRNA measurements at all.** It only measures genes.

So we cannot look at miR-29 directly. Instead we look for its **footprint**.

> **The footprint idea, by analogy.** You cannot see the wind, but you can see which way the
> trees lean. If miR-29's brake weakens, then *all 1,158 of its target genes should drift
> upward together* — a coordinated nudge, too consistent to be chance. We detect the brake by
> the shadow it leaves on its targets.

### The direction rule — memorise this, it is easy to get backwards

The output of our main test is a number called **NES** (explained in Part 3).

| What we see | What it means |
|---|---|
| **NES > 0** — targets bunched at the **UP** end | brake is **weaker** → microRNA **LESS active** |
| **NES < 0** — targets bunched at the **DOWN** end | brake is **stronger** → microRNA **MORE active** |

It feels inverted because a brake and the thing it brakes move in *opposite* directions.
Brake off → targets up. So a microRNA that genuinely goes **up** should give **NES < 0**.

### The core question

Is BMD's microRNA footprint **(a)** a weaker copy of DMD's, **(b)** a different pattern
altogether, or **(c)** nothing we can detect with only 11 BMD patients?

---

## Part 2 — What we had to build before any science could happen

### R was not installed

The entire analysis plan is built on R and Bioconductor. Neither was on the machine.

Installed R 4.6.2 and Bioconductor 3.23. Every package had to be **compiled from source**
(Homebrew's R cannot use the pre-built Mac binaries), which took most of an hour in the
background while other work continued. Full inventory in `results/sessionInfo.txt`.

`GEOquery` failed to compile. It turned out not to matter — see next.

### NCBI blocked our downloads

The data live at NCBI's **GEO** archive (accession **GSE109178**). The standard route is NCBI's
FTP server. From this machine it returns **HTTP 403 Forbidden** for every file.

Worked around using GEO's two official web endpoints, which serve identical content:

- metadata → `ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109178&targ=gsm&form=text&view=brief`
- raw data → `ncbi.nlm.nih.gov/geo/download/?acc=GSE109178&format=file` (217 MB)

Logged as a deviation in `docs/methods_log.md`. The parsed fields are exactly the ones the
standard R route would have exposed.

### Hamza's input files did not exist yet, so we built all six

CLAUDE.md §4 says: if an input is missing, build it yourself to the same format, put it in
`inputs/auto/`, log it, and keep moving. Hamza's versions will be diffed against ours later.

| File | What it is | Where it came from |
|---|---|---|
| `mirna_families.tsv` | 2,064 microRNA families with their 7-letter seeds | TargetScanHuman 8.0 |
| `mirna_targets.tsv` | 115,308 microRNA→gene targeting predictions | TargetScanHuman 8.0, conserved sites |
| `utr_lengths.tsv` | 3′UTR length for 19,429 genes | TargetScan UTR sequences, alignment gaps stripped |
| `fibrosis_sets.tsv` | 275 core scar-matrix genes + 54 TGF-β signalling genes | MSigDB (Naba et al. 2012; Hallmark) |
| `fibroblast_markers.tsv` | 372 markers of the scar-producing cell | MSigDB C8, Rubenstein human skeletal-muscle single-cell atlas |
| `sample_sheet.csv` | who is in which group | parsed from GEO |

**A judgement call you should know about.** CLAUDE.md asked us to filter to microRNAs actually
present in muscle, and to pick a published source ourselves. We used the **miRNA Tissue Atlas
2025** (Saarland University), taking its *limb skeletal muscle* samples — the same kind of
muscle as the thigh biopsies in our dataset — with a cutoff of **mean RPMM ≥ 10**.

> **That cutoff is ours, not the spec's.** It is logged as a parameter to sensitivity-test if
> anyone challenges it.

**A sanity check that the atlas data are sound:** the most abundant microRNAs it reports in
skeletal muscle are miR-486, miR-10b, **miR-133a, miR-1, miR-206**. The last three are the
textbook muscle-specific "myomiRs" — microRNAs famous for existing almost nowhere else. Their
appearing at the top is strong evidence the table is behaving correctly.

miR-29 in muscle: miR-29a-3p at **1,895 RPMM**, miR-29b-3p at 650, miR-29c-3p at 439. Clearly
present.

**Step 15 (pig comparison) is BLOCKED.** It requires one specific table from Reschke et al.
2026, which we could not obtain. CLAUDE.md explicitly forbids building that input from any
other source, so we did not invent one.

---

## Part 3 — The analysis, step by step

### Step 1 — Who is in the dataset

49 muscle biopsies, all from the **vastus lateralis** (a thigh muscle), all on the same
Affymetrix chip type.

| Group | n | Keep? |
|---|---|---|
| Healthy control | **6** | ✅ |
| BMD | **11** | ✅ |
| DMD | **17** | ✅ |
| Limb-girdle dystrophy (8 dysferlin + 7 FKRP) | 15 | ❌ different disease, dropped |

**This matched the specification exactly** — the script refuses to continue otherwise.

**Two confounds we found and cannot fix:**

1. **Control ages are simply not in GEO.** All 6 records say `age: NA`. We were asked to record
   them; they do not exist. The age confound is therefore unresolvable and goes in Limitations.
2. **BMD patients are about a decade older than DMD patients** — median **12.2 years**
   (range 2.5–60.8) versus **3.5 years** (range 0.9–8.0). DMD children are diagnosed and
   biopsied young; BMD is milder and often caught later. This means our BMD-vs-DMD comparison
   mixes *disease severity* with *age*. It cannot be adjusted away, because correcting for age
   would also delete the disease signal. **This must be stated in the paper.**

### Step 2 — Turning chip images into numbers, and checking the chips

Each chip has ~54,675 spots. Converting raw brightness into comparable numbers uses a standard
procedure called **RMA**, which does three things: removes background glow, forces all 34 chips
onto a common scale so they can be compared, and combines the ~11 spots measuring each gene
into one value.

Result: a **54,675 × 34** table.

**A QC trap worth understanding** (this cost real time and is a good story for the professor).

The plan said: flag any sample whose correlation with all the other samples is unusually low.
Correlation asks "do these two chips rank genes the same way?" — a broken chip agrees with
nobody.

Applied literally, it flagged **three of our six controls**. All six controls were the six
lowest-scoring samples.

That is not six bad chips. It is the metric breaking:

> **Worked example — GSM2934814, one of the flagged controls.**
>
> | Compared against | median correlation |
> |---|---|
> | the other **5 controls** | **0.9854** ← excellent |
> | the 11 BMD patients | 0.9620 |
> | the 17 DMD patients | 0.9437 |
> | **median over all 33** | **0.9496** ← the flagged number |
>
> The chip is fine. But 28 of its 33 comparisons are against *diseased* muscle, so the median
> lands inside the healthy-vs-diseased comparisons. With 28 diseased chips and only 6 healthy
> ones, **every control was guaranteed to score lowest.**

Following that rule would have deleted healthy controls for being healthy — in a study whose
entire point is comparing patients against healthy muscle, and with only 6 controls to spare.

**The fix, in two parts:**

1. **Compare like with like.** Score each sample only against others in its own group. The
   controls flip from lowest (0.9529) to effectively highest (**0.9839**).
2. **Measure the chips directly** using the two standard Affymetrix quality measures, which
   look *inside* a chip rather than comparing chips to each other, so biology cannot leak in:
   - **RLE** — does this chip read systematically high or low? Healthy ≈ 0, flag beyond ±0.10.
     **Our 34 chips: −0.047 to +0.057.** All clean.
   - **NUSE** — do the ~11 spots measuring each gene agree with each other? Typical = 1.00,
     flag above 1.05. **Our 34 chips: 0.981 to 1.039.** All clean.

**Decision: keep all 34** (confirmed by Aziz). Two samples remain *mildly* flagged on the
within-group measure — GSM2934803 and GSM2934813 — and are kept, because the measures that
actually test chip quality are clean and because **step 13 tests those exact samples directly**
by re-running the analysis 17 times, leaving one sample out each time.

### Step 3 — From chip spots to genes

Chip spots are not genes; several spots can measure the same gene, and some match nothing.

54,675 spots → 44,647 have a known gene → 1,524 ambiguous ones (matching several genes) dropped
→ 43,123 kept → collapsed to **20,823 genes**, keeping the brightest-reading spot per gene.

### Step 4 — Which genes differ between groups

Used **limma**, the standard tool for this. It makes three comparisons: BMD vs control, DMD vs
control, and BMD vs DMD directly.

> **Why limma rather than a plain t-test.** With 6 controls, a gene can look dramatic purely
> because its handful of measurements happened to land close together. limma borrows
> information across all 20,823 genes to estimate how much noise is normal, and corrects such
> genes down. Here it added the equivalent of **4.2 extra samples'** worth of stability.

**The sanity check — do the known scar genes behave?**

| Gene | BMD vs control | DMD vs control | BMD vs DMD |
|---|---|---|---|
| COL1A1 | **+2.52** (5.7×) | **+4.06** (16×) | −1.54 |
| COL1A2 | +2.04 (4.1×) | +2.94 (7.7×) | −0.90 |
| COL3A1 | +2.31 (5.0×) | +3.26 (9.6×) | −0.95 |
| COL5A1 | +1.92 (3.8×) | +2.83 (7.1×) | −0.91 |
| LOX | +2.61 (6.1×) | +3.57 (11.9×) | −0.97 |

*(Numbers are log2 fold changes; +1 means doubled, +2 means 4×. All adjusted p < 0.01.)*

Every scar gene is up in **both** diseases, and **less** up in BMD than in DMD. Exactly the
textbook picture: both diseases scar, BMD scars less. The pipeline is working.

These values also sit comfortably beside the published Eisenberg 2007 figures (~6× for COL3A1,
~5× for COL1A2 in DMD) — same direction, same order of magnitude.

### Step 5 — Building the microRNA target lists

2,064 human microRNA families → keep only those present in muscle (253) → keep only those with
at least 15 target genes on our array → **128 families** to test.

miR-29 carries **1,158** target genes.

### Step 6 — The smoke-alarm test ✅ **PASSED**

**The logic:** before trusting this method on BMD, where nobody knows the answer, check it
reproduces the *known* answer in DMD. miR-29 is known to fall in DMD. If our method cannot see
that, the method is broken and the paper becomes a methods paper instead.

**How the test works.** Line up all 20,823 genes left-to-right, most-increased-in-DMD first.
Walk along the line. Every time you step on a miR-29 target, the score goes **up**; every time
you step on a non-target, it drifts **down**, with step sizes set so that randomly scattered
targets would finish at zero.

> **Worked example with 10 genes, 3 of them targets (T).**
>
> Order `T T T n n n n n n n` → each T adds +1/3, each n subtracts 1/7:
> 0.33 → 0.67 → **1.00** → 0.86 → 0.71 → … → 0. **Peak 1.00, reached early.**
>
> Order `n n n n n n n T T T` → −0.14 → −0.29 → … → **−1.00**, then back to 0.
> Same size, **negative**, reached late.
>
> Order `n T n n T n n T n n` (scattered) → wanders near 0. **No footprint.**
>
> That peak, rescaled so sets of different sizes can be compared, is the **NES**.

**Result:**

| | |
|---|---|
| miR-29, DMD vs control | **NES = +2.22** |
| adjusted p-value | **5.4 × 10⁻³⁵** |
| rank among 128 families | **3rd** |
| Verdict | **GO** |

The plain-numbers version of the same claim, for anyone who distrusts rank tests: miR-29's
1,158 targets rose by **+0.304** on average, versus **+0.036** for the other 19,665 genes.
The gap — **+0.269**, about a 20.5% rise — is the "ruler" we use later to measure what size
of footprint we can detect.

NES > 0 means targets are up, which means the brake is weaker, which means **miR-29 is
predicted less active in DMD**. That is the published direction. **The method recovers known
biology.**

> ### ⚠️ The caveat that must travel with this result
>
> **103 of 128 families reached significance, and 99 of them pointed the same way.**
>
> This is normal and expected: microRNA target lists overlap heavily, and in DMD the whole
> transcriptome shifts, so almost every target list drifts upward together. A test that
> compares a gene set against *the rest of the genome* will light up for nearly everything.
>
> **So the smoke test proves the method works. It does NOT prove miR-29 is special.**
> Specificity has to be earned elsewhere — and it was. See Part 4.

### Step 7 — All 128 families, all three comparisons

| Comparison | Families significant |
|---|---|
| BMD vs control | 104 of 128 |
| DMD vs control | 103 of 128 |
| **BMD vs DMD (direct)** | **7 of 128** |

**miR-29 across the three:**

| Comparison | NES | adjusted p | rank |
|---|---|---|---|
| BMD vs control | +2.17 | 1.4 × 10⁻³¹ | 6 / 128 |
| DMD vs control | +2.22 | 1.5 × 10⁻³⁴ | 3 / 128 |
| **BMD vs DMD** | **−1.59** | **4.7 × 10⁻¹⁰** | **1 / 128** |

**This is the headline.** See Part 4.

---

## Part 4 — What we think we have found

### The result

When BMD and DMD are each compared against the shared healthy controls, almost everything looks
significant and miR-29 is merely one of a crowd.

**But when BMD and DMD are compared directly against each other** — a comparison in which the
controls play no part at all, so no quirk of those 6 healthy people can influence anything —
the crowd vanishes. Only **7 of 128** families remain, and **miR-29 is first, by roughly five
orders of magnitude** (4.7 × 10⁻¹⁰ versus 1.9 × 10⁻⁴ for the runner-up).

Its sign is negative, which by the direction rule means miR-29's targets are **lower in BMD
than in DMD** → the brake is **more active in BMD** → **less weakened in BMD than in DMD**.

**In plain words: BMD shows the same miR-29/collagen disturbance as DMD, but milder** — which is
exactly what a milder form of the same disease should show, and exactly what exon-skipping
therapy is trying to achieve.

### Three independent specificity signals, none of them a p-value

**1. The leading-edge genes are collagen.** The "leading edge" is the subset of genes actually
driving a family's score. miR-29's is **6.5% fibrosis-set genes — the highest of any top
family** (runner-up 4.1%; most sit at 1–3%). The list reads as a collagen roll-call:

> COL1A1, COL1A2, COL3A1, COL4A1, COL4A2, COL4A4, COL5A1, COL5A2, COL5A3, COL6A2, COL6A3,
> COL6A6, COL15A1, COL16A1, COL21A1, COL27A1, ELN, FBN1, LAMA2, LAMC1, NID1, SPARC, MFAP2, SKI…

This is content-based evidence. It does not depend on any statistical threshold.

**2. It is not merely scar tissue** (step 10a). A fair objection: diseased muscle is full of
scar, so maybe we are just rediscovering collagen. We deleted **all 275 core scar-matrix genes**
from both the ranked list and the target sets, and re-ran. miR-29 retained **96.7%** of its
score (+2.224 → +2.151, p = 1.4 × 10⁻²⁹).

**3. It is not merely fibroblast content** (step 12). Second objection: diseased biopsies simply
contain more scar-producing cells (fibroblasts), so we might be measuring cell composition. We
built a per-sample fibroblast score from 316 fibroblast markers that are *not* themselves scar
genes, and put it into the model.

> Fibroblast score: control **−0.79**, BMD **+0.04**, DMD **+0.25** (p = 1.2 × 10⁻⁸).
> BMD sits between healthy and DMD — the expected ordering.

This adjustment cut the significant families from **103 to 51** — and miR-29 **improved its
relative standing**, keeping **79%** of its score where other top families kept 55–60%, and
afterwards leading by eleven orders of magnitude (1.7 × 10⁻¹⁶ vs 3.1 × 10⁻⁵).

> **Caveat required by the spec:** scar load genuinely does differ between groups, so this
> covariate also absorbs some real disease signal. Read it as *"does the footprint survive?"* —
> never as the headline estimate.

---

## Part 5 — What we deliberately do NOT yet claim

### The trap we have not escaped yet: shared controls

BMD-vs-control and DMD-vs-control use **the same 6 control samples**.

> **Worked example of why that is dangerous.** Take two genes with *no real difference* between
> any group, where by chance the 6 controls read slightly oddly:
>
> | Gene | Control avg | BMD avg | DMD avg | BMD−Ctrl | DMD−Ctrl |
> |---|---|---|---|---|---|
> | 1 | 9.70 | 10.05 | 9.98 | **+0.35** | **+0.28** |
> | 2 | 10.30 | 9.96 | 10.02 | **−0.34** | **−0.28** |
>
> Both columns move together — not because BMD and DMD share biology, but because they share
> the same flawed baseline. **With our exact group sizes, pure noise alone produces a
> correlation of about 0.7 with a slope just under 1** — which is precisely the "BMD is a
> weaker copy of DMD" picture, manufactured from nothing.

Our measured agreement between BMD and DMD footprints is **r = 0.9123**, slope **0.9769**.

**We cannot interpret those numbers yet**, and they must never be judged against zero. They are
only meaningful against a **shuffled-label baseline**: scramble the group labels at random
1,000 times, keeping 6/11/17, rerun the entire pipeline each time, and see what correlation
appears when nothing real is left.

**That is step 8.** It was started in this session but **did not finish** — it hit a
performance bug (see below) and was stopped with no result. Until it is re-run and lands, the
claim "BMD is a weaker copy of DMD" is **not supported** and must not be written.

> **Why step 8 did not finish — worth knowing, because it is a good engineering story.**
> The job was projected at 98 minutes. After 99 minutes it had done about 19% of the work.
> The cause: our shuffle loop runs 7 parallel workers, and the `fgsea` tool *inside* each
> worker defaults to grabbing every core on the machine. So 7 workers each tried to use all 8
> cores — up to 56 threads fighting over 8. The machine was spending **52% of its effort on
> switching between tasks** instead of computing, and each worker sat idle ~80% of the time.
> Fixed by telling the inner tool to stay single-threaded. Logged in `docs/methods_log.md`.

### Everything else still outstanding

| Step | Question it answers | Why it matters |
|---|---|---|
| **8** *(must be re-run)* | Does r = 0.91 beat the shared-control baseline? | Gates the "weaker copy" claim entirely |
| **9** | Does miR-29's BMD/DMD difference survive shuffling patient labels? | Confirms the headline |
| **11** | Is this a long-3′UTR artefact? | Genes with long tails appear on many target lists *and* behave alike in disease — a plain random comparison is biased in our favour |
| **10b** | How big must a footprint be for us to detect it at n = 11? | Without this, "we found nothing" is uninterpretable — it could mean "nothing there" or "too small to see" |
| **13** | Does any single sample carry the result? | Also settles the two QC-flagged samples |
| **14** | Is BMD's weaker signal just a sample-size effect? | BMD has 11 patients, DMD has 17 |
| **15** | Pig cross-species comparison | **BLOCKED** — required source table unavailable |

---

## Part 6 — Honest limitations

1. **We never measured a microRNA.** Everything here is inferred from how target genes behave.
   Correct language throughout: *predicted*, *candidate*, *may contribute*. Not *causes*,
   not *biomarker*.
2. **n = 11 BMD.** Every result is hypothesis-generating. Step 10b will tell us exactly how
   small a real effect we would have missed.
3. **Control ages are unknown.** Not recorded in GEO at all. Unresolvable.
4. **BMD patients are ~10 years older than DMD patients**, so the direct comparison — the one
   carrying our headline — confounds disease with age. Cannot be corrected without also
   removing the disease signal.
5. **Target predictions are predictions.** TargetScan infers targeting from sequence matching;
   it does not prove any particular gene is really braked by miR-29 in human muscle.
6. **One dataset, one muscle, one chip platform.** No independent replication (the planned
   replication dataset was cut earlier for other reasons).

---

## Part 7 — Reference: where everything lives

**Before running anything:** `export PATH="/opt/homebrew/bin:$PATH"` — R is not on the default
path.

| Path | Contents |
|---|---|
| `scripts/00–14` | one script per step, each runs standalone |
| `docs/methods_log.md` | every parameter, default, version, deviation, cut — with reasons |
| `docs/results_log.md` | every number you might quote, with the script that produced it |
| `docs/figure_notes.md` | per figure: axes, a worked example point, what patterns mean |
| `docs/fact_check.md` | confirmed vs `[UNVERIFIED]` claims |
| `docs/handoffs.md` | the auto-built inputs, for Hamza to diff against his |
| `docs/SESSION_01_HANDOFF.md` | how to resume, in order |
| `results/tables/`, `results/figures/` | outputs |

**Figures produced so far:** QC boxplot · QC within-group correlation · QC NUSE ·
**miR-29 enrichment plot (the ranked walk)** · NES heatmap across the three comparisons ·
fibroblast score by group.
