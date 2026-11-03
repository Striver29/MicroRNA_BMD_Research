# BMD microRNA–Fibrosis Project — Complete Analysis Write-up

**Dataset:** GSE109178 · **Sessions:** 1–2 Nov 2026 · **Computation freeze:** Friday 6 Nov
**Status:** All analyses complete except the detection-limit curve (step 10b, handed to Hamza).

> **How to read this page.** It assumes no biology and no statistics. Every term is defined
> where it first appears, and every statistical claim comes with a worked example using real
> numbers from our own data. If something isn't explained, that's a bug — tell me.

---

# PART 0 — The result in one minute

We asked whether a molecule called **miR-29** behaves differently in two related muscle
diseases, BMD and DMD.

We found a pattern that looked spectacular. Then we tested it the correct way, and it
**dissolved**. Our honest conclusion is:

> With 6 healthy people, 11 BMD patients and 17 DMD patients, we **cannot tell** whether the
> microRNA differences we see are real or are coincidences of which people happened to be in
> which group.

This is a **null result**. It is not a failed project. Three things make it worth publishing:

1. **We showed the method points the right way.** miR-29's target genes behave exactly as the
   published biology predicts, and that survives every bias check we threw at it.
2. **We found that the standard software is dramatically overconfident.** The tool everyone
   uses for this (`fgsea`) gave a p-value of 10⁻³⁴ where the honest answer is 0.052. That is an
   overstatement by **34 orders of magnitude**, and it is a warning other researchers need.
3. **We caught a trap in the public data archive.** A dataset that looked like independent new
   samples turned out to be *the same physical experiments re-uploaded under a different
   name*. Merging it would have faked a significant result.

---

# PART 1 — The biology, from scratch

## 1.1 What the diseases are

Muscle cells are long fibres that contract. Inside each fibre is a protein called **dystrophin**
that acts as a shock absorber — it stops the fibre tearing itself apart when it contracts.

**Duchenne muscular dystrophy (DMD)** happens when the dystrophin gene is broken so badly that
no working protein is made. Every contraction damages the fibre. Fibres die and are replaced by
**collagen scar tissue** — a process called **fibrosis**. Muscle turns, slowly, into scar.
Boys with DMD typically stop walking in their early teens.

**Becker muscular dystrophy (BMD)** is *the same gene*, broken less badly. The body still makes
a shortened dystrophin that partly works. Same disease process, much slower. Some people with
BMD walk into their sixties — our dataset includes a BMD patient aged **60.8**.

So BMD and DMD are not two diseases. They are **one disease at two severities**.

## 1.2 Why the comparison matters

The leading DMD therapy is called **exon skipping**. It works by editing the genetic message so
that a DMD-type mutation produces a BMD-type protein. It deliberately **converts DMD into
something like BMD**.

So "what does BMD actually look like molecularly?" is the same question as "what is this
therapy aiming at?" And remarkably, **nobody had tested whether the microRNA–fibrosis link
known in DMD also runs in human BMD.** That gap is our project.

## 1.3 What a microRNA is

A gene is an instruction. The cell copies it into a working message (**mRNA**), then reads that
message to build a protein.

A **microRNA** is a very short piece of RNA — about 22 letters long — that acts as a **brake**.
It finds messages and shuts them down, so less protein gets made.

How does it know which messages to grab? Through a 7-letter stretch called the **seed** (letters
2 through 8 of the microRNA). The seed pairs with a matching 7-letter site in the tail end of a
message — a region called the **3′UTR**. Match found → message destroyed or blocked.

Because one microRNA has one seed, it brakes **every** message carrying that matching site —
often hundreds or thousands of genes simultaneously.

> **Worked example.** miR-29's seed is `AGCACCA`. In our data, **1,158 genes** present on our
> chip carry a matching site. miR-29 is a single brake pressing on all 1,158 at once. A large
> number of them are **collagen genes** — the scar proteins.

**The known fact we built on:** in DMD, miR-29 levels go **down**. The brake on collagen comes
off, collagen accumulates, scarring accelerates. This is established in the literature.

## 1.4 The problem: we can't see microRNAs

**Our dataset contains no microRNA measurements at all.** It measures genes only.

So we cannot look at miR-29 directly. We look for its **footprint**.

> **The analogy.** You cannot see wind. But you can see which way the trees lean. If miR-29's
> brake weakens, then all 1,158 of its target genes should drift **upward together** — a
> coordinated nudge, too consistent to be chance. We detect the brake by the shadow it casts.

## 1.5 The direction rule (easy to get backwards)

Our main test outputs a number called **NES** (explained fully in Part 3).

| What we see | What it means |
|---|---|
| **NES > 0** — targets bunched at the **UP** end | brake **weaker** → microRNA **LESS** active |
| **NES < 0** — targets bunched at the **DOWN** end | brake **stronger** → microRNA **MORE** active |

It feels inverted because a brake and the thing it brakes always move in **opposite**
directions. Brake off → targets up. So a microRNA that genuinely *increases* gives **NES < 0**.

---

# PART 2 — Building the machine

## 2.1 Nothing was installed

The entire analysis plan is built on R and Bioconductor. **Neither was on the machine.**
Installed R 4.6.2 and Bioconductor 3.23, every package compiled from source — about an hour,
run in the background while other work continued. One package (`GEOquery`) failed to build; it
turned out we didn't need it.

## 2.2 The data archive blocked us

The data live at NCBI's **GEO** archive, accession **GSE109178**. The standard download route is
NCBI's FTP server, which returns **403 Forbidden** from this machine. Worked around using GEO's
two official web endpoints, which serve identical content. Logged as a deviation.

## 2.3 Six input files that didn't exist yet, built from scratch

| File | What it is | Source |
|---|---|---|
| `mirna_families.tsv` | 2,064 microRNA families + their 7-letter seeds | TargetScanHuman 8.0 |
| `mirna_targets.tsv` | 115,308 microRNA→gene predictions | TargetScanHuman 8.0 |
| `utr_lengths.tsv` | 3′UTR length for 19,429 genes | TargetScan UTR sequences |
| `fibrosis_sets.tsv` | 275 scar-matrix genes + 54 TGF-β genes | MSigDB (Naba 2012) |
| `fibroblast_markers.tsv` | 372 markers of the scar-producing cell | MSigDB C8 (Rubenstein 2020) |
| `sample_sheet.csv` | who is in which group | parsed from GEO |

**A judgement call:** we had to filter to microRNAs actually present in muscle, choosing a
published source ourselves. We used the **miRNA Tissue Atlas 2025**, taking its *limb skeletal
muscle* samples — the same kind of muscle as our thigh biopsies — with a cutoff of mean
RPMM ≥ 10. **That cutoff is ours, not the spec's**, and is logged as a parameter to test if
challenged.

> **A check that the atlas is trustworthy.** The most abundant microRNAs it reports in skeletal
> muscle are miR-486, miR-10b, **miR-133a, miR-1, miR-206**. The last three are the textbook
> muscle-specific "myomiRs" — famous for existing almost nowhere else in the body. Their
> appearing at the top is strong evidence the data are behaving.

miR-29 in muscle: miR-29a-3p at **1,895 RPMM**, miR-29b-3p 650, miR-29c-3p 439. Clearly present.

---

# PART 3 — The analysis, step by step

## Step 1 — Who is in the dataset

49 muscle biopsies, all from the **vastus lateralis** (a thigh muscle), all on the same chip.

| Group | n | Keep? |
|---|---|---|
| Healthy control | **6** | ✅ |
| BMD | **11** | ✅ |
| DMD | **17** | ✅ |
| Limb-girdle dystrophy | 15 | ❌ different disease |

Matched the specification exactly. **Two confounds we found and cannot fix:**

1. **Control ages are simply absent from GEO.** All 6 records say `age: NA`. We were asked to
   record them; they do not exist. Unresolvable — goes in Limitations.
2. **BMD patients are ~a decade older than DMD patients.** Median **12.2 years** (range
   2.5–60.8) versus **3.5 years** (range 0.9–8.0). DMD is diagnosed and biopsied young; BMD is
   milder and caught later. So our BMD-vs-DMD comparison **mixes disease severity with age**,
   and cannot be corrected without also deleting the disease signal.

## Step 2 — From chip images to numbers, and a QC trap

Each chip has ~54,675 spots. Converting raw brightness into comparable numbers uses a standard
procedure called **RMA**: remove background glow, force all 34 chips onto a common scale,
combine the ~11 spots measuring each gene into one value. Result: a **54,675 × 34** table.

### The trap (a good story for the professor)

The plan said: flag any sample whose correlation with all the others is unusually low. A broken
chip agrees with nobody. Applied literally, this flagged **three of our six controls** — and in
fact all six controls were the six lowest-scoring samples.

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
> lands inside the healthy-vs-diseased comparisons. **With 28 diseased chips and only 6 healthy
> ones, every control was guaranteed to score lowest.**

Following that rule would have deleted healthy controls *for being healthy*, in a study whose
entire purpose is comparing patients against healthy muscle — with only 6 controls to spare.

**The fix, two parts:**

1. **Compare like with like.** Score each sample only against others in its own group. Controls
   flip from lowest (0.9529) to effectively highest (**0.9839**).
2. **Measure the chips directly** with the two standard Affymetrix quality measures, which look
   *inside* a chip rather than comparing chips to each other:
   - **RLE** — does this chip read systematically high or low? Healthy ≈ 0, flag beyond ±0.10.
     **Ours: −0.047 to +0.057.** All clean.
   - **NUSE** — do the ~11 spots measuring each gene agree? Typical = 1.00, flag above 1.05.
     **Ours: 0.981 to 1.039.** All clean.

**Decision: keep all 34.** Later confirmed by step 13, which showed no single sample matters.

## Step 3 — From chip spots to genes

54,675 spots → 44,647 have a known gene → 1,524 ambiguous dropped → 43,123 kept → collapsed to
**20,823 genes**, keeping the brightest-reading spot per gene.

## Step 4 — Which genes differ, and does the biology make sense?

We used **limma**, the standard tool.

> **Why limma and not a plain t-test.** With 6 controls, a gene can look dramatic purely because
> its handful of measurements happened to land close together. limma borrows information across
> all 20,823 genes to estimate how much noise is normal, and corrects such genes down. Here it
> added the equivalent of **4.2 extra samples'** worth of stability.

**The sanity check — do the known scar genes behave?**

| Gene | BMD vs control | DMD vs control | BMD vs DMD |
|---|---|---|---|
| COL1A1 | **+2.52** (5.7×) | **+4.06** (16×) | −1.54 |
| COL1A2 | +2.04 (4.1×) | +2.94 (7.7×) | −0.90 |
| COL3A1 | +2.31 (5.0×) | +3.26 (9.6×) | −0.95 |
| COL5A1 | +1.92 (3.8×) | +2.83 (7.1×) | −0.91 |
| LOX | +2.61 (6.1×) | +3.57 (11.9×) | −0.97 |

*(log2 fold changes: +1 = doubled, +2 = 4×. All adjusted p < 0.01.)*

Every scar gene is up in **both** diseases, and **less** up in BMD than DMD. Exactly the
textbook picture: both diseases scar, BMD scars less. **The pipeline works.** These values also
sit comfortably beside published figures (Eisenberg 2007: ~6× COL3A1, ~5× COL1A2 in DMD).

## Step 5 — Building the microRNA target lists

2,064 human families → keep those present in muscle (253) → keep those with ≥15 targets on our
chip → **128 families tested**. miR-29 carries **1,158** target genes.

**Note: miR-29 is not our only interest.** It is our *positive control* — the one case where the
answer is already known, so it tells us whether the method works. All 128 are tested every time.

## Step 6 — The smoke-alarm test

**The logic:** before trusting this method on BMD, where nobody knows the answer, check it
reproduces the *known* answer in DMD.

**How it works.** Line up all 20,823 genes left to right, most-increased-in-DMD first. Walk
along. Every time you step on a miR-29 target the score goes **up**; every non-target drifts it
**down**, with step sizes set so randomly scattered targets would finish at zero.

> **Worked example with 10 genes, 3 of them targets (T).**
>
> Order `T T T n n n n n n n` → each T adds +1/3, each n subtracts 1/7:
> 0.33 → 0.67 → **1.00** → 0.86 → 0.71 → … → 0. **Peak 1.00, reached early.**
>
> Order `n n n n n n n T T T` → −0.14 → −0.29 → … → **−1.00**, then back to 0.
> Same size, **negative**, reached late.
>
> Order `n T n n T n n T n n` → wanders near 0. **No footprint.**
>
> That peak, rescaled so different-sized sets are comparable, is the **NES**.

**Result:** miR-29, DMD vs control → **NES = +2.22**, fgsea padj = **5.4 × 10⁻³⁵**, rank **3 of
128**. On its own terms: a spectacular pass.

The plain-numbers version: miR-29's 1,158 targets rose **+0.304** on average versus **+0.036**
for the other 19,665 genes. The gap — **+0.269 log2, a 20.5% rise** — is the **"ruler"** used
later to measure what size of footprint we can detect.

> ### ⚠ The caveat that must travel with this
> **103 of 128 families reached significance, and 99 pointed the same way.** Target lists
> overlap heavily, and in DMD the whole transcriptome shifts, so almost every list drifts up
> together. **The smoke test shows the method works. It does NOT show miR-29 is special.**

## Step 7 — All 128 families, three comparisons

| Comparison | families fgsea called significant |
|---|---|
| BMD vs control | 104 / 128 |
| DMD vs control | 103 / 128 |
| **BMD vs DMD (direct)** | **7 / 128** |

miR-29: +2.17 (BMD vs ctrl), +2.22 (DMD vs ctrl), **−1.59 (BMD vs DMD, rank 1 of 128)**.

**At this point we thought we had our headline.** We were wrong, and steps 8–9 are why.

---

# PART 4 — The three bias checks that miR-29 passed

## Step 10a — Is it just scar tissue?

A fair objection: diseased muscle is full of scar, so maybe we're just rediscovering collagen.
We deleted **all 275 core scar-matrix genes** from both the ranked list and the target sets, and
re-ran.

**miR-29 retained 96.7% of its score** (+2.224 → +2.151). Not just scar.

## Step 12 — Is it just more scar-producing cells?

Second objection: diseased biopsies contain more **fibroblasts** (the cells that make scar), so
maybe we're measuring cell composition. We built a per-sample fibroblast score from 316
fibroblast markers that are *not* themselves scar genes, and added it to the model.

> Fibroblast score: control **−0.79**, BMD **+0.04**, DMD **+0.25** (p = 1.2 × 10⁻⁸).
> BMD sits between healthy and DMD — the expected ordering.

This cut significant families from 103 to 51 — and miR-29 **improved its relative standing**,
keeping **79%** of its score where rivals kept 55–60%.

## Step 11 — The long-tail trap

Genes with long 3′UTRs have more room for microRNA sites, so they appear on **many** target
lists — and they also behave alike in disease. A plain random comparison is therefore biased in
our favour. We drew 1,000 random gene sets per family **matched on 3′UTR length decile by
decile**.

**A UTR-matched random set already scores NES ≈ +1.6. Not zero.**

| | real NES | matched null | z | p |
|---|---|---|---|---|
| miR-29, DMD vs ctrl | 2.222 | **1.592** ± 0.096 | **6.56** | 0.000999 |
| miR-29, BMD vs ctrl | 2.172 | **1.669** ± 0.093 | **5.39** | 0.000999 |

**About 72% of miR-29's apparent enrichment is explained by 3′UTR length alone.** The footprint
is still real — the excess is 6.5 standard deviations above a fair baseline, and not one of
1,000 matched sets reached it — but **NES must never be quoted as though zero were the
baseline.**

## Step 13 — Does one patient carry everything?

We dropped each of the 11 BMD and 6 control samples in turn — 17 complete re-runs.

| | full data | range across 17 runs | significant in |
|---|---|---|---|
| miR-29, BMD vs control | +2.169 | +2.120 to +2.186 | **17/17** |
| miR-29, BMD vs DMD | −1.590 | −1.674 to −1.500 | **17/17** |

miR-29 stayed **rank 1 of 128 in every single BMD-vs-DMD run**. The two samples flagged back in
step 2 proved harmless, and the disputed-severity sample GSM2934804 changed nothing. **The
QC decision to keep all 34 is now evidence-backed, not a judgement call.**

---

# PART 5 — Where it all came apart

## 5.1 The idea: shuffle the labels

You can always find patterns — stare at clouds and you'll see a face. So the real question is:
**if there were nothing here at all, how often would I see a pattern this strong anyway?**

The standard way to answer that is to **scramble the labels and try again.**

> **Worked example with heights.** Six people:
>
> - **Group A:** 170, 172, 174 cm → average **172**
> - **Group B:** 176, 178, 180 cm → average **178**
>
> Group B is 6 cm taller. Real?
>
> Write all 6 heights on cards. Shuffle. Deal 3 into a fake "A" and 3 into a fake "B". Measure
> the gap. Do it 500 times.
>
> - Fake groups **almost never** reach 6 cm → your result is real.
> - Fake groups reach 6 cm **all the time** → you have nothing. Any 3-vs-3 split produces *some*
>   gap; yours isn't special.

That is the whole idea.

## 5.2 fgsea was shuffling the wrong thing

`fgsea` gave us p = 0.00000000000000000000000000000000015 — "this essentially never happens by
luck."

**But fgsea keeps the people fixed and shuffles the *genes*.** That answers "is this particular
list of 1,158 genes special compared to other gene lists?" — not our question, which is "is the
difference between these *people* real?"

So we did it properly: **shuffled the people.** Randomly reassigned who counts as healthy, BMD
and DMD, keeping the group sizes at 6/11/17, and re-ran the entire pipeline from scratch. 500
times. (This is steps 8 and 9, and it is where most of the computing time went.)

## 5.3 Result: everything dissolves

**Step 9 — the direct BMD vs DMD comparison:**

| Family | fgsea padj | honest (patient-shuffle) p |
|---|---|---|
| miR-24-3p | 1.9e−04 | 0.214 |
| **miR-29-3p** | **4.7e−10** | **0.239** |
| miR-133a/133b | 1.5e−03 | 0.449 |
| miR-30-5p | 1.9e−04 | 0.533 |
| miR-1-3p/206 | 4.1e−04 | 0.533 |

**Zero of 128 families survive.** miR-29's real NES of −1.59 sits **1.19 standard deviations**
from a null centred at 0.05 with SD **1.386**. Reshuffling which patients are "BMD" routinely
produces differences this large.

**Step 8 — the vs-control comparisons:**

| Family (DMD vs control) | fgsea padj | honest p | overstated by |
|---|---|---|---|
| miR-30-5p | 2.3e−36 | 0.070 | 10³⁴ |
| miR-1-3p/206 | 2.7e−36 | 0.080 | 10³⁴ |
| **miR-29-3p** | **1.5e−34** | **0.052** | **10³²** |
| miR-23-3p | 2.0e−28 | 0.058 | 10²⁶ |

**miR-29's famous 5.4 × 10⁻³⁵ becomes p = 0.052.** It misses significance. BMD vs control is
p = 0.084.

> **Honesty note on the ceiling.** With 500 shuffles the smallest possible p is 1/501 = 0.002.
> After correcting for 128 families, the smallest achievable adjusted value is 0.255 — so
> "0 families significant after correction" is partly a **resolution limit**, not purely a
> finding. These should be judged on raw p. But miR-29 sits at 0.052 / 0.084 / 0.239, so more
> shuffles would sharpen the estimate, not rescue it.

## 5.4 The shared-control trap

We had also noticed BMD and DMD footprints correlate at **0.91** out of 1.0 — seemingly strong
evidence they're the same process at two severities.

**But both were compared against the *same* 6 healthy people.** If those 6 are unusual in any
way — and with only 6, they will be — then **both** comparisons get knocked in the same
direction, making BMD and DMD look alike automatically.

| | measured | shuffled-label null | p |
|---|---|---|---|
| **r** | **0.9123** | mean 0.510, median 0.553, **95th pct 0.933** | **0.064** |
| **slope** | 0.9769 | mean 0.516, 95th pct 1.094 | **0.114** |

**The null's 95th percentile (0.933) is higher than our observed 0.912.** Shuffled labels —
where nothing real remains — routinely manufacture the "BMD is a weaker copy of DMD" picture.

## 5.5 Step 14 — and it isn't a sample-size artefact either

Could BMD merely *look* weaker because it has 11 patients and DMD has 17? We randomly cut DMD
down to 11 and re-ran, 100 times.

| | miR-29 NES vs control |
|---|---|
| full DMD (n = 17) | 2.224 |
| **DMD cut to 11** | **median 2.167** (2.074–2.258) |
| **actual BMD (n = 11)** | **2.169** |

**53.0% of the cut-down DMD draws fall at or below BMD's value** — a coin flip. BMD sits dead
centre of what DMD itself produces at the same sample size.

So: against controls, **BMD and DMD are indistinguishable**. Head to head, **also
indistinguishable**. Two independent routes to the same answer.

---

# PART 6 — The duplicate-dataset trap

When we realised the problem was sample size, the obvious move was to find more data. We
searched GEO properly. One dataset looked perfect: **GSE38417** — same chip, 6 controls and 16
DMD, apparently free extra power.

**It is not a new dataset. It is the same physical experiments re-uploaded under a different
accession number.**

We proved it at the raw-data level, comparing one control from each:

| Check | Result |
|---|---|
| probe cells compared | 1,354,896 |
| intensities identical | **TRUE** |
| correlation | **1.00000000** |
| cells differing | **0** |

Same lab, same internal array IDs (`MIDO-CTRL_Pos-1UP` … appears in **both**), only the
accession and file format differ.

**Had we merged it**, we'd have counted the same 6 controls twice. Standard errors would have
shrunk by roughly √2 on pure fiction, and miR-29's p = 0.052 would likely have crossed into
"significant" — a result manufactured entirely by duplication. **This is the single worst error
available in this project, because it would have looked like success.**

`scripts/check_duplicate_arrays.R` now exists as a mandatory gate before merging any accession.

**The other candidates:**
- **GSE13608** — 6 genuine controls, but biopsy site "unknown" (ours are all vastus lateralis),
  ages 26–85 versus DMD children aged 0.9–8, one sample is a **pool** of several people. Buys
  only 6→12 controls (~18% precision) for two new confounds. **Rejected.**
- **GSE291383** — 6 BMD, 12 DMD, **no controls**, RNA-seq. The only other human BMD muscle
  dataset in existence.

**Conclusion: 11 BMD patients is the binding constraint, and no public dataset fixes it.**

---

# PART 7 — What this means, honestly

## What we can claim

1. **The method recovers known biology directionally.** miR-29's targets rise in dystrophic
   muscle, the genes driving it are collagens (COL1A1, COL1A2, COL3A1, COL4A1/2/4, COL5A1/2/3,
   COL6A2/3/6, ELN, FBN1, LAMA2, SPARC…), and this survives deleting all scar genes, adjusting
   for fibroblast content, a 3′UTR-matched null, and dropping any single patient.
2. **miR-29 ranks 3rd of 128 even on the honest test** (p = 0.052), inside a cluster of families
   all hovering at the edge. That is the signature of **underpowered**, not **absent**.
3. **This is the first test of the microRNA–fibrosis link in human BMD.** The answer is "we
   cannot resolve it", which is itself new information.

## What we must NOT claim

- ❌ miR-29 is significantly altered in BMD or DMD *(honest p = 0.084 and 0.052)*
- ❌ BMD's footprint is a weaker copy of DMD's *(r p = 0.064; step 14 p ≈ 0.53)*
- ❌ BMD and DMD differ in any microRNA family *(0 of 128 survive; miR-29 p = 0.239)*
- ❌ Anything causal. We never measured a microRNA — only inferred it from target behaviour.

**Required language:** *predicted*, *candidate*, *may contribute*. Never *causes*, never
*biomarker*.

## The limitations list

1. **We never measured a microRNA.** Everything is inferred from target-gene behaviour.
2. **n = 11 BMD.** The binding constraint, and unfixable from public data.
3. **Control ages unknown** — absent from GEO entirely.
4. **BMD patients ~10 years older than DMD patients**, so the direct comparison confounds
   disease with age.
5. **Target predictions are predictions.** TargetScan infers targeting from sequence matching;
   it does not prove miR-29 really brakes any given gene in human muscle.
6. **One dataset, one muscle, one chip platform.** No independent replication.

---

# PART 8 — The one analysis still running

**Step 10b, the detection limit** — handed to Hamza (needs ~15.5 CPU-hours; see
`docs/HANDOFF_hamza_step10b.md`).

**Why it is now the most important result in the paper.** A null result only means something if
you can say **how big an effect you would have caught**. Otherwise "we found nothing" is
indistinguishable from "our method is blind."

**How it works.** Plant a fake footprint of known size into the 11 BMD samples, then see how
often our own pipeline finds it. The **ruler** is the real DMD miR-29 footprint: **+0.2686 log2
units, a 20.5% rise**.

> **Worked example.** At f = 50% we add 0.5 × 0.2686 = **+0.134** to each of ~1,158 planted
> genes in the 11 BMD samples only — about a **10% rise**. Then re-run everything and ask
> whether we catch it. Repeat 100 times. If 62 of 100 runs flag it, that point is at 62%.
>
> **At f = 0 nothing is planted**, so anything caught there is a **false alarm** — that point
> should sit near 0.05.

The **detection limit** is where the curve crosses **80%**. That single number turns "we found
nothing" into "here is exactly how large an effect would have had to be, and why 11 patients
isn't enough" — which is a far stronger paper.

---

# PART 9 — Figures produced

| Figure | Shows |
|---|---|
| `02_qc_boxplot.png` | all 34 chips line up after normalisation |
| `02_qc_within_group_cor.png` | the QC trap resolved — controls flip from worst to best |
| `02_qc_nuse.png` | per-chip noise, all within tolerance |
| `06_mir29_enrichment_DMDvsCtrl.png` | **the ranked walk** — miR-29's targets crowded at the up end |
| `07_nes_heatmap.png` | all three comparisons for the top 25 families |
| `08_shuffle_r.png` / `08_shuffle_slope.png` | **the shared-control trap** — our value inside the null |
| `09_mir29_bmd_vs_dmd_shuffle.png` | **miR-29 BMD-vs-DMD sitting inside the shuffle distribution** |
| `12_fibroblast_score.png` | fibroblast content: control < BMD < DMD |
| `13_mir29_leave_one_out.png` | no single patient carries the result |
| `14_subsample_dmd.png` | BMD sits dead centre of DMD-at-n=11 |
| *(pending)* `10b_detection_curve.png` | the detection limit |

Every figure has axes, a worked example point, and "which pattern means what" documented in
`docs/figure_notes.md`.

---

# PART 10 — Where everything lives

**Before running anything:** `export PATH="/opt/homebrew/bin:$PATH"` — R is not on the default
path.

| Path | Contents |
|---|---|
| `scripts/00`–`14` | one script per step, each runs standalone |
| `scripts/check_duplicate_arrays.R` | **mandatory** before merging any GEO accession |
| `docs/methods_log.md` | every parameter, version, deviation, cut — with reasons |
| `docs/results_log.md` | every number you might quote, with the script that produced it |
| `docs/figure_notes.md` | per figure: axes, worked example, what patterns mean |
| `docs/fact_check.md` | confirmed vs `[UNVERIFIED]`, plus the duplicate-dataset finding |
| `docs/HANDOFF_hamza_step10b.md` | the outstanding analysis |
| `results/tables/`, `results/figures/` | all outputs |

## Lessons that cost us time (so they don't recur)

1. **`fgsea` nested inside parallel workers defaults to grabbing every core.** 7 workers × 8
   cores on an 8-core machine = 52% of CPU lost to context-switching; a 98-minute job projected
   out to 6.8 hours. Fix: `nproc = 1` inside parallel loops.
2. **`setorder(dt, -abs(x))` is invalid in data.table.** It crashed one step *after* two hours
   of computing, and was latent in another.
3. **Scripts must save expensive results before doing anything optional.** Step 9 died on its
   final line reading a file that didn't exist yet, nearly losing a 2-hour run.
4. **Never run two heavy R jobs at once** on this machine. They starve each other.
5. **Runtime projections from serial timing loops are unreliable.** Measure a small *parallel*
   batch instead.
