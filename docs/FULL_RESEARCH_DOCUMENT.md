# The BMD microRNA–Fibrosis Study: Complete Reference Document

**Dataset:** GEO accession GSE109178 · **Analysis sessions:** 1–6 November 2026
**Computation frozen:** Friday 6 November 2026 · **Status:** steps 0–14 complete, step 15 dropped
**Written for:** Aziz, as the single source to draft the Methods, Results, Discussion and
Limitations from. Every number here traces to a script in this repository that actually ran.

> ### How to read this document
>
> It assumes **no biology and no statistics**. Every term is defined the first time it is used,
> and every statistical claim is followed by a worked example with real numbers from our own
> data. It is long on purpose — it is a reference, not an essay. Read Parts 1–4 once to build
> the vocabulary, then use Part 5 as a per-step lookup while writing.
>
> **If a sentence here conflicts with an older note in `docs/`, this document wins.** It was
> written after the computation freeze, against the final numbers.

---

# PART 1 — The biology, from zero

## 1.1 Muscle, and what goes wrong in muscular dystrophy

Skeletal muscle is the meat that moves your skeleton. It is built from long cylindrical cells
called **muscle fibres**. Each fibre contracts by sliding internal protein filaments past one
another. To survive the mechanical stress of contracting, a fibre needs its flimsy outer
membrane anchored firmly to the scaffolding outside the cell. The protein that does the
anchoring is called **dystrophin**.

Dystrophin is encoded by the *DMD* gene on the X chromosome — the largest gene in the human
genome. If that gene is broken, fibres have no shock absorber. Every contraction tears them a
little. They die, regenerate, die again, and eventually the regeneration machinery is exhausted.

**This is muscular dystrophy: progressive muscle destruction caused by a faulty structural
protein.** It is genetic, present from birth, and gets worse with time.

## 1.2 Duchenne versus Becker — the same gene, two severities

Both diseases come from mutations in the same *DMD* gene. What separates them is the **reading
frame**.

DNA is read in groups of three letters, each group specifying one amino acid. Think of the gene
as a sentence made of three-letter words:

```
THE BIG RED CAT ATE THE RAT        <- the healthy gene
```

Delete a chunk whose length is *not* a multiple of three and everything after the deletion
shifts, turning into nonsense:

```
THE BIG RED ATE THE RAT            <- delete "CAT" (3 letters: frame preserved)
THE BIG RED CTA TET HER AT         <- delete "A" (1 letter: frame destroyed)
```

- **Duchenne muscular dystrophy (DMD)** = the frame is destroyed. The protein is truncated into
  uselessness. **Essentially no working dystrophin.** Symptoms start around age 2–5, wheelchair
  by about 12, historically fatal in early adulthood.
- **Becker muscular dystrophy (BMD)** = the frame is preserved. The protein is shortened but
  still partly functional. **Some working dystrophin.** Milder, slower, hugely variable — some
  patients are diagnosed in childhood, others stay active into their sixties.

> **This is the single most important fact for understanding why the project exists.** BMD is
> not a different disease from DMD. It is the *same* disease with a partly working protein. In
> our data the BMD patients range from 2.5 to 60.8 years old; the DMD patients are all children,
> 0.9 to 8.0 years old. That age gap is a real problem we will have to confess to (Part 8).

## 1.3 Fibrosis: muscle turning into scar

When a muscle fibre dies, the body patches the hole. The patch is not new muscle — it is
**connective tissue**, the same material as a scar: dense, tough, inert, non-contractile. The
main ingredient is **collagen**, a rope-like structural protein.

**Fibrosis** is the accumulation of that scar. In dystrophy it is not a side effect; it is a
large part of why patients lose strength. Muscle is literally replaced by scar.

The genes you will see over and over in this project are the collagen genes — `COL1A1`,
`COL1A2`, `COL3A1`, `COL5A1` — plus `LOX`, an enzyme that cross-links collagen fibres to make
the scar stiffer. **When these genes are switched on harder, the muscle is scarring harder.**

A broader term you will meet: the **matrisome**, the full catalogue of proteins that make up the
material *between* cells (the extracellular matrix, or **ECM**). The "core matrisome" is the
structural part of it: collagens, glycoproteins, proteoglycans. We use a published list of 275
such genes.

## 1.4 Genes, expression, and what a microarray measures

A **gene** is a stretch of DNA holding instructions for one protein. Every cell carries all
~20,000 human genes, but uses only some of them. To use a gene, the cell copies it into a
disposable working copy made of RNA, called a **transcript** or **messenger RNA (mRNA)**.

**Gene expression** = how many mRNA copies of a gene are floating in the cell. High expression
means the cell is making a lot of that protein right now.

**A microarray measures expression for every gene at once.** Physically it is a small glass chip
covered in millions of microscopic spots. Each spot is coated with short single-stranded DNA
sequences (**probes**) that match one specific gene. You extract RNA from a muscle biopsy,
label it with a fluorescent dye, and wash it over the chip. RNA sticks to the spot that matches
it. Brighter spot = more of that RNA = higher expression.

Our chip is the **Affymetrix HG-U133 Plus 2.0** (GEO platform code **GPL570**). It carries
**54,675 probesets**. A "probeset" is a group of ~11 probes that all target the same transcript;
the ~11 are averaged to make one number, which is more reliable than any single spot.

Microarrays are older technology than RNA sequencing. That matters here only in one way: it
means our data cannot be merged with modern RNA-seq datasets, which closes off one route to
getting more patients (Part 8).

## 1.5 MicroRNAs: the brakes

Here is the part the project is actually about.

Not all RNA becomes protein. Some RNA molecules are tiny — about **22 letters long** — and their
job is to *switch other genes off*. These are **microRNAs**.

**How a microRNA works, mechanically:**

1. An mRNA transcript has a coding part (the actual protein instructions) and, after it, an
   untranslated tail called the **3′UTR** ("three-prime untranslated region"). The tail carries
   no protein instructions. It is a regulatory noticeboard.
2. A microRNA is loaded into a protein machine and used as a search string.
3. Letters **2 through 8** of the microRNA — seven letters, called the **seed** — scan 3′UTRs
   looking for a complementary match.
4. Where the seed matches, the machine clamps on and the transcript is destroyed or blocked.

```
mRNA 3'UTR:  ... A C G U G G C U G ...
                   | | | | | | |
microRNA seed:     U G C A C C G       <- letters 2-8, pairing with the UTR site
```

**So a microRNA is a brake.** And because the seed is only seven letters long, the same brake
lands on *hundreds or thousands* of different genes. One microRNA tunes an entire programme at
once.

**Seed families.** Different microRNAs that happen to share the same 7-letter seed hit the same
targets, so they are grouped into a **family** and treated as one unit. Our star family is
**miR-29-3p**, seed `AGCACCA`, whose members are hsa-miR-29a-3p, hsa-miR-29b-3p and
hsa-miR-29c-3p. ("hsa" = *Homo sapiens*.)

## 1.6 Why miR-29 specifically

**miR-29 is the brake on collagen.** Its targets are heavily enriched for collagen and other
matrix genes. The logic chains together:

| miR-29 level | Brake strength | Collagen genes | Muscle |
|---|---|---|---|
| Normal | Firm | Held down | Healthy |
| **Reduced** | **Weak** | **Released, rise** | **Scars** |

In Duchenne, miR-29 is **reduced**. Eisenberg et al. 2007 (*PNAS*) reported miR-29c about
**6-fold down** in DMD muscle, alongside COL3A1 about **6-fold up** and COL1A2 about **5-fold
up**. Brake off, collagen up, scar forms. This is established, published, and it is the result
our method has to be able to rediscover before anyone should believe anything else it says.

> **The same paper is also why BMD is an open question.** Eisenberg had only **6 BMD samples**.
> In BMD, only miR-146b and miR-221 passed their fold-change cutoff, and their microRNA–gene
> linkage analysis reached significance only in DMD and in Miyoshi myopathy — **not in BMD**.
> So nobody has shown the miR-29–fibrosis axis operates in human BMD. That gap is our project.

## 1.7 Why anyone should care

Modern therapies for Duchenne — **exon skipping** drugs — work by deliberately deleting an extra
piece of the broken gene so that the reading frame is restored. The goal is explicitly to
**convert a DMD patient into a BMD-like state**: swap no dystrophin for partial dystrophin.

If those therapies succeed, the long-term question becomes: *what happens to fibrosis in a
BMD-like muscle?* Does the scarring programme switch off, or merely slow down? Understanding
what the microRNA–fibrosis axis looks like in naturally occurring BMD is a direct proxy for
what a successfully treated DMD patient might look like years later.

## 1.8 The central trick: footprints

**Our dataset contains no microRNA measurements.** The HG-U133 Plus 2.0 array measures mRNAs.
MicroRNAs are too short for it. We cannot see miR-29 directly.

So we infer it. This is the conceptual core of the entire project:

> **If a brake weakens, everything it was holding down drifts up together.**
> We cannot see the brake. We can see the things it holds. If the whole group of miR-29 targets
> moves up as a group, that is the **footprint** of a weakened miR-29.

Tiny worked example. Suppose miR-29 has 5 targets and we have 3 non-targets:

| Gene | Is a miR-29 target? | Change in disease (log2) |
|---|---|---|
| COL1A1 | yes | +2.5 |
| COL3A1 | yes | +2.3 |
| COL5A1 | yes | +1.9 |
| COL4A1 | yes | +1.2 |
| ELN | yes | +0.8 |
| ACTB | no | +0.1 |
| GAPDH | no | −0.2 |
| TUBB | no | 0.0 |

The targets average +1.74; the non-targets average −0.03. The targets rose **as a group**. That
group drift is the footprint, and it is evidence — not proof — that the miR-29 brake is weaker.

**The direction rule, which you must apply consistently everywhere:**

| What we observe | What it implies about the microRNA |
|---|---|
| Targets bunched at the **UP** end | Brake **weaker**, microRNA **less active** |
| Targets bunched at the **DOWN** end | Brake **stronger**, microRNA **more active** |

Note the inversion: a microRNA that is measured as going **up** should produce targets going
**down**. Getting this backwards flips every conclusion, so check it in every sentence you write.

---

# PART 2 — The research question

## 2.1 What we set out to ask

> **Is the microRNA–fibrosis link that is known to operate in Duchenne also operating in Becker?
> And is Becker's version a weaker copy of Duchenne's, a different pattern, or nothing we can
> detect?**

Three outcomes were possible, and all three were publishable:

| Outcome | Meaning |
|---|---|
| **Weaker copy** | Same microRNAs, smaller effect — BMD is DMD dialled down |
| **Different pattern** | Different microRNAs — BMD is qualitatively its own thing |
| **Nothing detectable** | We cannot tell, at this sample size |

**We landed on the third.** How we got there, and why that is still a real result, is the rest
of this document.

## 2.2 The binding constraint, stated up front

**We have 11 BMD patients.** That is the whole problem. Everything in Part 6 traces back to it.

We checked whether more data existed and the answer was no:

- **GSE38417** — turned out to be a **duplicate re-deposit of GSE109178**, the same raw arrays
  posted under a second accession (Part 6.3). Merging it would have doubled the apparent sample
  size with literally identical numbers and manufactured significance out of nothing.
- **GSE13608** — adds muscle-site and age confounds for roughly 18% more precision. Not worth it.
- **GSE291383** — 6 BMD patients but **no healthy controls**, and it is RNA-seq, which cannot be
  merged with microarray data without introducing a platform effect larger than the signal.

**No public dataset fixes the n = 11 problem.** Say this plainly in the Discussion.

---

# PART 3 — The data

## 3.1 The dataset

| Property | Value |
|---|---|
| Accession | **GSE109178**, NCBI Gene Expression Omnibus |
| Platform | Affymetrix HG-U133 Plus 2.0 (**GPL570**), 54,675 probesets |
| Tissue | **Vastus lateralis** (a thigh muscle) needle biopsies — all samples, same muscle |
| Total samples | **49** |
| Samples we used | **34** |

## 3.2 Who is in it

| Group | n | What it is |
|---|---|---|
| **Control** | **6** | Healthy muscle |
| **BMD** | **11** | Becker |
| **DMD** | **17** | Duchenne |
| LGMD | 15 | **Excluded** — limb-girdle muscular dystrophy, a different disease family (8 LGMD2B / dysferlinopathy + 7 LGMD2I / FKRP). Not dystrophin-related, so irrelevant to our question. |

The 6 / 11 / 17 / 15 counts came out of GEO exactly as the project specification predicted,
which is the first sign the labels were parsed correctly.

## 3.3 Ages and severity — the confound, recorded honestly

| Group | Age median | Age range | Severity breakdown |
|---|---|---|---|
| BMD | **12.2 y** | **2.5 – 60.8 y** | 6 mild, 3 moderate, 2 severe |
| DMD | **3.5 y** | **0.9 – 8.0 y** | 4 moderate, 13 severe |
| Control | **not reported in GEO** `[UNVERIFIED]` | — | — |

**Read this table carefully, because it limits what the paper can claim.** The BMD patients are
roughly a decade older than the DMD patients, and we do not know how old the controls are. So
in the direct BMD-versus-DMD comparison, *disease severity and age are entangled*. Any
difference we found there could be a difference between a 12-year-old and a 3-year-old rather
than between Becker and Duchenne. We **state** this; we do not try to correct it, because with
11 patients an age covariate would consume degrees of freedom we cannot spare.

As it happens we found no BMD-versus-DMD differences at all, so the confound did not get the
chance to mislead us — but it must still appear in Limitations.

## 3.4 One sample with a disputed label

**GSM2934804** is a BMD sample that GEO labels `severe` (`BMD-12-[6.4]-severe`, Deletion Exon 49,
age 6.4), which was questioned. It was carried as a robustness item and tested in step 13.
Dropping it moves miR-29's NES from +2.169 to +2.156. **Immaterial. Item closed.**

## 3.5 Where the input gene lists came from

None of these are things we invented; all are published resources, with versions recorded.

| Input | Source and version | Size |
|---|---|---|
| microRNA seed families | **TargetScanHuman 8.0**, `miR_Family_Info.txt` (dated 2018-04-06), filtered to human (species 9606). Seed = the `Seed+m8` column, positions 2–8. | 2,064 families |
| microRNA → target genes | **TargetScanHuman 8.0**, `Predicted_Targets_Info.default_predictions.txt`, human. These are the **conserved-site default predictions**. | 115,308 rows, 219 families, 13,034 genes |
| 3′UTR lengths | **TargetScanHuman 8.0**, `UTR_Sequences.txt`, human rows; alignment gap characters stripped, longest transcript per gene | 19,429 genes |
| Which microRNAs are in muscle | **miRNA Tissue Atlas 2025** (Center for Bioinformatics, Saarland). Rows with `Organ_system = limb_muscle`. Threshold **mean RPMM ≥ 10** — *our* threshold, not the spec's. | 353 of 2,656 microRNAs pass |
| Fibrosis gene sets | MSigDB 2026.1.Hs via `msigdbr` 26.1.1: **`NABA_CORE_MATRISOME`** (Naba et al. 2012, PMID 22159717) = 275 genes; **`HALLMARK_TGF_BETA_SIGNALING`** = 54 genes; 2 genes overlap | 329 rows |
| Fibroblast markers | MSigDB C8: `RUBENSTEIN_SKELETAL_MUSCLE_FAP_CELLS` + `..._FBN1_FAP_CELLS` (Rubenstein et al. 2020), union | 372 markers; 326 not in a fibrosis set; **316 present on our array** |

**Sanity check that the muscle atlas is real data and not noise:** the highest-expressed
microRNAs it reports in limb muscle are miR-486-5p, miR-10b-5p, miR-133a-3p, miR-1-3p,
miR-22-3p, miR-378a-3p and miR-206. miR-1, miR-133a and miR-206 are the textbook muscle-specific
"myomiRs". The table is behaving exactly as a muscle sample should. miR-29 family levels in
muscle: miR-29a-3p **1,895 RPMM**, miR-29b-3p 650, miR-29c-3p 439 — comfortably expressed.

---

# PART 4 — The statistics, from zero

You need six ideas. Each is defined and then demonstrated on our own numbers.

## 4.1 log2 — why every number is a logarithm

Expression values span an enormous range, and what matters biologically is *ratios*, not
differences. So everything is expressed as **log base 2**.

The only rule you need: **log2 fold change (logFC) of +1 means "doubled".**

| logFC | Linear meaning |
|---|---|
| +1 | 2× |
| +2 | 4× |
| +3 | 8× |
| **+4.058** (our COL1A1 in DMD) | **2^4.058 ≈ 16.7×** |
| **+2.517** (our COL1A1 in BMD) | **2^2.517 ≈ 5.7×** |
| 0 | unchanged |
| −1 | halved |

Why logs: a gene going from 100 to 200 and a gene going from 1,000 to 2,000 are the same
biological event (a doubling), and logs make them the same number (+1). Raw differences would
call them 100 and 1,000 and let the loud gene drown the quiet one.

## 4.2 The t-statistic — "is this difference bigger than the noise?"

Comparing group averages is not enough, because every measurement is noisy. The **t-statistic**
divides the difference by the noise:

```
                difference between group averages
    t  =  ----------------------------------------------
           how much samples bounce around within groups
```

Big t = the difference is large relative to the scatter. Small t = it could easily be chance.
**The sign of t carries the direction** (positive = up in disease).

Worked example with invented numbers. Two genes, both up by exactly 1.0:

| | Difference | Within-group scatter | t |
|---|---|---|---|
| Gene A | 1.0 | 0.2 | **5.0** — convincing |
| Gene B | 1.0 | 2.0 | **0.5** — meaningless |

Same difference, completely different credibility. This is why we rank by t and never by fold
change alone.

## 4.3 The moderated t — borrowing strength across genes (what `limma` does)

With only 6 controls, the scatter estimate for any single gene is itself unreliable. A gene can
look spectacular purely because its 6 control values happened to land close together.

**limma** fixes this by **moderation**: it looks at the scatter of all 20,823 genes, works out
what a typical scatter looks like, and pulls each gene's individual estimate partway toward
that typical value. A gene claiming implausibly tiny noise gets corrected upward and loses its
unearned significance.

The strength of that pull is reported as the **prior degrees of freedom**. Ours is **4.23**,
against **31 residual degrees of freedom** (34 samples − 3 groups). Interpretation: the
moderation is worth roughly 4 extra samples' worth of information per gene — a meaningful but
not overwhelming stabilisation.

**The output we carry forward is the moderated t for every gene in every comparison.** It is the
ranking that everything downstream is built on.

## 4.4 Multiple testing and the FDR

Test 20,823 genes at the usual p < 0.05 threshold and **about 1,041 will pass by pure chance**
even if nothing is happening (20,823 × 0.05). So a raw p-value is uninterpretable at this scale.

The fix is the **false discovery rate (FDR)**, via the Benjamini–Hochberg procedure. It converts
each p-value into an **adjusted p-value** (written `adj.P` or `padj`) with a new meaning:

> **padj < 0.05 means: among all the genes I am calling significant, about 5% are expected to be
> false positives.**

It controls the *proportion* of your claimed findings that are wrong, rather than the chance of
making any error at all. That is the right trade-off for genome-scale screening.

## 4.5 Enrichment and the NES — scoring a whole gene set at once

We do not care about single genes. We care whether miR-29's **1,158 targets, as a group**, sit
unusually high. The tool is **gene set enrichment analysis**, implemented by the `fgsea` package.

**The mechanic, step by step:**

1. Rank all 20,823 genes by moderated t, most-up at the top, most-down at the bottom.
2. Walk down that list from top to bottom carrying a running score.
3. Every time you hit a miR-29 target, the score goes **up**. Every time you hit a non-target,
   it drifts **down** slightly.
4. If targets are concentrated near the top, the score climbs to a high peak early. If they are
   scattered evenly, it wanders around zero.
5. The peak height is the **enrichment score**. Normalising it for set size gives the
   **normalised enrichment score (NES)**.

ASCII sketch of the walk:

```
score
  ^
  |        .-""-.                    <- high peak: targets bunched at the top (NES > 0)
  |      ."      "-.
  |    ."            "--.
  0 --+--------------------"-------->  position in the ranked list
      top                        bottom
```

**Reading the NES:**

| NES | Meaning | What it implies about the brake |
|---|---|---|
| **> 0** | Targets bunched at the UP end | Brake **weaker** |
| ≈ 0 | Targets scattered randomly | No signal |
| **< 0** | Targets bunched at the DOWN end | Brake **stronger** |

Our headline value: **miR-29 NES = +2.224 in DMD versus control.** Targets are strongly bunched
at the up end, so miR-29 is predicted to be *less active* in DMD — which is exactly the
published direction. The method rediscovered known biology.

**The leading edge.** fgsea also reports *which* genes were in the pile that drove the peak.
That list is biologically interpretable: for miR-29 in DMD it is dominated by collagens
(Part 5, step 7).

## 4.6 The most important idea in the entire project: what you permute

A p-value answers "how often would chance alone produce something this extreme?" To answer it,
you need a definition of **chance** — and there are two completely different ones available here.
**They give answers that differ by thirty-plus orders of magnitude.**

### Null A — shuffle the *genes* (what fgsea does by default)

Keep the patients fixed. Keep the ranked list fixed. Pick 1,158 **random genes** instead of
miR-29's targets and ask how often they score NES ≥ 2.224.

**The hidden assumption: that genes are independent of one another.** They are emphatically not.
Collagen genes are co-regulated — they are switched on together by the same signals, so they
rise and fall as a bloc. A random set of 1,158 genes is a scattered sample of unrelated
behaviours; miR-29's targets are a *correlated bloc*. Comparing a correlated bloc against a
scattered-set null makes the bloc look impossibly unlikely.

This null answers: *"Is this gene set special **compared to other gene sets**?"*

### Null B — shuffle the *patients* (what our steps 8 and 9 do)

Keep the gene set fixed. **Randomly reassign the patient labels** — take the 34 samples and deal
out 6 "control", 11 "BMD", 17 "DMD" tickets at random, ignoring who actually has what. Then rerun
the *entire* pipeline: limma, ranking, fgsea. Record the NES. Repeat 500 times.

Now you have 500 NES values produced by a world in which the disease labels are meaningless. If
our real NES sits comfortably inside that pile, **our result is indistinguishable from noise.**

This null answers: *"Could the **disease** have produced this, or would random groups of
patients do just as well?"* — which is the question the paper is actually asking.

### Worked example of why this matters so much

Null A, on our data, says miR-29 in DMD has **padj = 1.5 × 10⁻³⁴**.
Null B, on the *same data, same gene set, same contrast*, says **p = 0.052**.

```
Null A (shuffle genes)    p = 0.00000000000000000000000000000000015
Null B (shuffle patients) p = 0.052
                              ^ not significant at the conventional 0.05 threshold
```

**That is a gap of about 32 orders of magnitude.** Null B is the honest one, because it is the
only one that tests the thing we claim: that *disease status* drives the footprint.

> **Memorise this, because it is the paper.** fgsea's p-values answer a question we were not
> asking, and on correlated biological gene sets they are wrong by an almost unimaginable margin.
> Everything in Part 6 follows from this single distinction.

## 4.7 The shared-control trap

One more design hazard, and it is subtle.

BMD-versus-control and DMD-versus-control use **the same 6 control samples**. So any fluke in
those 6 contaminates *both* comparisons **in the same direction at once**.

Worked example. Two genes with no real disease effect at all, where control sample noise alone
creates the pattern:

| Gene | Ctrl avg | BMD avg | DMD avg | BMD−Ctrl | DMD−Ctrl |
|---|---|---|---|---|---|
| 1 | 9.70 | 10.05 | 9.98 | **+0.35** | **+0.28** |
| 2 | 10.30 | 9.96 | 10.02 | **−0.34** | **−0.28** |

Gene 1's controls happened to read low, so it looks up in *both* diseases. Gene 2's controls
happened to read high, so it looks down in *both*. Plot BMD−Ctrl against DMD−Ctrl and you get a
beautiful straight line through the origin — **manufactured entirely by 6 noisy controls.**

**Consequence:** the correlation between BMD's and DMD's footprints is **guaranteed** to look
high. It can never be judged against zero. It must be judged against what shuffled labels
produce — which is step 8, and which is where this project's headline came from.

**And the clean escape hatch:** compare BMD against DMD **directly**, patient to patient,
throwing the controls out of the test entirely. Then no control fluke can possibly create a
difference. That is step 9, and it is the comparison that actually answers "is BMD different
from DMD?"

---

# PART 5 — Every step, in order

Each step below gives: **what it does**, **why it exists**, **what came out**, and **how to read
it**. Scripts are in `scripts/`, outputs in `results/`.

---

## Step 0 — Setup (`00_setup.R`)

Installed the software stack and recorded it. Written to `results/sessionInfo.txt`.

| Component | Version |
|---|---|
| R | **4.6.2** (2026-10-29) |
| Bioconductor | **3.23** |
| limma | 3.68.5 |
| fgsea | 1.38.0 |
| msigdbr | 26.1.1 (MSigDB 2026.1.Hs) |
| BiocParallel | 1.46.0 |
| data.table | 1.18.6.1 |
| ggplot2 | 4.0.3 |
| others | oligo, pd.hg.u133.plus.2, hgu133plus2.db, AnnotationDbi |

**Deviations to disclose:** `GEOquery` failed to build on this machine and was never installed.
It was not needed — see step 1. Every package was compiled **from source**, because Homebrew's
R build does not use CRAN's macOS binaries.

---

## Step 1 — Download and labels (`01_download_labels.R`)

**What it does.** Fetches the raw array files and the sample metadata from GEO, then builds the
sample sheet that says who is in which group.

**Deviation, and it is worth a Methods sentence.** The standard route, `getGEOSuppFiles()`, was
unavailable (no GEOquery), and NCBI's FTP server returns **HTTP 403** for file downloads from
this machine. We used GEO's two official web endpoints instead:

- metadata: `www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109178&targ=gsm&form=text&view=brief`
- raw archive: `www.ncbi.nlm.nih.gov/geo/download/?acc=GSE109178&format=file`

These serve identical content, and the fields parsed are exactly the ones `pData()` would expose.
Only the 34 included CEL files were extracted from the 49-file archive.

**Result.** 6 control / 11 BMD / 17 DMD / 15 LGMD — matching the specification exactly. Ages and
severities as in Part 3.3. **Control ages are absent from GEO** and are marked `[UNVERIFIED]`.

**Output:** `results/tables/01_sample_sheet.tsv`

---

## Step 2 — Normalisation and quality control (`02_rma_qc.R`)

**What it does.** Converts raw fluorescence intensities into comparable expression values, then
checks that no array is technically broken.

**Why normalisation is needed.** Arrays are processed on different days with slightly different
amounts of material, so one chip can be globally brighter than another for purely technical
reasons. **RMA** (Robust Multi-array Average) corrects this in three moves: subtract background,
force every array onto the same overall intensity distribution (quantile normalisation), then
summarise each probeset's ~11 probes into one value.

**A deliberate choice worth defending:** we ran RMA on **only the 34 included samples**, not all
49. Quantile normalisation lets every sample influence every other sample's values, so including
the 15 LGMD arrays would have let a different disease shape our numbers.

**Output matrix: 54,675 probesets × 34 samples.**

### The QC, and a methodological disagreement we resolved

The specification asked for "each sample's median correlation with all the others", flagging
anything unusually low. **That metric fails on this dataset**, and understanding why is a nice
Methods paragraph.

28 of our 34 arrays are dystrophic muscle; only 6 are healthy. A healthy sample is therefore
being correlated against 28 samples that are biologically unlike it and only 5 that are like it.
**All 6 controls are guaranteed to score lowest — by group size, not by quality.** On that
metric three controls looked like outliers and would have been thrown out, which would have been
a serious error.

So we added two better measures:

1. **Within-group median correlation** — each sample against others *in its own group*, removing
   the group-size effect.
2. **RLE and NUSE** — the standard Affymetrix array-quality statistics, computed from
   `oligo::fitProbeLevelModel`. Crucially these are calculated **within each probeset's own ~11
   probes**, so patient biology cannot inflate them. They test the *chip*, not the patient.

| Metric | Result | Flag threshold | Verdict |
|---|---|---|---|
| Median RLE, all 34 arrays | within **±0.057** | 0.10 | none flagged |
| Median NUSE, all 34 arrays | **0.981–1.039** | 1.05 | none flagged |
| Overall-correlation by group | control 0.9529, BMD 0.9704, DMD 0.9803 | — | controls lowest — the artefact |
| **Within-group** by group | control **0.9839**, BMD 0.9757, DMD 0.9855 | — | **controls highest** |

Look at the last two rows together. On the flawed metric controls rank **lowest**; on the
corrected metric the very same controls rank **highest**. The "outliers" were an artefact of
counting.

> ### ✅ DECISION (Aziz, 2 Nov 2026): keep all 34 samples. This STOP point is closed.
>
> RLE and NUSE — the measures that actually test array quality — flag nothing. Two samples
> remain mildly flagged on the within-group metric (**GSM2934803**, BMD, 0.9676 vs a 0.9757
> group centre; **GSM2934813**, control, 0.9784 vs 0.9839) and were **kept**, with step 13
> testing directly whether either was carrying the result. Neither was.

**Outputs:** `results/figures/02_qc_boxplot.png`, `02_qc_nuse.png`, `02_qc_within_group_cor.png`,
`02_qc_median_cor.png`; `results/tables/02_qc_metrics.tsv`

---

## Step 3 — Probes to genes (`03_annotate_collapse.R`)

**What it does.** Translates probeset IDs (`1007_s_at`) into gene symbols (`COL1A1`), and
reduces the many probesets that measure the same gene down to one row per gene.

**Why.** Our microRNA target lists are written in gene symbols, so the two sides must speak the
same language. And a gene measured by 5 probesets would otherwise get 5 votes in the enrichment
walk instead of 1.

**The funnel:**

```
54,675 probesets on the chip
   ->  44,647  have a gene symbol in hgu133plus2.db
   ->  -1,524  dropped: map to more than one gene (ambiguous, unusable)
   ->  43,123  kept
   ->  20,823  genes after collapsing
```

11,043 genes had more than one probeset. **Collapse rule: keep the probeset with the highest
mean expression.** Rationale: among probesets for the same gene, the brightest is generally the
one that hybridises properly; dim ones are often badly placed or cross-hybridising, and sit in
the noise floor where ratios are unstable.

**No low-expression filter was applied** (the project default). The resulting t-statistics span
−15.0 to +17.7, a healthy range with no sign of noise domination.

**Output: a 20,823 genes × 34 samples matrix.**

---

## Step 4 — Differential expression with limma (`04_limma.R`)

**What it does.** For every one of the 20,823 genes, in each of three comparisons, computes how
much it changed and how trustworthy that change is.

**Design:** `~0 + group`, giving one coefficient per group, then three contrasts:

| Contrast | Question |
|---|---|
| **BMD − Control** | What does Becker do to this gene? |
| **DMD − Control** | What does Duchenne do to this gene? |
| **BMD − DMD** | Do the two diseases differ on this gene? (controls cancel out) |

Pipeline: `lmFit` → `contrasts.fit` → `eBayes`. Residual df **31**, prior df **4.23**.

### How many genes changed

| Contrast | Up at FDR < 0.05 | Down at FDR < 0.05 |
|---|---|---|
| BMD − control | 4,549 | 4,659 |
| DMD − control | 5,547 | 6,654 |
| BMD − DMD | 4,564 | 4,039 |

**Roughly 9,000–12,000 of 20,823 genes change.** That is an enormous fraction, and it is not a
bug — dystrophic muscle is grossly remodelled tissue, full of scar, immune infiltrate and
regenerating fibres. **But it is a warning**, and it comes back to bite us in step 6: when half
the transcriptome is moving, almost *any* gene set will look enriched.

### The fibrosis sanity check — does the method see known biology?

Before trusting anything, confirm the collagen genes behave as a century of pathology says they
must. `results/tables/04_sanity_fibrosis_genes.tsv`, in log2 fold change:

| Gene | BMD−Ctrl | adj.P | DMD−Ctrl | adj.P | BMD−DMD | adj.P |
|---|---|---|---|---|---|---|
| COL1A1 | +2.517 | 9.4e−08 | **+4.058** | 5.6e−13 | −1.542 | 4.6e−06 |
| COL1A2 | +2.043 | 1.4e−07 | +2.940 | 7.5e−12 | −0.897 | 4.2e−04 |
| COL3A1 | +2.311 | 8.3e−09 | +3.260 | 4.6e−13 | −0.949 | 1.3e−04 |
| COL5A1 | +1.919 | 1.4e−06 | +2.828 | 7.5e−11 | −0.910 | 8.5e−04 |
| LOX | +2.607 | 2.4e−07 | +3.573 | 3.7e−11 | −0.966 | 3.1e−03 |

**Three things to read off this table:**

1. **All five are up in both diseases.** Both BMD and DMD muscle is scarring. Good — expected.
2. **All five are up *less* in BMD than in DMD**, and significantly so. COL1A1 is **2^4.058 ≈
   16.7×** in Duchenne but **2^2.517 ≈ 5.7×** in Becker. Becker scars, but roughly a third as hard.
3. **It matches the literature.** Eisenberg 2007 reported COL3A1 about 6× and COL1A2 about 5× up
   in DMD; we get 2^3.260 ≈ 9.6× and 2^2.940 ≈ 7.7×. Same direction, same order of magnitude,
   different cohort and platform.

> **This is one of the few hard, defensible, per-gene results in the whole project, and it
> survives everything.** It is a limma t-test on individual genes — it involves **no permutation
> of anything** and so is untouched by the fgsea problem in Part 6.1. The BMD-versus-DMD collagen
> difference is **real and reportable**. What is *not* reportable is the claim that the
> *microRNA footprint* differs (step 9). Keep those two statements rigorously apart.

**Outputs:** `results/tables/limma_BMDvsCtrl.tsv`, `limma_DMDvsCtrl.tsv`, `limma_BMDvsDMD.tsv`,
`04_sanity_fibrosis_genes.tsv`. **This replaces the old preliminary GEO2R claim, which must not
be cited anywhere.**

---

## Step 5 — Building the microRNA target sets (`05_mirna_inputs.R`)

**What it does.** Turns TargetScan's raw files into the 128 gene sets we will test.

**The filtering funnel, with the reason for each cut:**

```
2,064  human TargetScan 8.0 seed families
  ->  253  expressed in limb muscle (miRNA Tissue Atlas 2025, mean RPMM >= 10)
  ->  128  with at least 15 targets present on our array
```

**Why filter to muscle.** A microRNA that is not present in muscle cannot be regulating anything
in muscle. Testing it would add a meaningless hypothesis and dilute the multiple-testing
correction, making real signals harder to see.

**Why require ≥ 15 targets.** Enrichment scores on tiny sets are wildly unstable — with 4 genes,
one outlier dominates. 15 is the project's agreed floor. **No upper limit** was imposed, which is
why some sets are very large.

**Result: 128 families, sizes 15 to 1,446, median 371.**

**miR-29-3p specifically:**

| Property | Value |
|---|---|
| Seed (positions 2–8) | **`AGCACCA`** |
| Members | hsa-miR-29a-3p, hsa-miR-29b-3p, hsa-miR-29c-3p |
| TargetScan conservation class | 2 |
| Muscle expression | 1,895 RPMM (miR-29a-3p) |
| **Targets present on our array** | **1,158 genes** |

> **Note the set sizes and remember them.** miR-29 has 1,158 targets out of 20,823 genes — that
> is **5.6% of the measured transcriptome in a single "gene set"**. Sets this large, drawn from
> overlapping prediction lists, are exactly the conditions under which gene-permutation p-values
> break down. This is foreshadowing.

**Outputs:** `results/rds/family_sets.rds`, `results/tables/05_family_sizes.tsv`

---

## Step 6 — The smoke-alarm test (`06_smoke_test.R`)

**What it does.** The go/no-go gate for the whole project: run fgsea on DMD-versus-control and
check that **miR-29 comes out positive**. miR-29's involvement in Duchenne fibrosis is
established. If our method cannot rediscover a known result in the group where it is known to be
true, the method does not work and nothing downstream could be believed.

### ✅ VERDICT: GO

| Measure | Value |
|---|---|
| miR-29-3p NES | **+2.224** |
| p | 3.4e−36 |
| padj | **5.4e−35** |
| Rank among 128 families | **3rd by padj** (2nd by NES — see the note below) |
| Direction | NES > 0 → targets up → brake weaker → miR-29 *less active* in DMD — **the published direction** |
| Mean log2FC of the 1,158 targets | **+0.304** |
| Mean log2FC of the other 19,665 genes | **+0.036** |
| **The "ruler"** (footprint size; used by step 10b) | **+0.269 log2** ≈ **20.5%** on the linear scale |

**How to read the ruler.** miR-29's targets rose, on average, 0.269 log2 units *more* than
everything else. Converting: 2^0.269 ≈ 1.205, so the targets moved about **20.5% further up**
than the typical gene. That single number becomes the measuring stick in step 10b.

### ⚠ The caveat that must travel with this result everywhere

**103 of 128 families reached padj < 0.05 in this comparison, and 99 of them point the same way
(NES > 0).**

That is not a subtle problem. Essentially *every* muscle microRNA family looks "significantly
de-repressed" in Duchenne. The reason is step 4's finding: when ~11,000 of 20,823 genes shift
upward, a *competitive* rank test will report almost any gene set as enriched, because the
whole list it is walking down is tilted.

> **So the smoke test proves the method *can* recover known biology. It does not, on its own,
> show that miR-29 is special.** Those are different claims and only the first is supported here.
> The specificity evidence comes from steps 9, 10a, 11 and 12, where miR-29 is compared against
> properly matched baselines rather than against zero.

**Outputs:** `results/tables/06_smoke_fgsea_DMDvsCtrl.tsv`,
`results/figures/06_mir29_enrichment_DMDvsCtrl.png` (the ranked-walk plot from Part 4.5)

**A small discrepancy you will notice and must be able to explain:** step 6 reports padj =
5.4e−35 and step 7 reports 1.5e−34 for the *same* comparison. fgsea estimates p-values by random
permutation, so two runs with different random draws give slightly different numbers at these
extremes. Both are "astronomically small and meaningless" for the reason in Part 6.1. Quote
step 7's value for consistency with the main tables.

---

## Step 7 — All families, all three comparisons (`07_fgsea_all.R`)

**What it does.** Runs the same enrichment test for all 128 families across all three contrasts,
producing the master results table.

### The result, and the shape that drives the paper

| Contrast | padj < 0.05 | of which NES > 0 | NES < 0 |
|---|---|---|---|
| BMD − control | 104 / 128 | 99 | 5 |
| DMD − control | 103 / 128 | 99 | 4 |
| **BMD − DMD** | **7 / 128** | **0** | **7** |

**Look at that third row.** Against controls, ~103 of 128 families are "significant". In the
direct patient-to-patient comparison, only 7 are. **Removing the shared controls collapses the
significance count by a factor of fifteen.** That is a measure of how much of the
against-controls signal was a global tilt rather than microRNA-specific biology.

### miR-29 across the three

| Contrast | NES | padj | Rank **by padj** | Rank **by NES** |
|---|---|---|---|---|
| BMD − control | +2.169 | 1.4e−31 | 6 / 128 | 11 / 128 |
| DMD − control | +2.224 | 1.5e−34 | 3 / 128 | 2 / 128 |
| **BMD − DMD** | **−1.590** | **4.7e−10** | **1 / 128** | — |

> **⚠ Always say which ranking you mean.** These two orderings disagree, because padj depends on
> set size as well as effect size. In DMD-versus-control miR-29 is **3rd by padj but 2nd by NES**;
> in BMD-versus-control it is **6th by padj but 11th by NES**. The verified top 3 by padj in
> DMD-versus-control are miR-30-5p (2.3e−36), miR-1-3p/206 (2.7e−36), miR-29-3p (1.5e−34).

In the direct comparison miR-29 is **rank 1**, about **five orders of magnitude** ahead of the
runner-up (miR-30-5p at 1.9e−04; miR-1-3p/206 at 4.1e−04).

**Direction check:** NES < 0 for BMD − DMD means miR-29's targets sit *lower* in BMD than in DMD,
so the miR-29 brake is **less weakened in Becker than in Duchenne** — the milder-disease
direction, exactly as biology predicts. *(It will fail the honest test in step 9. Do not get
attached to it.)*

### The leading edge — what is actually driving the signal

The genes in miR-29's leading edge for DMD-versus-control that are also in the fibrosis sets:

> COL15A1, COL16A1, COL1A1, COL1A2, COL21A1, COL27A1, COL3A1, COL4A1, COL4A2, COL4A4, COL5A1,
> COL5A2, COL5A3, COL6A2, COL6A3, COL6A6, CRISPLD1, ELN, FBN1, FRAS1, HMCN1, LAMA2, LAMC1,
> MFAP2, MFAP3, NID1, PXDN, SKI, SPARC, SVEP1

**That is a roll-call of collagens.** The enrichment is not being driven by a random assortment
of genes that happen to carry a seed match — it is being driven by precisely the scar programme
miR-29 is published to control.

**miR-29's leading edge is 6.5% fibrosis-set genes — the highest of any top family** (next
highest 4.1%; most are 1–3%).

> **This is specificity evidence of a completely different kind, and you should lean on it in the
> Discussion.** It is about *content*, not p-values. It cannot be inflated by permutation
> problems, correlated gene structure, or sample size, because it is simply a statement about
> which genes are in the list. **It survives the entire collapse described in Part 6.**

### The agreement between BMD and DMD footprints

| Measure across 128 families | Value |
|---|---|
| Pearson r (BMD NES vs DMD NES) | **0.9123** |
| Slope (BMD NES on DMD NES) | **0.9769** |

An r of 0.91 looks like overwhelming evidence that BMD and DMD have the same footprint, and a
slope near 1 looks like they have it at the same strength.

> **🚫 Do not quote either number against zero.** Both are measured against the same 6 controls,
> so Part 4.7's trap applies in full. Their meaning is decided entirely by step 8.

**Outputs:** `results/tables/07_fgsea_wide.tsv`, `07_leading_edge.tsv`,
`07_leading_edge_summary.tsv`; `results/figures/07_nes_heatmap.png`

---

## Step 8 — The shuffled-label baseline (`08_shuffle_baseline.R`) — **the step that changed the paper**

**What it does.** Establishes what r and the slope look like when the disease labels are
**meaningless**, so that step 7's r = 0.9123 can be judged against the right baseline instead of
against zero.

**The procedure, one iteration:**

1. Take the 34 samples and randomly deal out 6 "control", 11 "BMD", 17 "DMD" tickets, ignoring
   who actually has what.
2. Rerun limma on those fake groups, both control contrasts.
3. Rerun fgsea over all 128 families.
4. Record the resulting r and slope.

Repeat **500 times**. (The plan said 1,000; the measured throughput made that 5.7 hours for this
step alone, and the project rules explicitly permit 500. Logged. The empirical p floor is
therefore 1/501 = 0.002, ample resolution for our purposes.) **Runtime: 206.7 minutes.**

### The headline result

| | Real | Null mean | Null median | Null SD | **Null 95th pct** | Null max | **Empirical p** |
|---|---|---|---|---|---|---|---|
| **Pearson r** | **0.9123** | 0.510 | 0.553 | 0.345 | **0.933** | 0.984 | **0.0639** |
| **Slope** | **0.9769** | 0.516 | 0.555 | 0.404 | 1.094 | 1.311 | **0.1138** |

**Read the r row very slowly, because this is the moment the project turned over.**

Randomly shuffled, biologically meaningless labels produce a correlation of **0.933 or higher on
5% of runs**, and reached **0.984** at worst. Our real value is **0.9123**.

```
      null distribution of r from 500 shuffles
      |
  0.0 -+----------########################-------+- 1.0
                   ^        ^       ^     ^
                   |        |       |     null max 0.984
                   |        |       null 95th pct 0.933
                   |        REAL r = 0.9123   <- inside the noise
                   null median 0.553
```

**Our real r does not even clear the 95th percentile of pure noise.** The empirical p is
**0.064**. The slope is worse: p = **0.114**.

> ### 🚫 `r = 0.91` is not evidence of anything.
> Neither is a slope near 1. **Both are simply what this study design produces when nothing is
> going on.** The project specification predicted r ≈ 0.7 from noise alone with these group
> sizes; the real null median is 0.553 and its upper tail reaches 0.984, so the warning was
> qualitatively right and if anything understated the problem.

### The free by-product: a patient-shuffle p for every family

Each shuffle also produces a null NES for every family, which gives each family an honest
**patient-shuffle p** to sit beside its fgsea padj. **This is where the project's central finding
became visible:**

| | fgsea padj | **Patient-shuffle p** |
|---|---|---|
| miR-29, **DMD** − control | 1.46e−34 | **0.0522** |
| miR-29, **BMD** − control | 1.44e−31 | **0.0843** |

**Same data. Same gene set. Same contrast. Two nulls. About 32 orders of magnitude apart.**
Neither honest p-value clears 0.05.

### Counting families, against what chance would give

With 128 families tested at p < 0.05, **6.4 would be expected to pass by chance alone**.

| Contrast | Families with shuffle p < 0.05 | Versus 6.4 expected |
|---|---|---|
| DMD − control | **1** of 128 | **below** chance |
| BMD − control | **10** of 128 | barely above chance |

The 10 in BMD: miR-409-3p (0.018), miR-376-3p (0.026), miR-382-5p (0.032), miR-196-5p (0.034),
miR-136-5p (0.034), miR-193-3p (0.038), miR-186-5p (0.042), miR-425-5p (0.042), miR-145-5p
(0.046), miR-423-5p (0.050).

> **🚫 Do not present that list as a finding.** None survives any correction for 128 tests, and
> 10 out of 128 is not meaningfully more than the 6.4 that noise delivers for free. *(Several of
> the first few are described in the literature as members of the 14q32 / DLK1-DIO3 imprinted
> cluster, which would be interesting if it were real — but this is `[UNVERIFIED]`, unchecked
> against any source, and the counts are at chance level. Leave it out.)*

### Where miR-29 ranks on the honest test — **a precision point you must get right**

| | Rank by NES | Rank by patient-shuffle p |
|---|---|---|
| **DMD** − control | 2nd of 128 | **3rd of 128** |
| **BMD** − control | 11th of 128 | **45th of 128** |

> **⚠ The claim "miR-29 ranks 3rd of 128 even on the honest test" is true for DMD-versus-control
> ONLY.** In BMD-versus-control it ranks **45th**. The DMD ranking is what supports the
> interpretation "underpowered rather than absent". The BMD ranking does not support it. Never
> let "3rd of 128" sit unqualified next to a sentence about Becker — a reviewer will check.

**Outputs:** `results/tables/08_family_shuffle_p.tsv`, `08_shuffle_log.txt`;
`results/figures/08_shuffle_r.png`, `08_shuffle_slope.png`

---

## Step 9 — The direct BMD-versus-DMD comparison (`09_bmd_vs_dmd.R`, `09b_classify.R`)

**What it does, and why it is the comparison that matters.** Steps 7 and 8 both measure BMD and
DMD against the same 6 controls, so the shared-control trap applies. **This step throws the
controls away entirely** and compares the 28 patients against each other, shuffling the BMD/DMD
label among them keeping 11/17.

With no controls in the test at all, **a fluke in the controls cannot possibly create a
BMD–DMD difference.** This is the clean comparison, and it is where "is BMD different from DMD?"
is actually answered.

### The answer is no — for every single family

| | Value |
|---|---|
| Families with fgsea padj < 0.05 | **7** of 128 |
| Families with **patient-shuffle p < 0.05** | **0** |
| Families passing **both** | **0** |

**miR-29 specifically:**

| | Value |
|---|---|
| NES, BMD − DMD | **−1.590** |
| fgsea padj | 4.70e−10 |
| **Patient-shuffle p** | **0.2394** |
| Shuffle q (BH-corrected over 128) | 0.951 |
| Null NES mean / SD | 0.053 / 1.386 (from 496 usable shuffles; 4 returned non-finite NES) |

On this single family, fgsea is overconfident by a factor of roughly **2.3 × 10⁸**.

**Note the direction, then note that it does not matter.** NES < 0 means miR-29 targets sit
lower in BMD, i.e. the brake is predicted *more active* in Becker — the milder-disease direction
you would hope for. **It is not significant (p = 0.24).** Report the direction only as a
non-significant observation, if at all.

### Shared / BMD-only / DMD-only — **the single cleanest table in the project**

> **Run note, for the Methods.** The classification block at the end of `09_bmd_vs_dmd.R` needs
> step 8's output file. Step 9 ran on 3 Nov at 01:49; step 8 did not finish writing that file
> until 05:46. The block therefore died with a file-not-found error and **no classification table
> was ever produced** — this went unnoticed until 6 Nov. The 500 shuffles themselves had
> completed and were saved, so `09b_classify.R` rebuilt the classification from the saved tables
> with **no recomputation of anything**. Logged in `docs/methods_log.md`.

Because **0 of 128 families differ** between BMD and DMD, no family can possibly be BMD-only or
DMD-only. Every family lands in "shared" or "neither" — and *which one* depends entirely on
which yardstick you use for "altered versus controls":

| "Altered vs controls" judged by | shared | neither | BMD-only | DMD-only |
|---|---|---|---|---|
| **(a) fgsea padj < 0.05** | **101** | 27 | 0 | 0 |
| **(b) patient-shuffle p < 0.05** | **1** (miR-196-5p) | 127 | 0 | 0 |

> ### This one table is the whole paper in six numbers.
>
> Asked *"how many microRNA families are altered in dystrophic muscle?"*, on **identical data**:
> - **fgsea answers 101 of 128.**
> - **The patient shuffle answers 1 of 128** — and 1 is *fewer* than the 6.4 you would expect
>   from chance alone, so even that single survivor is unremarkable.

miR-29 is classified **"shared"** under (a) and **"neither"** under (b).

### The sentence to write

> Comparing BMD and DMD patients directly, with healthy controls excluded so that control-sample
> variation cannot contribute, no microRNA family's target footprint differed between the two
> groups under patient-label permutation (0 of 128 families, all q > 0.95; miR-29 p = 0.24). We
> therefore cannot distinguish BMD's predicted microRNA activity from DMD's at these sample
> sizes, and we report no shared, BMD-specific or DMD-specific families.

**Outputs:** `results/tables/09_bmd_vs_dmd.tsv`, `09_family_classification.tsv`, `09_log.txt`,
`09_classify_log.txt`; `results/figures/09_mir29_bmd_vs_dmd_shuffle.png`

---

## Step 10a — Is the footprint just the scar genes? (`10a_ecm_removed.R`)

**The worry.** miR-29's targets are full of collagens. Dystrophic muscle is full of scar. So
maybe the "miR-29 footprint" is nothing more than "this biopsy contains a lot of scar tissue",
with no microRNA regulation involved at all.

**The test.** Delete all **275 core-matrisome genes** from *both* the ranked list **and** the
family sets, then rerun. (Deleting from both matters — otherwise set sizes would be wrong and
the scores incomparable. Families dropping below 15 genes after deletion were excluded.)

| | NES before | NES after | padj after | Survives? |
|---|---|---|---|---|
| miR-29, DMD − control | +2.224 | **+2.151** (96.7% retained) | 1.4e−29 | ✅ |
| miR-29, BMD − control | +2.169 | **+2.109** | 1.0e−27 | ✅ |

Across all families the median NES barely moved (1.679 → 1.684).

> **Conclusion: the footprint is not merely the scar genes.** Strip out every structural matrix
> gene and **96.7% of miR-29's enrichment remains**. Whatever is driving the signal, it is
> distributed across miR-29's non-matrisome targets too. This is genuine specificity evidence
> and it survives Part 6 intact.

**Output:** `results/tables/10a_ecm_removed.tsv`

---

## Step 11 — The 3′UTR-length-matched null (`11_utr_matched.R`) — **the most sobering step**

**The problem, and it is a deep one.** Genes with long 3′UTRs have more room for seed matches, so
they appear on *many* microRNA target lists at once. They also tend to behave alike in disease.
So comparing miR-29's targets against *random* genes is rigged in our favour: we would be
comparing long-tailed genes against a mixture of long and short, and long-tailed genes drift up
anyway.

**The fix.** For each family, bin all genes into **deciles of 3′UTR length**, then draw 1,000
random sets that copy that family's per-bin counts exactly. If miR-29 has 200 targets in the
longest decile, every random comparison set also gets 200 genes from the longest decile. Score
them identically.

### The result: a UTR-matched random set already scores NES ≈ +1.6. Not zero.

| | Real NES | UTR-matched null mean | Null SD | **z** | Empirical p |
|---|---|---|---|---|---|
| miR-29, DMD vs control | 2.222 | **1.592** | 0.096 | **6.56** | 0.000999 |
| miR-29, BMD vs control | 2.172 | **1.669** | 0.093 | **5.39** | 0.000999 |

*(p = 0.000999 is the floor, 1/1001 — **not one** of the 1,000 matched sets ever reached the real
value. The real NES differs very slightly from step 7's because it is recomputed on the
UTR-restricted gene universe, so that it is comparable to the nulls.)*

**How to read this, in two halves — and you must write both halves.**

**The bad half:** about **72%** of miR-29's apparent enrichment (1.59 out of 2.22) is explained
by **3′UTR length alone**, with no microRNA biology involved whatsoever. The raw NES of 2.22
badly overstates the real effect, because the implicit baseline of "zero" is wrong — the correct
baseline is about 1.6.

**The good half:** the remaining excess is **6.5 standard deviations** above a properly matched
baseline, and **not one of 1,000 matched random sets came close.** The footprint is real.

> **🚫 Never quote an NES as though zero were the baseline.** For a gene set of this size and UTR
> profile, the honest baseline is ≈ +1.6.

**All 25 top families survived at q < 0.05**, so this step validates the method broadly rather
than singling out miR-29. By z-score miR-29 ranks **2nd** in DMD-versus-control, behind
miR-1/206 (the canonical muscle microRNAs — reassuring, since muscle loss is the other thing
unambiguously happening in these biopsies).

**Output:** `results/tables/11_utr_matched.tsv`

---

## Step 12 — Adjusting for fibroblast content (`12_fibroblast_cov.R`)

**The worry, which is different from step 10a's.** A dystrophic biopsy contains fewer muscle
cells and more **fibroblasts** (the cells that make scar). So part of any expression difference
is not "muscle cells changed their behaviour" but simply "this sample contains different cells".
That is a **tissue composition** effect, not a regulation effect.

**The test.** Build a per-sample fibroblast score as the mean row-z-score of **316 array-present
fibroblast markers that are *not* themselves fibrosis-set genes** (using fibrosis genes would be
circular). Then refit limma as `~0 + group + fibro` and rerun everything.

### The score behaves exactly as it should

| Group | Mean fibroblast score |
|---|---|
| Control | **−0.791** |
| BMD | **+0.040** |
| DMD | **+0.253** |

ANOVA on group: **p = 1.2e−08**. Healthy muscle has the fewest fibroblasts, Duchenne the most,
**and Becker sits neatly in between** — a quantitative confirmation that BMD is intermediate
tissue, obtained independently of any microRNA analysis.

### What survives the adjustment

| Contrast | Significant before | After | Survive |
|---|---|---|---|
| BMD − control | 104 | 80 | 74 |
| DMD − control | 103 | **51** | 45 |
| BMD − DMD | 7 | 3 | 3 |

**miR-29 survives all three:**

| Contrast | NES raw → adjusted | padj adjusted |
|---|---|---|
| BMD − control | +2.169 → +1.596 | 1.1e−10 |
| DMD − control | +2.224 → **+1.762** (**79% retained**) | **1.7e−16** |
| BMD − DMD | −1.590 → −1.300 | 7.1e−03 |

> **Second specificity signal.** Adjusting for fibroblast content **halves** the number of
> significant families (103 → 51), but miR-29 retains **79%** of its NES where other top families
> retain only ~55–60% — and afterwards it leads the field by **eleven orders of magnitude**
> (1.7e−16 versus 3.1e−05 for the next family). Whatever miR-29 is picking up, it is
> substantially *not* explained by how much scar tissue is in the biopsy.
>
> **Required caveat (and do include it):** scar load genuinely differs by group, so this
> covariate also absorbs real disease signal along with the composition effect. Read this step
> as *"does the footprint survive?"* — not as the headline effect estimate.

**Outputs:** `results/tables/12_fibro_adjusted_fgsea.tsv`, `12_fibroblast_scores.tsv`;
`results/figures/12_fibroblast_score.png`

---

## Step 13 — Leave-one-out (`13_leave_one_out.R`) — 7.3 min

**What it does.** With 11 BMD and 6 control samples, a single unusual patient could in principle
be carrying the entire result. So: drop each BMD sample in turn (11 runs) and each control in
turn (6 runs), **17 complete re-runs of the whole pipeline**, and see whether the answer moves.

| | Full data | Range across 17 re-runs | Significant in |
|---|---|---|---|
| miR-29, BMD vs control | +2.169 | **+2.120 to +2.186** | **17/17** |
| miR-29, BMD vs DMD | −1.590 | **−1.674 to −1.500** | **17/17** |

miR-29 remained **rank 1 of 128 in every single BMD-versus-DMD run**. The top-20 family list
overlapped the full-data list by **18–20 of 20** in every run.

**This closes the step-2 QC question definitively.** The two samples that were mildly flagged on
the within-group correlation metric are harmless:

- **GSM2934803** (BMD) — the largest single-sample effect of any, and *still* 18/20 overlap with
  miR-29 at rank 1.
- **GSM2934813** (control) — 19/20 overlap.

**GSM2934804** (the disputed-severity sample from Part 3.4) is also immaterial: dropping it gives
NES +2.156 (versus +2.169) and −1.507 (versus −1.590), rank 1 either way. **Robustness item
closed.**

> **What this step does and does not prove.** It proves the result is not an artefact of one
> patient — it is a property of the group. **It does not rescue the result from step 8**, because
> every one of those 17 re-runs used fgsea's gene-permutation p-value. Stability and validity are
> different things: a biased measurement repeated 17 times gives 17 consistent biased answers.

**Outputs:** `results/tables/13_leave_one_out.tsv`, `13_mir29_leave_one_out.tsv`;
`results/figures/13_mir29_leave_one_out.png`

---

## Step 14 — Is BMD weaker only because n is smaller? (`14_subsample_dmd.R`) — 31.1 min

**The question.** BMD's miR-29 NES is 2.169; DMD's is 2.224. Is BMD's footprint genuinely
weaker, or does it only *look* weaker because we have 11 Becker patients and 17 Duchenne ones?

**The test.** Randomly cut DMD down to 11 patients, rerun, repeat **100 times**. That produces
the distribution of NES values Duchenne itself yields *at Becker's sample size*. Then see where
Becker's real value falls in it.

| | miR-29 NES vs control | Families significant |
|---|---|---|
| Full DMD (n = 17) | 2.224 | 103 / 128 |
| **DMD cut to 11**, 100 draws | **median 2.167** (range 2.074–2.258) | median 102 (98–107) |
| **Actual BMD (n = 11)** | **2.169** | 104 / 128 |

**53.0% of the DMD-at-n=11 draws land at or below BMD's value.** That is a coin flip. BMD sits
almost exactly at the median of the distribution Duchenne produces at matched sample size.

```
   DMD subsampled to n = 11, 100 draws
   2.07 |-----------########################------| 2.26
                            ^
                            BMD's actual value, 2.169
                            53% of DMD draws are at or below this
```

> ### 🚫 Therefore: measured against healthy controls, BMD's miR-29 footprint is NOT detectably
> weaker than DMD's. The apparent gap (2.169 vs 2.224) is entirely sample-size noise.

### Why this does not contradict the collagen genes in step 4

This looks paradoxical at first — step 4 showed COL1A1 at +2.52 in BMD versus +4.06 in DMD, a
large and significant difference. How can the footprints be identical when the genes clearly are
not?

**Because NES is a *rank* measure, and rank measures saturate.** NES asks *"are the targets
bunched at the up end of the list?"* — not *"how far up did they go?"* In both diseases miR-29's
targets are firmly at the up end, so both score about 2.2. Becker's collagen genes rose less far,
but they still rose further than most other genes, which is all a rank test can see.

```
DMD:  targets at ranks 1..1158 of 20823, having moved +4.06   -> NES 2.22
BMD:  targets at ranks 1..1158 of 20823, having moved +2.52   -> NES 2.17
                 ^ same position in the ranking, very different magnitude
```

**The magnitude difference is real and is visible only in the per-gene limma results (step 4) or
in the direct BMD-versus-DMD contrast — never in the vs-control NES values.**

### Required phrasing

> **🚫 Do not write:** "BMD's miR-29 footprint is weaker than DMD's", on the basis of the
> vs-control scores.
>
> **✅ Write:** BMD shows the same miR-29 footprint pattern as DMD, and the enrichment strength
> measured against healthy controls is indistinguishable between the two groups at matched
> sample size.

> **⚠ A sentence in an earlier draft of the logs was wrong and was corrected on 6 Nov 2026.** It
> claimed the direct BMD-versus-DMD comparison showed miR-29 targets "significantly less elevated
> in BMD". That rested on fgsea's padj of 4.7e−10. Step 9's patient-label shuffle gives
> **p = 0.239** for that same contrast. **The footprint difference is not significant and must
> not be claimed.** The per-*gene* collagen magnitudes (step 4) genuinely do differ and are safe
> to report, because limma permutes nothing. Keep the two claims apart.

**Outputs:** `results/tables/14_subsample_dmd.tsv`; `results/figures/14_subsample_dmd.png`

---

## Step 10b — The detection limit (`10b_detection_limit.R`, `10b2_calibrated_catch.R`)

**The question this step exists to answer.** Having found nothing, we owe the reader a
quantitative statement of *how small a signal we could have found*. "We found nothing" is only
informative if you also say what "nothing" means in units of effect size.

**It ran on different hardware**, which must be disclosed: Ubuntu 26.04, AMD EPYC 9R14, 16
physical cores, 30 GB RAM, 3 Nov 2026. The Mac (4 cores, 8 GB) needed ~15.5 CPU-hours and kept
sleeping. **Every relevant package version was identical** (limma 3.68.5, fgsea 1.38.0,
BiocParallel 1.46.0, data.table 1.18.6.1, ggplot2 4.0.3). Runtime **76.2 minutes** for the full
specification — 6 × 100 simulations, `nPermSimple = 10000`, no reductions.

### The method: plant a signal of known size and see if we catch it

1. **The ruler.** The DMD miR-29 footprint from step 6: **+0.2686 log2** (targets minus
   non-targets). This is "100%" — a full-strength, Duchenne-sized footprint.
2. **The plant.** Take a random gene set the same size as miR-29's (**1,158 genes**), matched on
   3′UTR-length profile. Add *f* × ruler to those genes **in the 11 BMD samples only**.
3. **The rerun.** Redo BMD-versus-control limma and fgsea, testing the planted set among all 128
   families.
4. **A "catch"** = the planted set comes out significant.
5. **The grid.** *f* ∈ {0, 10, 25, 50, 75, 100}%, **100 repeats each** = 600 simulations.

**At *f* = 0 nothing is planted, so the catch rate there is the false-alarm rate. It should be
about 5%.**

### ⚠ The specified catch rule failed its own sanity check

| | Value |
|---|---|
| Unplanted (*f* = 0) NES over 100 runs | mean **1.679**, SD 0.087, range 1.455–1.867 |
| Spec rule (NES > 0 **and** padj < 0.05) at *f* = 0 | **100 / 100 caught** |
| The *least* significant unplanted set | padj **2.1e−06** |

**Read that again.** With **nothing planted at all** — 1,158 random genes, unmodified — fgsea
called the set significant in **every single one of 100 runs**, and the weakest of those 100 had
an adjusted p-value of two in a million.

**The false-alarm rate is 1.00, not 0.05.** A detection curve that reads 100% at zero signal
carries no information whatsoever, and the script's own printout of "80% detection reached at
*f* = 0%" is an artefact of this. **Do not quote that line.**

Note also that the unplanted mean of **1.679** independently reproduces step 11's UTR-matched
null of 1.669. Two different scripts, two different methods, same answer: **a random UTR-matched
gene set scores NES ≈ 1.67.**

### The fix: a rule whose false-alarm rate is actually 5%

`10b2_calibrated_catch.R` re-scores the 600 **already-saved** runs — no limma or fgsea rerun —
using a bar set at the **95th percentile of the unplanted runs: NES > 1.812**.

| *f* (% of DMD footprint) | Planted log2 | Median NES | z vs unplanted | Catch, spec rule | **Catch, calibrated** |
|---|---|---|---|---|---|
| 0 | 0.000 | 1.67 | 0.00 | 1.00 | **0.05** (by construction) |
| 10 | 0.027 | 1.87 | 2.28 | 1.00 | **0.73** |
| 25 | 0.067 | 2.17 | 5.59 | 1.00 | **1.00** |
| 50 | 0.134 | 2.59 | 10.52 | 1.00 | 1.00 |
| 75 | 0.201 | 2.94 | 14.63 | 1.00 | 1.00 |
| 100 | 0.269 | 3.23 | 17.81 | 1.00 | 1.00 |

**Worked example of one row.** At *f* = 10 we add 0.1 × 0.2686 = **+0.027 log2** (about 1.9% on
the linear scale) to 1,158 genes in the 11 BMD samples. The median NES rises from 1.67 to 1.87.
**73 of the 100 runs** beat the 1.812 bar, so the catch rate is 0.73.

**80% detection is crossed at ≈ 14% of the DMD footprint (+0.037 log2)**, by straight-line
interpolation between the 10% and 25% grid points. *(The grid is coarse there; "about 14%" is
the honest precision.)*

**Important honesty note:** the 5% at *f* = 0 is **true by construction**, because the bar was
set using those same 100 runs. It is not an independent estimate of the false-alarm rate.

### How much of the "ruler" is just 3′UTR length?

(1,000 UTR-matched draws per contrast, same decile profile.)

| Contrast | miR-29 gap (targets − others) | UTR-matched random set gap | **Share explained by UTR length alone** |
|---|---|---|---|
| DMD − control | +0.2686 (the ruler) | +0.1124 | **42%** |
| BMD − control | +0.1923 (= 72% of the ruler) | +0.0973 | **51%** |

So BMD's miR-29-*specific* excess is 0.1923 − 0.0973 = **+0.095 log2**, which is **≈ 35% of the
ruler** on the planting scale.

### ⚠ The contradiction, and why it is the most interesting thing in the step

Line these three facts up:

1. The calibrated curve says we can catch a footprint of **14%** of DMD's, 80% of the time.
2. BMD's **real** miR-29 footprint is **≈ 35%** — comfortably above 14%, and at *f* = 35% the
   simulation catches it **100% of the time**.
3. Yet the honest patient-label shuffle (step 8) gives miR-29 in BMD **p = 0.084** — we do *not*
   detect it.

**How can the simulation say "easily detectable" while the real test says "not detected"?**

Because the planting is an **idealisation**, and every one of its assumptions runs in our favour:

- It adds **exactly the same shift to every single BMD patient.** Real biology varies enormously
  from patient to patient — our Becker patients range from 2.5 to 60.8 years old and from mild
  to severe.
- The planted genes are **random**, so they do not move together across patients. Real target
  sets are correlated blocs.
- The null **varies which genes are chosen**, not which patients are labelled. It is Null A
  wearing a disguise.

> ### 🔑 The conclusion, and it is a genuine finding in its own right
>
> **The binding limit at n = 11 is not the size of the gene-level signal. It is
> patient-to-patient variability** — and a planting simulation that applies a constant shift to
> random genes **cannot simulate that.** Step 10b, as specified, cannot measure the limit that
> actually matters for our null result. The calibrated curve is a **best-case floor**: the real
> detection limit is *at least* this large and probably much larger.

> **🚫 Do NOT write:** "our method detects footprints down to 14% of DMD's at n = 11."
>
> **✅ Write:** Under idealised planting, a footprint ≈ 14% of the DMD miR-29 footprint was
> detectable against a 3′UTR-matched null. The observed BMD miR-29 footprint exceeds this, yet
> does not survive patient-label permutation (p = 0.084). At n = 11 the detection limit is set by
> between-patient variability, which planting simulations of this kind do not capture.

**Outputs:** `results/tables/10b_detection_curve.tsv`, `10b_detection_runs.tsv` (all 600 runs),
`10b_log.txt`, `10b2_detection_curve_calibrated.tsv`, `10b2_ruler_utr_share.tsv`, `10b2_log.txt`;
`results/figures/10b_detection_curve.png`, `10b2_detection_curve_calibrated.png`

---

## Step 15 — Pig overlap — **DROPPED**

**The plan.** Reschke et al. 2026 (*J Cachexia Sarcopenia Muscle*) reported microRNA changes in a
pig model engineered to carry a BMD-like dystrophin deletion. Mapping their microRNAs to human
names and then to our families would have given an independent, cross-species check, with the
prediction: pig microRNA **up** → our BMD−Ctrl NES **< 0**.

**Why it was dropped.** **Table 1 of that paper could not be obtained.** The project rules state
that `pig_mirnas.tsv` may be built *only* from that table and from no other source, precisely so
that nobody is tempted to reconstruct directions from an abstract or a figure caption. It was
never built, the step never ran, and **no pig comparison appears anywhere in the paper.**

**Status changed from "blocked" to "dropped" at the 6 Nov computation freeze.** Mention it in
Limitations as an intended cross-species validation that could not be sourced — never as a
pending result.

---

# PART 6 — The three findings

The project did not answer its original question. It produced three other things, and the first
two are publishable on their own merits.

---

## 6.1 Finding 1 — fgsea's p-values are catastrophically overconfident on correlated target sets

**This is the strongest result in the project.** It is a methods finding, it is quantified, it is
reproduced four independent ways, and it has consequences for a very large published literature.

### The evidence, assembled

| Claim | fgsea (shuffle genes) | **Honest (shuffle patients)** | Gap |
|---|---|---|---|
| miR-29, DMD vs control | padj 1.5e−34 | **p = 0.052** | ~32 orders |
| miR-29, BMD vs control | padj 1.4e−31 | **p = 0.084** | ~30 orders |
| miR-29, BMD vs DMD | padj 4.7e−10 | **p = 0.239** | ~9 orders |
| Families differing BMD vs DMD | 7 of 128 | **0 of 128** | — |
| **Families altered vs controls in both groups** | **101 of 128** | **1 of 128** (6.4 expected by chance) | — |
| r between BMD and DMD footprints | 0.9123 | **p = 0.064** (null 95th pct 0.933) | — |

### Four independent demonstrations, each stronger than the last

1. **The direct comparison** (step 9): fgsea says 7 of 128 families differ between BMD and DMD;
   the patient shuffle says **0**.
2. **The classification table** (step 9b): fgsea says **101 of 128** families are altered in both
   diseases; the patient shuffle says **1** — below the 6.4 expected by chance.
3. **The UTR-matched null** (step 11): a random gene set matched only on 3′UTR length already
   scores NES ≈ **+1.6**. The implicit baseline of zero is simply wrong, and **72%** of miR-29's
   apparent enrichment is explained by UTR length alone.
4. **The clincher — step 10b's false-alarm test.** This is the one to lead with, because here
   **we know the ground truth is nothing.** We took 1,158 random genes, changed **absolutely
   nothing** about them, and ran the test 100 times. fgsea called them significant **100 times
   out of 100**, and the *least* significant of those had **padj = 2.1e−06**.

> **Demonstration 4 is not an inference, an estimate, or an argument. It is a controlled
> experiment with a known answer, and the method got it wrong every single time.**

### Why it happens

fgsea's default null asks *"would 1,158 **random genes** score this high?"* That question
silently assumes genes are independent. They are not. Collagen genes are co-regulated and move
as a bloc; so do the ribosomal genes, the immune genes, the mitochondrial genes. A **correlated
bloc** of 1,158 genes behaves nothing like a **scattered sample** of 1,158 genes, so the null is
far too narrow and p-values fall off a cliff.

Two things make our data an especially severe case, and both should be stated so the finding is
not overgeneralised:

- **Enormous gene sets.** miR-29's 1,158 targets are 5.6% of the measured transcriptome. Overlap
  between microRNA target lists is extensive.
- **A globally tilted ranked list.** ~11,000 of 20,823 genes shift in dystrophic muscle
  (step 4), so a *competitive* test finds almost any set enriched. 103 of 128 families came out
  "significant" in step 6, **99 of them in the same direction**.

### The honest scope of the claim

> **✅ Claim:** On microRNA target sets in a tissue with a globally shifted transcriptome,
> fgsea's gene-permutation adjusted p-values were overconfident by up to ~32 orders of magnitude
> relative to patient-label permutation, with a measured false-alarm rate of 1.00 where 0.05 was
> expected.
>
> **🚫 Do not claim:** "fgsea is broken", or that this invalidates GSEA generally. The tool is
> answering the question it was designed to answer — *"is this set special compared to other
> sets?"* — correctly. The failure is in **using that answer for a different question**: *"does
> disease status drive this?"* Frame it as a **mismatch between the null and the hypothesis**,
> not a bug. That framing is both more accurate and much harder to attack.

---

## 6.2 Finding 2 — the Becker result itself: a null, with a measured floor

**What we can say:**

1. **No microRNA family's target footprint differs between BMD and DMD** under patient-label
   permutation (0 of 128, all q > 0.95). The comparison excludes controls entirely, so this is
   not a shared-control artefact.
2. **miR-29 is not significantly altered** in either disease versus controls on the honest test
   (DMD p = 0.052, BMD p = 0.084) — though in DMD it is close enough to the threshold, and ranks
   **3rd of 128**, that "underpowered" is a far better description than "absent".
3. **BMD's footprint is not detectably a weaker copy of DMD's.** Against controls the two are
   indistinguishable at matched sample size (step 14: 53% of DMD-at-n=11 draws fall at or below
   BMD — a coin flip). The correlation between them, r = 0.91, is **inside the shuffled-label
   noise** (p = 0.064).
4. **The limiting factor is between-patient variability, not signal size** (step 10b).

**Why this is still worth publishing.** The field's default assumption is that Becker is
Duchenne-dialled-down. We tested that directly, with the correct null, and **could not confirm
it at n = 11** — and we can state quantitatively what we could and could not have seen. A null
result with a measured detection floor is a genuine contribution; a null result without one is
just a failed experiment.

**The honest framing of what the paper now is** — and this was the fallback the project plan
anticipated from the start:

> **"Can microRNA target footprints recover known microRNA biology in dystrophic muscle?
> Detection limits at n = 11."**

---

## 6.3 Finding 3 — GSE38417 is a duplicate re-deposit of GSE109178

While looking for more Becker samples, we found that **GSE38417 contains the same raw arrays as
GSE109178**, re-deposited under a second accession.

| Check | Result |
|---|---|
| Probe cells compared | **1,354,896** |
| Correlation of raw intensities | **r = 1.000000** |
| Cells differing | **0** |

Not "similar". **Bit-for-bit identical raw intensities.**

**Why this matters enough to report.** Merging the two accessions would have doubled the apparent
sample size using literally the same measurements. Every standard statistic would have treated
the duplicates as independent patients, halving standard errors and **manufacturing significance
from nothing**. Given that n = 11 is this field's binding constraint, the temptation to merge
Becker accessions is strong and the trap is easy to fall into.

**Practical output:** `scripts/check_duplicate_arrays.R`. **Run it before merging any accession.**

---

# PART 7 — What you may and may not write

A sentence bank. Use it directly.

## 7.1 🚫 Sentences that are now forbidden

| Forbidden | Why |
|---|---|
| ❌ "miR-29 is significantly altered in BMD / DMD" | Honest p = 0.084 / 0.052. Neither clears 0.05. |
| ❌ "BMD's footprint is a weaker copy of DMD's" | Step 14: coin flip. Step 8: r p = 0.064. |
| ❌ "BMD and DMD differ in microRNA family X" | Step 9: **0 of 128**, all q > 0.95. |
| ❌ "The direct comparison shows miR-29 significantly less elevated in BMD" | fgsea padj 4.7e−10, but honest p = **0.239**. Corrected 6 Nov. |
| ❌ "We detect footprints down to 14% of DMD's" | Best case only; the real limit is patient variability. |
| ❌ "r = 0.91 shows the footprints agree" | Shuffled labels reach 0.933 on 5% of runs. |
| ❌ Quoting any NES against a baseline of zero | The matched baseline is ≈ **+1.6** (step 11). |
| ❌ "80% detection reached at f = 0%" | A script artefact of the broken catch rule. |
| ❌ Anything citing the old preliminary GEO2R result | Superseded by step 4's limma tables. |
| ❌ "miR-29 ranks 3rd of 128 on the honest test" *near a BMD sentence* | True for **DMD only**; BMD is **45th**. |

## 7.2 ✅ Sentences that are supported

| Supported | Evidence |
|---|---|
| ✅ Collagen genes are up in both diseases, and significantly less so in BMD | Step 4; limma per-gene, no permutation involved |
| ✅ COL1A1 ≈ 16.7× in DMD, ≈ 5.7× in BMD | Step 4, log2 +4.058 / +2.517 |
| ✅ Consistent with Eisenberg 2007 in direction and order of magnitude | Step 4 vs the checked literature |
| ✅ The method recovers the known miR-29 direction in DMD | Step 6, NES +2.224 |
| ✅ The footprint is not merely scar genes | Step 10a: 96.7% of NES retained after removing all 275 matrisome genes |
| ✅ The footprint is not merely fibroblast content | Step 12: 79% retained, leads by 11 orders of magnitude after adjustment |
| ✅ The footprint exceeds a 3′UTR-matched baseline by 6.5 SD | Step 11; 0 of 1,000 matched sets reached it |
| ✅ miR-29's leading edge is the most fibrosis-dense of any family (6.5%) | Step 7 |
| ✅ The result is not driven by any single patient | Step 13: 17/17 runs, rank 1 throughout |
| ✅ BMD tissue is compositionally intermediate | Step 12: fibroblast score −0.791 / +0.040 / +0.253, ANOVA p = 1.2e−08 |
| ✅ No family differs between BMD and DMD | Step 9 |
| ✅ fgsea's false-alarm rate here is 1.00 where 0.05 was expected | Step 10b, 100/100 on unplanted sets |

## 7.3 Required hedging language

Per the project's standing rules, and because every one of these is an **inference from target
behaviour, not a measurement of a microRNA**:

- Use **"predicted"**, **"candidate"**, **"may contribute"**, **"consistent with"**.
- Never **"causes"**, **"proves"**, **"demonstrates"**, **"validated biomarker"**.
- State explicitly: **we infer microRNA activity from target footprints; we do not measure
  microRNAs.** The array cannot see them.
- With 11 BMD patients, **every result is hypothesis-generating.** Say so.

---

# PART 8 — Limitations, written out

Draft these; they are all real and a reviewer will find them anyway.

1. **Sample size.** 11 BMD, 6 controls. This is the binding constraint on everything, and no
   public dataset fixes it (Part 2.2).
2. **Six controls shared between two comparisons.** Quantified rather than hand-waved: shuffled
   labels produce r ≥ 0.933 on 5% of runs (step 8).
3. **Age confound.** BMD median 12.2 y (range 2.5–60.8); DMD median 3.5 y (range 0.9–8.0).
   Disease and age are entangled in the direct comparison. **Control ages are not reported in
   GEO at all** `[UNVERIFIED]`.
4. **No microRNA measurements.** Everything is inferred from predicted targets. A real miR-29
   measurement in BMD muscle would settle in one experiment what we could not settle in fourteen
   analysis steps.
5. **Predicted targets, not validated ones.** TargetScan conserved-site predictions contain false
   positives and miss real sites.
6. **Tissue composition.** Dystrophic biopsies contain fewer muscle cells and more fibroblasts.
   Step 12 adjusts for it, but the covariate also absorbs real disease signal, so neither the
   adjusted nor the unadjusted number is a clean estimate.
7. **3′UTR-length bias.** ~72% of miR-29's raw enrichment is attributable to UTR length alone
   (step 11). Any study quoting NES against zero on target sets is subject to this.
8. **Microarray, not RNA-seq.** Older platform, compressed dynamic range, and it cannot be merged
   with modern sequencing datasets.
9. **The detection-limit simulation is a best case.** It plants a constant shift on uncorrelated
   genes and cannot model between-patient variability — the very thing that limits us (step 10b).
10. **The muscle-expression threshold (mean RPMM ≥ 10) is ours**, not taken from any publication.
    It was never sensitivity-tested. Disclose it.
11. **No cross-species validation.** Step 15 was dropped because the source table could not be
    obtained (Part 5, step 15).
12. **500 label shuffles, not 1,000**, for runtime reasons. The empirical p floor is 1/501 =
    0.002, which is ample at the precision we need, but it is a deviation and is logged.

---

# PART 9 — Figures

14 figures exist in `results/figures/`. Per-figure reading notes are in `docs/figure_notes.md`.

## 9.1 Recommended for the paper

| Figure | File | What it shows |
|---|---|---|
| **1** | `06_mir29_enrichment_DMDvsCtrl.png` | The ranked walk. The method recovering known biology. |
| **2** | `07_nes_heatmap.png` | NES for top families across all three comparisons. Shows the vs-control tilt and its collapse in BMD−DMD. |
| **3** | `08_shuffle_r.png` | **The pivotal figure.** Real r = 0.9123 marked inside the shuffled-label distribution. |
| **4** | `10b2_detection_curve_calibrated.png` | Detection curve under both catch rules. The flat grey line at 100% *is* the finding. |
| **5** | `14_subsample_dmd.png` | BMD sits at the median of DMD-at-n=11. |
| Supp. | `02_qc_boxplot.png`, `02_qc_nuse.png`, `02_qc_within_group_cor.png` | QC, including the group-size artefact |
| Supp. | `09_mir29_bmd_vs_dmd_shuffle.png` | miR-29 in the direct comparison vs its shuffle null |
| Supp. | `12_fibroblast_score.png`, `13_mir29_leave_one_out.png` | Robustness |

## 9.2 Figure 3 is the one to get right

It is the figure that justifies the entire reframing of the paper, and it is the one a reviewer
will stare at longest. It must show:

- the histogram of 500 shuffled r values;
- the real r = 0.9123 marked with a vertical line;
- the 95th percentile at 0.933 marked, so the reader can see **the real value sits below it**.

**Caption to adapt:** *"Correlation between BMD and DMD microRNA footprints (r = 0.9123, red)
against 500 permutations of patient labels preserving group sizes (grey). Shuffled labels
produce r ≥ 0.933 on 5% of permutations. Empirical p = 0.064. The apparent agreement between
the two footprints is not distinguishable from what the shared-control design produces by
chance."*

**🚫 Avoid scatter plots** unless specifically asked for. If one is unavoidable, it must be
printed **next to** the shuffle histogram — a scatter plot of BMD NES against DMD NES looks like
overwhelming evidence and is, on its own, misleading for exactly the reason in Part 4.7.

---

# PART 10 — Reproducibility

## 10.1 Repository map

```
CLAUDE.md                     project specification + status block
scripts/                      00_setup.R ... 14_subsample_dmd.R, one per step
  09b_classify.R              step 9 classification (completed 6 Nov)
  10b2_calibrated_catch.R     step 10b recalibration
  check_duplicate_arrays.R    the GSE38417 duplicate detector
inputs/           Hamza's input files
inputs/auto/      input files auto-built to the same format (see docs/handoffs.md)
data/raw/         CEL files (gitignored, not redistributed)
data/processed/   RMA matrix, gene-level matrix (.rds)
results/tables/   36 .tsv/.txt outputs
results/figures/  14 .png
results/rds/      14 saved intermediate objects
docs/             methods_log.md, results_log.md, handoffs.md, figure_notes.md,
                  fact_check.md, NOTION_session_01.md, this document
```

Every script runs standalone: `Rscript scripts/NN_*.R`. Each reads from disk, writes to disk, and
sets `set.seed()` wherever randomness enters.

## 10.2 Seeds

| Steps | Seed |
|---|---|
| 1–4 | 1 |
| 6, 7, 10a, 12 (fgsea) | 42 |
| 8 | 20261101 |
| 9 | 20261102 |
| 11 | 20261103 |
| 10b | 20261104 |
| 10b2 | 20261105 |

## 10.3 Key parameters, with their justification

| Parameter | Value | Why |
|---|---|---|
| Probe collapse | highest mean expression | Project default |
| Low-expression filter | none | Project default; t spans −15.0 to +17.7, healthy |
| Target predictions | TargetScan 8.0 default (conserved), human | Project default |
| Muscle threshold | mean RPMM ≥ 10, limb muscle | **Ours.** Not from any publication — disclose |
| Family size | ≥ 15 array-present targets, no maximum | Project default |
| `nPermSimple` | **10,000**, `eps = 0` | At the default 1,000, **13 of 128 families returned NA NES**. Raised to remove them. |
| Label shuffles | **500** (steps 8, 9) | 1,000 would have taken 5.7 h for step 8 alone. Permitted fallback; logged |
| UTR-matched draws | 1,000 (step 11) | p floor 1/1001 |
| Detection repeats | 100 per grid point (steps 10b, 14) | — |

## 10.4 Runtimes, for the Methods

| Job | Time |
|---|---|
| Step 8, 500 shuffles, 7 cores | **206.7 min** |
| Step 9, 500 shuffles | included above |
| Step 10b, 600 simulations, 16 cores | **76.2 min** |
| Step 13, 17 leave-one-out runs | 7.3 min |
| Step 14, 100 subsample draws | 31.1 min |

## 10.5 Engineering lessons worth a sentence in a supplement

These cost real time and are genuinely useful to anyone reproducing this:

- **`fgsea` inside `bplapply` must be given `nproc = 1`.** Otherwise every parallel worker grabs
  every core, and 52% of CPU time is lost to context-switching. A 98-minute job was projecting
  to 6.8 hours before this was found.
- **Project runtimes from a parallel batch, never from a serial timing loop.** Step 8's serial
  loop projected ~98 minutes; the real parallel run took 206.7.
- **`setorder(dt, -abs(x))` is invalid in `data.table`** and crashed a step after two hours of
  work. Compute the column first.
- **`saveRDS` expensive objects *before* any optional block.** Step 9's classification block died
  on a missing file and would have taken the 500 shuffles with it had they not been saved first.
- **Check for duplicate accessions before merging anything** (Part 6.3).

---

# PART 11 — Open items

## 11.1 Decisions that are yours to make

1. **Which catch rule to headline in step 10b.** Both are in the repo and neither was
   silently switched. **Recommendation: headline the calibrated curve, and present the flat
   100% line of the spec rule as a *result* — it is the single cleanest demonstration of
   Finding 1 — rather than as a failed analysis.**
2. **How much weight to give Finding 1 versus Finding 2.** Finding 1 (fgsea overconfidence) is
   the stronger, more certain and more broadly useful result. Finding 2 (the Becker null) is
   what you set out to do. **Recommendation: lead with the biological question, and let the
   methods finding emerge as the reason the biological answer is null.** That is the honest
   narrative order and it reads better than the reverse.

## 11.2 Still unverified — check before submission

From `docs/fact_check.md` and the logs. **Do not cite any of these until checked against the
source:**

- **Control ages** — absent from GEO entirely `[UNVERIFIED]`.
- **The exact citation for the miRNA Tissue Atlas 2025 release** — the lineage (Ludwig et al.
  2016 NAR → Keller et al. 2022 NAR → 2025 update) is recorded, but the 2025 citation itself was
  never verified against the paper.
- **Eisenberg "5 BMD"** and **"GSE3307 = the same dataset"** — not found; the checked figure is
  **6 BMD**.
- **Xie 2024 gene counts** — not in the abstract.
- **The "13 genes" claim** — undefined, origin unknown.
- **Dadgar 2014 as the source of GSE109178** — unverified.
- **The pig "18 upregulated" claim and its directions** — unverified, and step 15 is dropped.
- **GSM2934804's severity** — disputed; immaterial to results either way (step 13).
- **GSE157668's contents** — never established.

## 11.3 Checked facts you may rely on

- **Eisenberg et al. 2007 (PNAS):** 6 BMD samples; only miR-146b and miR-221 passed their
  fold-change cutoff in BMD; microRNA–gene links were significant only in DMD and Miyoshi
  myopathy; in DMD, miR-29c ≈ 6-fold down, COL3A1 ≈ 6-fold up, COL1A2 ≈ 5-fold up. *(Units are
  not stated in the source, so write "about".)*

---

# PART 12 — The one-paragraph summary

> We asked whether the microRNA–fibrosis axis known to operate in Duchenne muscular dystrophy
> also operates in Becker, using 11 BMD, 17 DMD and 6 control muscle biopsies from GSE109178 and
> inferring microRNA activity from the collective behaviour of predicted target genes. The
> approach recovered the expected biology directionally: miR-29 targets were coordinately
> elevated in Duchenne muscle (NES +2.22), its leading edge was dominated by collagen genes, and
> the signal survived removal of all matrisome genes, adjustment for fibroblast content, a
> 3′UTR-length-matched null and leave-one-out resampling. **However, all of these used
> gene-permutation p-values.** When significance was instead assessed by permuting *patient
> labels* — the only null that tests whether disease status drives the footprint — the result did
> not survive: miR-29 reached p = 0.052 in Duchenne and p = 0.084 in Becker, no microRNA family
> differed between the two diseases (0 of 128), and the apparent correlation between their
> footprints (r = 0.91) fell inside the shuffled-label distribution (p = 0.064). The discrepancy
> is extreme and reproducible: on gene sets of this size and correlation structure,
> gene-permutation adjusted p-values were overconfident by up to ~32 orders of magnitude, and in
> a controlled test on **unmodified random gene sets** the method returned a significant result in
> **100 of 100 trials**. Simulation indicates that the limiting factor at n = 11 is
> between-patient variability rather than the magnitude of the gene-level signal. We therefore
> report a null result for the Becker comparison with a measured detection floor, and a
> methodological caution for the interpretation of microRNA target-enrichment analyses generally.

---

*Compiled 6 November 2026 from `docs/methods_log.md`, `docs/results_log.md`,
`docs/handoffs.md`, `docs/figure_notes.md` and the script outputs in `results/`, after the
computation freeze. Every figure quoted here traces to a script in this repository that ran to
completion.*
