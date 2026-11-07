# Fact-check register

Rule from CLAUDE.md §9: never invent a number, citation or claim. Anything unverified is
marked `[UNVERIFIED]` with the reason.

## Cleared in this session

| Claim | Status | Evidence |
|---|---|---|
| GSE109178 has 49 samples on GPL570, vastus lateralis | **CONFIRMED** | GEO sample records; `results/tables/01_sample_sheet.tsv`. All 49 are `GPL570`, all `vastus lateralis`. |
| 6 control / 11 BMD / 17 DMD / 15 LGMD | **CONFIRMED** | `scripts/01_download_labels.R` hard-asserts this and passed. LGMD splits 8 dysferlin (DYS) + 7 FKRP, matching "8 LGMD2B and 7 LGMD2I". |
| BMD ages run 2.5 to 60.8 years | **CONFIRMED** | min 2.5, max 60.8, median 12.2, all 11 known. |
| DMD patients are children | **CONFIRMED** | ages 0.9–8.0, median 3.5, all 17 known. |
| GSM2934804 severity | **RESOLVED as far as GEO goes** | GEO records it as `pathology: severe`, description `BMD-12-[6.4]-severe`, mutation `Deletion Exon 49`, age 6.4. The dispute is about whether that label is right, which GEO cannot settle. Kept as the step-13 robustness item. |
| miR-29 seed is `AGCACCA` | **CONFIRMED** | TargetScan 8.0 `miR_Family_Info.txt`: family `miR-29-3p`, seed `AGCACCA`, human members hsa-miR-29a-3p / 29b-3p / 29c-3p, conservation code 2 (highly conserved). |
| Eisenberg 2007: DMD COL3A1 ~6-fold up, COL1A2 ~5-fold up | **CONSISTENT with our data** | ours: COL3A1 logFC 3.26 = 9.6×, COL1A2 2.94 = 7.7×. Same direction, same order of magnitude. Different platform/normalisation, so "about" is the right word. |

## Citations checked 7 Nov 2026 (Hamza, with Claude)

Each one was read in its source text; the "Checked against" column says which text. Journal,
volume, pages and DOI for all seven were confirmed against Crossref or the article record on
7 Nov 2026. Quotes are verbatim. Page numbers are the journal's printed pages where given.

### The four §6.3 statistics citations (WRITEUP_PLAN §6.3, §10): all CLEARED

| Citation | Checked against | What it supports (verbatim) | ⚠ Caution for the write-up |
|---|---|---|---|
| **Goeman JJ, Bühlmann P.** Analyzing gene expression data in terms of gene sets: methodological issues. *Bioinformatics* 2007;**23**(8):980–987. doi:10.1093/bioinformatics/btm051. PMID 17303618 | Published PDF (OUP), supplied by Hamza | p.980 (abstract): gene-sampling P-values *"can be wildly anti-conservative"*. p.981: defines competitive vs self-contained tests, and gene-sampling vs subject-sampling methods. **p.984: *"the gene-sampling P-value does not say anything about biological replications of the experiment using different subjects."*** p.984: *"A significant P-value excludes random variation at the subject level as an explanation for the associations found, and therefore gives confidence that the same associations will be found for a new sample of subjects."* p.986: GSEA uses a gene-sampling statistic with a subject-sampling P-value, so it is *"similar to the hybrid method"*. | Their simulation (pp.984–985) tests the **2×2-table (hypergeometric-type) method, not GSEA or fgsea**. Cite them for the concept and the warning, not as a demonstration on our tool. |
| **Wu D, Smyth GK.** Camera: a competitive gene set test accounting for inter-gene correlation. *Nucleic Acids Res* 2012;**40**(17):e133. doi:10.1093/nar/gks461. PMID 22638577. PMC3458527 | Full published text (PMC) | *"Existing gene set tests that rely on gene permutation are shown here to be extremely sensitive to inter-gene correlation."* *"even quite modest correlations can dangerously inflate the apparent false discovery rate (FDR)"*. *"CAMERA is freely available as a function in the limma software package"*. | They also say sample permutation of a competitive statistic *"inevitably alters the hypothesis that is being tested"*, and call GSEA *"a hybrid test for which the null and alternative hypotheses are difficult to characterize"*. **Our label shuffle is that hybrid** (fgsea NES recomputed on shuffled labels). Say so ourselves in Methods §4.7, citing Goeman & Bühlmann and Wu & Smyth. |
| **Subramanian A, Tamayo P, Mootha VK, et al.** Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *PNAS* 2005;**102**(43):15545–15550. doi:10.1073/pnas.0506580102. PMID 16199517. PMC1239896 | Published PDF, hosted by the Broad Institute (gsea-msigdb.org) | *"the permutation of class labels preserves gene-gene correlations and, thus, provides a more biologically reasonable assessment of significance than would be obtained by permuting genes."* On gene permutation: *"This approach is not strictly accurate: because it ignores gene-gene correlations, it will overestimate the significance levels and may lead to false positives."* | The same passage adds that gene permutation *"can be useful for hypothesis generation"*. Fair to acknowledge; it matches our "hypothesis-generating" framing. |
| **Efron B, Tibshirani R.** On testing the significance of sets of genes. *Ann Appl Stat* 2007;**1**(1):107–129. doi:10.1214/07-AOAS101 | arXiv math/0610667v2, the IMS reprint (same text as published, different pagination) | On "row randomization" (= gene sampling): *"Its defect concerns correlation between the genes"*; *"row randomization destroys gene-wide correlations"*. *"Permutation methods deal nicely with the correlation problem."* Proposes restandardization and the maxmean statistic (R package GSA). | They also show the **reverse failure**: in a simulation where *"the first 10 genes in every gene-set are 2.5 units higher in the treatment group"* (half of all genes shifted, nothing special about any one set), sample-label permutation alone makes *"all of the gene-sets look to be significant"*. That is why they propose restandardization, and it closely matches our globally tilted list (~half the transcriptome moving). §6.3 should say neither null is perfect. No journal page numbers from the reprint, so no pin cites. |

### Age-dependence of muscle ECM/collagen expression (WRITEUP_PLAN §5.2, §7.3): PARTIAL

**No source was found showing that human muscle collagen gene expression is higher at ~3.5 y than
at ~12 y** (our DMD vs BMD medians). The best available evidence supports a broader statement,
*"ECM gene expression in muscle is higher early in life and falls with maturation"*. Use #1 + #2.

| Citation | Checked against | Supports (verbatim) | Does not support |
|---|---|---|---|
| **1. Schaiter A, Hentschel A, Kleefeld F, et al.** (senior authors Roos A, Schänzer A). Molecular composition of skeletal muscle in infants and adults: a comparative proteomic and transcriptomic study. *Sci Rep* 2024;**14**:22965. doi:10.1038/s41598-024-74913-4. PMID 39362957. PMC11450201 | Full text (Europe PMC) | **Human.** RNA-seq on 4 infants (7–28 months) vs 7 adults (30–56 y): *"transcriptomic analyzes showed down-regulation of GO cellular component terms 'actin cytoskeleton' and 'extracellular matrix organization' in adults (Supplemental Fig. 3)"*. The ECM protein tenascin-X (TNXB) was also *"significantly down-regulated in our adult muscle samples"*. | **Never names a collagen gene.** Very small n. Infant biopsies all quadriceps; adult biopsies mixed (quadriceps 6, deltoid 2, gastrocnemius 1, intercostal 1, of the full 10). All biopsies were diagnostic, for suspected neuromuscular or cardiac disease, and judged normal. |
| **2. Listrat A, Jousse C, Tournayre J, Boby C, Lee K, Béchet D, Wang H, Goh KL.** Transcriptional profiling of male mouse muscle across aging stages: a gene ontology analysis of the muscle matreotype. *BMC Genomics* 2026;**27**:380. doi:10.1186/s12864-026-12677-z. PMID 41803688. PMC13085706 | Full text (Europe PMC) | **Mouse** gastrocnemius at 2, 11 and 25 months, 18 male C57BL/6: *"The expression of all collagen isoforms (a total of 16 of the 43 described in the Mus musculus matrisome) decreased during both maturation and aging, but no collagens showed a specific decrease during aging (11 to 25 Mo)"*. The only direct collagen statement found. | Mouse only. Starts at 2 months, a young-adult mouse, not a juvenile. |
| **3. Simon A, Djeddi S, Bournon P, Reiss D, Thompson J, Laporte J.** Transcriptomic characterization of postnatal muscle maturation. *Dis Model Mech* 2025;**18**(2):dmm052098. doi:10.1242/dmm.052098. PMID 39945189. PMC11911633 | Full text (Europe PMC) | **Mouse**, E18.5 / 2 wk / 7 wk, n = 4 per age (tibialis anterior at 2 and 7 wk). For 2 → 7 wk (juvenile to adult): *"During late maturation, there was mainly downregulation of genes associated with cell cycle- and extracellular matrix-related processes."* | Mouse only. Gene-category (GO) level; no named collagens. Optional third citation. |
| Kovanen V, Suominen H 1989, *Eur J Appl Physiol*: rat; collagen-synthesis enzyme activity in muscle reportedly highest in the youngest animals, falling by ~2 months | `[UNVERIFIED]`: search snippets only, no DOI confirmed | Awaiting Hamza's PDF | — |
| ~~Su et al. 2015, *Skeletal Muscle*: "ECM genes up with adult aging"~~ | Full text (PMC4600214) | **Do not cite for this.** The line is Su describing Zahn et al. 2006, not Su's own result. Claude misattributed it in chat on 7 Nov; corrected here. | — |

**Suggested caveat wording (from #1 and #2):**
> Age could not be separated from disease severity: DMD patients were younger (median 3.5 y)
> than BMD patients (median 12.2 y), and control ages are not reported. Extracellular-matrix gene
> expression in skeletal muscle is higher early in life and declines with maturation: all
> collagen genes fall in maturing mouse muscle [Listrat 2026], and ECM-organisation genes are
> lower in adult than in infant human muscle [Schaiter 2024]. Part of the smaller collagen
> increase in BMD may therefore reflect age rather than milder disease.

### Input-source citations (WRITEUP_PLAN §4.5, §4.9, §10): both CLEARED, one CORRECTED

| Citation | Checked against | Result |
|---|---|---|
| **miRNA Tissue Atlas 2025** (muscle-expression filter) → **Rishik S, Hirsch P, Grandke F, Fehlmann T, Keller A.** miRNATissueAtlas 2025: an update to the uniformly processed and annotated human and mouse non-coding RNA tissue atlas. *Nucleic Acids Res* 2025;**53**(D1):D129–D137. doi:10.1093/nar/gkae1036 | The atlas home page (`web.ccb.uni-saarland.de/mirnatissueatlas_2025`) and Crossref | **CLEARED.** The site says: *"If you use miRNATissueAtlas please cite the latest version below: Rishik S, Hirsch P, Grandke F, Fehlmann T, Keller A. miRNATissueAtlas 2025 …"*. The other two DOIs on the page are the earlier versions: v1 = Ludwig et al. 2016, *NAR* 44(8):3865–3877 (gkw116); v2 = Keller et al. 2022, *NAR* 50(D1):D211–D221 (gkab808). Cite gkae1036 only. The site gives the year as 2024 (online); the issue is Jan 2025, so cite 2025. **Download check (7 Nov): REPRODUCED EXACTLY.** Our file comes from the downloads page's **"miRNA avg std human"** link, which serves `all_human_miRNA_rpmm_harmonized_meta_data.parquet` (12,940,357 bytes). Running the unchanged `scripts/build_inputs_muscle.py` on today's download gives *"2656 miRNAs; >=10 RPMM: 353"*, as logged on 1 Nov. Against `inputs/auto/mirna_families.tsv`, **2,064 of 2,064 families have identical `max_muscle_rpmm`, `n_members_in_atlas` and `expressed_in_muscle`** (253 muscle-expressed; miR-29-3p 1,894.69). *Correction:* an earlier note here said the site now serves only `.h5ad`. That was wrong: the `.h5ad` (1.0 GB) is the per-sample matrix behind a different link, "miRNA expression human". **Methods:** name the "miRNA avg std human" download. |
| **Rubenstein AB, Smith GR, Raue U, et al.** Single-cell transcriptional profiles in human skeletal muscle. *Sci Rep* 2020;**10**:229. doi:10.1038/s41598-019-57110-6. **PMID 31937892**. PMC6959232 | MSigDB gene-set pages for both `RUBENSTEIN_SKELETAL_MUSCLE_FAP_CELLS` (M41745) and `…_FBN1_FAP_CELLS` (M41750), plus the PubMed record via Europe PMC | **CORRECTED.** The PMID we had recorded, **32042044, is wrong**: it is Byun et al. 2020, *Nat Commun* 11:807, a liver-autophagy paper. MSigDB lists PMID 31937892 for both sets, so the earlier note "taken from MSigDB metadata" was also wrong; where 32042044 came from is unknown. Content fits our use: the abstract reports scRNA-seq of mononuclear cells *"from human vastus lateralis"* (the same muscle as GSE109178) and *"two fibro-adipogenic progenitor (FAP) cell subtypes"*. MSigDB's exact source is Supplementary Table S1; dataset GSE130646. Fixed on 7 Nov in `scripts/build_inputs.R`, the `source` column of `inputs/auto/fibroblast_markers.tsv` and `docs/WRITEUP_PLAN.md`. **No rerun needed:** no script reads the `source` column. |

## Still unverified — do not put in the paper yet

| Item | Why it is unverified |
|---|---|
| Control ages and sexes in GSE109178 | `[UNVERIFIED]` — **not reported in GEO at all.** All 6 control records have `age: NA` and no sex. CLAUDE.md §2 asked us to record them; they do not exist in the source. The age confound therefore **cannot be resolved** and must be stated as a limitation. See §Age confound below. |
| Reschke et al. 2026 pig microRNA list and directions | `[UNVERIFIED]` — Table 1 not obtained. Step 15 BLOCKED. |
| Eisenberg "5 BMD", "GSE3307 = same dataset" | still not found. Do not use. |
| Xie 2024 gene counts; the "13 genes" claim; Dadgar 2014 as the source of GSE109178; the pig "18 upregulated"; GSE157668 contents | unchanged from CLAUDE.md §11. Do not use. |

## Age confound (must appear in Limitations)

- Control ages: **unknown**, not in GEO.
- BMD median 12.2 y (range 2.5–60.8). DMD median 3.5 y (range 0.9–8.0).
- So the **direct BMD-vs-DMD comparison confounds disease severity with age**: the BMD group is
  roughly a decade older. Any BMD-vs-DMD difference may be partly a child-vs-older-patient
  difference. We cannot adjust for this, because adjusting would also remove the disease
  signal and because control ages are missing entirely.
- State it; do not try to fix it (CLAUDE.md §2).

---

# GEO duplicate-accession discovery (2026-11-03)

**GSE38417 is NOT an independent dataset. It is a re-deposit of arrays already in GSE109178.**

Verified at the raw-data level, not by metadata:

| Check | Result |
|---|---|
| `GSM941830` (GSE38417 control) vs `GSM2934813` (GSE109178 control) | |
| probe cells compared | 1,354,896 |
| intensities identical | **TRUE** |
| Pearson correlation | **1.00000000** |
| cells differing | **0** |

Both are HG-U133_Plus_2. Both carry the same lab-internal array IDs
(`MIDO-CTRL_Pos-1UP` … `Pos-6UP`, `MIDO-DMD_Pos-11UP` …). The biopsies are Hoffman's (CNMC),
but **the two accessions were deposited by different groups** (checked against GEO on 7 Nov):
GSE38417 on 1 Jun 2012 by Susan G. Dorsey and Christopher W. Ward, University of Maryland
Baltimore (its summary credits Hoffman's histology at CNMC; no PubMed link); GSE109178 on
14 Jan 2018 from CNMC (contributor Eric P. Hoffman). The earlier note "Same lab (Hoffman, CNMC)"
was wrong. GSM941830's CEL file is `GSM941830_MIDO-CTRL_Pos-1UP.CEL`, the same internal ID as
`GSM2934813_MIDO-CTRL_Pos-1UP.CEL` (confirmed on GEO). GSE38417 has 22 samples.
Only the GEO accession and the CEL file encoding (text vs binary) differ.

**Consequence if this had been missed:** merging the two accessions would have double-counted
6 controls and ~16 DMD samples. Standard errors would have fallen by roughly sqrt(2) with no
new information, and miR-29's DMD-vs-control patient-shuffle p of 0.052 would probably have
crossed 0.05 — a "significant" result manufactured entirely by duplication.
`[SPECULATION — not computed; do not cite. WRITEUP_PLAN §9 forbids this sentence.]`

**Rule for this project:** before merging ANY GEO accession, compare raw CEL intensities
against the samples already held. Matching `Sample_title` / filename stems are a red flag;
identical intensities are proof. Script: `scripts/check_duplicate_arrays.R`.

## GSE13608 — genuine, but rejected on design grounds

Its 6 "Normal adult" samples ARE independent (best correlation with our controls 0.849 and
0.786, with >1.3 million cells differing — ordinary between-person variation).

Rejected anyway, because:
- **biopsy site recorded as "unknown"** for all 6; every GSE109178 sample is vastus lateralis;
- ages **26, 43, 43, 43–46, 61, 85**, versus DMD patients aged 0.9–8;
- one sample (`GSM343080`) is a **pooled M/F sample**, not one person;
- it would buy only 6 -> 12 controls (BMD-vs-control precision 0.508 -> 0.417, ~18%), which
  cannot rescue p = 0.052, in exchange for a muscle-site confound and an age confound.

## Other BMD sources in GEO

`GSE291383` (RNA-seq, GPL21697) holds **6 BMD + 12 DMD and NO controls**. It is the only other
human BMD muscle dataset found. Usable only as an independent check of the BMD-vs-DMD null,
and cross-platform (RNA-seq vs microarray).

**Conclusion: 11 BMD is the binding constraint on this project, and no public dataset fixes it.**
