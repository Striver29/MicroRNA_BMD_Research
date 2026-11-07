# Pre-submission verification: every number and citation (7 Nov 2026)

**Asked by Hamza:** "confirm all the data and citations to be true." Done with Claude, 7 Nov.
Nothing was recomputed; computation is still frozen. Every check reads the frozen outputs in
`results/tables/`, `inputs/auto/` and the step logs, or a primary source (GEO, Crossref,
PubMed/Europe PMC, the paper's own text).

**Bottom line**

| Check | Result |
|---|---|
| Numbers in `docs/WRITEUP_PLAN.md` | **95 checks against the result tables.** **4 numbers were wrong** (#2–#5 below), and reading the text turned up **3 wrong statements** (#1, #6, #7). All 7 are fixed; after the fixes, **95/95 pass**. Re-run any time: `python3 scripts/check_writeup_numbers.py` |
| Sample labels vs GEO | **49 of 49 match** GEO's own records on group, title, age, severity, platform GPL570, vastus lateralis and CEL file name. 0 mismatches. This replaces the "spot-check 2–3 GSMs" item. |
| Citations | **20 verified** against Crossref or PubMed, and against the paper's text wherever it is used for a claim. **One PMID was wrong** (Rubenstein; fixed earlier today). **One claim was wrong** (GSE38417 "same lab"; fixed). **Two are still pending** source text: Eisenberg 2007 plus its correction, and Kovanen 1989. |

---

## 1. Errors found and fixed

| # | Where | Was | Correct | Evidence |
|---|---|---|---|---|
| 1 | Abstract; §5.1 | Leave-one-out (17/17) listed as support for the **DMD** vs control miR-29 result | Step 13 reran **BMD vs control and BMD vs DMD only**. It never retested DMD vs control. Reworded so leave-one-out supports the BMD result. | `13_leave_one_out.tsv` contrasts; `scripts/13_leave_one_out.R` lines 25–27 |
| 2 | §5.1, `results_log.md`, `FULL_RESEARCH_DOCUMENT.md` | Leading edge "highest of any family (next 4.1%)" | Highest **of the 25 top families** (leading edges were extracted only for these); **next is 5.3%** (miR-381-3p); 4.1% is fourth | `07_leading_edge_summary.tsv` |
| 3 | same three files | Fibroblast adjustment: other top families keep "~55–60%" | The other top-10 families keep **48–74% (median 64%)**; miR-29's 79% is still the highest | `12_fibro_adjusted_fgsea.tsv` |
| 4 | same three files | After adjustment miR-29 leads "by eleven orders of magnitude, 1.7e−16 vs 3.1e−05 next" | Next family is **miR-144-3p at 5.5e−08**, about eight orders of magnitude. 3.1e−05 is miR-30-5p, the top family *before* adjustment | `12_fibro_adjusted_fgsea.tsv` |
| 5 | §5.2 table (all three files) | LOX BMD−Ctrl adj.P 2.4e−07 | **2.3e−07** (exact 2.349e−07; the old value was rounded twice, via 2.35e−07) | `limma_BMDvsCtrl.tsv` |
| 6 | §5.5; `fact_check.md` | GSE38417 and GSE109178 "from the same group (Hoffman, CNMC)" | **Different depositors:** GSE38417 by Dorsey & Ward, Univ. Maryland Baltimore, 1 Jun 2012 (summary credits Hoffman's histology; no PubMed link); GSE109178 from CNMC (Hoffman), 14 Jan 2018 | GEO series records |
| 7 | §5.4 | "at f = 35% the simulation catches it 100% of the time" | No 35% point was run. The calibrated rule already catches **100% at f = 25%**. | `10b2_detection_curve_calibrated.tsv` |

**Also added to Methods (§4.7, §4.10), as omissions rather than errors:**
- The shuffle p uses **only shuffles where fgsea returned a finite NES** for that family. miR-29
  keeps **497 of 500 in step 8** (p = 26/498 for DMD, 42/498 for BMD) and **496 in step 9**.
  106 of 128 families (DMD vs ctrl) and 110 (BMD vs ctrl) keep all 500.
  Source: `08_shuffle_baseline.R` line 88; the denominators are recoverable exactly from the
  p-values.
- **Seeds** corrected against `set.seed()` in each script: 1 (steps 2–5, `build_inputs.R`; step 1
  has none); 42 (6, 7, 10a, 12); 20261101 (8); 20261102 (9, 9b); 20261103 (11); 20261104 (10b);
  20261105 (10b2, **13**); **20261106 (14)**.
- §11 schedule dates corrected (Fri **6**, Sat **7** Nov).

**⚠ Not fixed, for Aziz:** the `CLAUDE.md` status block (line 45) still lists "17/17
leave-one-out runs" among the DMD-recovery evidence (error #1). It is the project instruction
file, so it was left for Aziz to edit.

## 2. Numbers verified correct as written

All checked against the frozen tables; full list in `scripts/check_writeup_numbers.py` output.
- **Samples:** 49 total, 34 used, 6/11/17, 15 LGMD (8 dysferlin + 7 FKRP); BMD median 12.2
  (2.5–60.8), 6 mild/3 moderate/2 severe; DMD median 3.5 (0.9–8.0), 4 moderate/13 severe;
  control ages all NA.
- **QC:** group medians 0.9529/0.9704/0.9803 overall and 0.9839/0.9757/0.9855 within group;
  |RLE| ≤ 0.057; NUSE 0.981–1.039; 0 RLE/NUSE flags; 3 controls flagged on the spec metric;
  GSM2934803 0.9676 and GSM2934813 0.9784.
- **Probes and genes:** 54,675 → 44,647 → −1,524 → 43,123 → 20,823. DE counts for all three
  contrasts. Collagen table (except LOX, fixed) and fold changes 16.7×/9.6×/7.7×/5.7×.
- **Inputs:** 2,064 families → 253 muscle → 128 (15–1,446, median 371); 115,308 / 219 / 13,034;
  19,429 UTRs; 275/54/2 fibrosis genes; 372/326/316 markers; miR-29 seed, members, class 2,
  1,158 targets, 1,895 RPMM.
- **Results:** ruler 0.2686 (20.5%); targets 0.304 vs 0.036; NES 2.224 (2nd by NES, 3rd by padj),
  103/99 families; ECM 96.7%; fibroblast 79%, 103 → 51 families; UTR z 6.56, null means
  1.592/1.669, 72%; shuffle p 0.052/0.084/0.239, q 0.951, ranks 3rd/45th, 1 and 10 families,
  7 → 0 of 128, 101/27 vs 1/127 (miR-196-5p); r and slope nulls (from the step-8 log); 10b f = 0
  100/100, weakest padj 2.1e−06, NES 1.679 ± 0.087 (1.455–1.867), bar 1.812, calibrated
  0.05/0.73/1/1/1/1, 14% = +0.037; excess 0.095 = 35%, UTR shares 42%/51%; leave-one-out
  ranges, rank 1, GSM2934804 2.156, top-20 overlap 18–20; subsample median 2.167
  (2.074–2.258), 53.0%; fibroblast scores −0.791/+0.040/+0.253, ANOVA p 1.2e−08; GSE13608
  precision 0.508 → 0.417.

**Confirmed only from logs or saved objects, not recomputable here** (no R or CEL files in this
container):
- eBayes **prior df 4.23**: the only plausible scalar in `results/rds/limma.rds` is 4.226.
- **Atlas muscle filter, re-run from a fresh download:** the "miRNA avg std human" parquet plus
  the unchanged `build_inputs_muscle.py` gives **353 of 2,656** microRNAs, and 2,064/2,064
  families identical to `inputs/auto/mirna_families.tsv` (253 muscle-expressed).
- Runtimes 206.7 / 76.2 / 7.3 / 31.1 min: the step logs.
- **Logged only, with no saved output to re-check:**
  - 11,043 genes with more than one probeset
  - 13 of 128 NA NES at fgsea's default permutation count
  - the GSE38417 bit-identity (1,354,896 cells, r = 1). The matching internal ID
    `MIDO-CTRL_Pos-1UP` *was* confirmed on GEO.
  - Bioconductor 3.23: not named in `sessionInfo.txt`; all listed package versions match.
  - the step-10b hardware

## 3. GEO facts

| Claim | Status |
|---|---|
| GSE109178: 49 samples, GPL570, all vastus lateralis, 6/11/17/7/8 | ✅ GEO series and all 49 sample records |
| **Source paper = Dadgar et al. 2014** (was `[UNVERIFIED]` in CLAUDE.md §11) | ✅ GEO lists PMID 25313409 (and 32592467, Dang et al. 2020). The paper describes exactly this set: *"49 patient muscle biopsies analyzed on HG-U133 plus 2.0 microarrays… 6 normal volunteers, 17 DMD, 11 BMD, 7 LGMD2I… and 8 LGMD2B"* |
| GSE38417 shares arrays with GSE109178 | ✅ IDs match on GEO; 22 samples; **depositor corrected (error #6)** |
| GSE13608: 6 normal adults, site unknown, ages 26–85, one pooled sample | ✅ ages 26, 43, 43, 43–46, 61, 85; site "unknown" for all 6; GSM343080 is an M/F pool |
| GSE291383: 6 BMD + 12 DMD, no controls, RNA-seq, GPL21697 | ✅ |
| "Eisenberg 5 BMD / GSE3307" confusion | Likely origin: Dadgar 2014's text says GSE3307 is the **Bakay et al. 2006** set, with *"BMD (n = 5)"*. Still do not use. |

*Observation, not a claim:* GEO's pathology note for **GSM2934805** (a BMD sample, age 5.8,
"severe") reads *"Biopsy not great, freeze artefact, mid stage DUchenne?"*. It is not in the
leave-one-out list of concerns, and step 13 shows no single BMD sample moves the result.

## 4. Citations

**Verified** (bibliographic details against Crossref or PubMed; claims against the text where we
use them; quotes in `docs/fact_check.md`):

| Use in paper | Citation |
|---|---|
| Data source | Dadgar S, Wang Z, Johnston H, et al., Hoffman EP. Asynchronous remodeling is a driver of failed regeneration in Duchenne muscular dystrophy. *J Cell Biol* 2014;207(1):139–158. doi:10.1083/jcb.201402079. PMID 25313409 |
| §6.3 | Goeman JJ, Bühlmann P. *Bioinformatics* 2007;23(8):980–987. doi:10.1093/bioinformatics/btm051 (published PDF, page-cited) |
| §6.3 | Wu D, Smyth GK. *Nucleic Acids Res* 2012;40(17):e133. doi:10.1093/nar/gks461 |
| §6.3 | Subramanian A, Tamayo P, Mootha VK, et al. *PNAS* 2005;102(43):15545–15550. doi:10.1073/pnas.0506580102 |
| §6.3 | Efron B, Tibshirani R. *Ann Appl Stat* 2007;1(1):107–129. doi:10.1214/07-AOAS101 |
| Age caveat | Schaiter A, et al. *Sci Rep* 2024;14:22965. doi:10.1038/s41598-024-74913-4 · Listrat A, et al. *BMC Genomics* 2026;27:380. doi:10.1186/s12864-026-12677-z · (optional) Simon A, et al. *Dis Model Mech* 2025;18(2):dmm052098 |
| Muscle filter | Rishik S, Hirsch P, Grandke F, Fehlmann T, Keller A. *Nucleic Acids Res* 2025;53(D1):D129–D137. doi:10.1093/nar/gkae1036 |
| Fibroblast markers | Rubenstein AB, et al. *Sci Rep* 2020;10:229. doi:10.1038/s41598-019-57110-6. **PMID 31937892** (corrected) |
| Matrisome set | Naba A, Clauser KR, Hoersch S, Liu H, Carr SA, Hynes RO. *Mol Cell Proteomics* 2012;11(4):M111.014647. doi:10.1074/mcp.M111.014647. PMID 22159717 |
| Hallmark TGF-β set | Liberzon A, Birger C, Thorvaldsdóttir H, Ghandi M, Mesirov JP, Tamayo P. *Cell Syst* 2015;1(6):417–425. doi:10.1016/j.cels.2015.12.004 |
| TargetScan 8.0 (the site asks for either) | McGeary SE, et al., Bartel DP. *Science* 2019;366(6472):eaav1741. doi:10.1126/science.aav1741 · Agarwal V, Bell GW, Nam JW, Bartel DP. *eLife* 2015;4:e05005. doi:10.7554/eLife.05005 |
| RMA | Irizarry RA, Hobbs B, Collin F, Beazer-Barclay YD, Antonellis KJ, Scherf U, Speed TP. *Biostatistics* 2003;4(2):249–264. doi:10.1093/biostatistics/4.2.249 |
| oligo | Carvalho BS, Irizarry RA. *Bioinformatics* 2010;26(19):2363–2367. doi:10.1093/bioinformatics/btq431 |
| limma | Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK. *Nucleic Acids Res* 2015;43(7):e47. doi:10.1093/nar/gkv007 · Smyth GK. *Stat Appl Genet Mol Biol* 2004;3:Article 3. doi:10.2202/1544-6115.1027 |
| fgsea | Korotkevich G, Sukhov V, Budin N, Shpak B, Artyomov MN, Sergushichev A. Fast gene set enrichment analysis. *bioRxiv* 060012. doi:10.1101/060012 (preprint) |

**Pending, needing source text that only Hamza's university login can reach:**
- **Eisenberg I, et al., Kunkel LM.** Distinctive patterns of microRNA expression in primary muscular
  disorders. *PNAS* 2007;104(43):17016–17021. doi:10.1073/pnas.0708115104. PMID 17942673.
  **It has a published correction, *PNAS* 2008;105(1):399 (PMC2224226)**, which no project doc
  mentions. The four Eisenberg facts (6 BMD; miR-146b/miR-221; DMD/Miyoshi linkage; miR-29c and
  collagen folds) have **no recorded source-text check**. Read both PDFs before citing any of them.
- **Kovanen & Suominen 1989**, *Eur J Appl Physiol*: optional age-caveat support; snippets only.

**Do not cite:** Xie 2024; Reschke et al. 2026 (step 15 dropped); Su et al. 2015 for "ECM up
with age" (that line describes Zahn et al. 2006).
