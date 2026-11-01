# Handoffs: auto-built inputs

Per CLAUDE.md §4: every file Hamza owns that was not in `inputs/` was rebuilt to the same
format in `inputs/auto/` so the pipeline could keep moving. Loading rule in every script is
`inputs/<file>` if it exists, else `inputs/auto/<file>`, and the script prints which it used.

**When Hamza's version arrives:** diff it against the auto version, report the differences to
Aziz, switch to Hamza's version, rerun the affected steps unless results are frozen (§7).

Built by `scripts/build_inputs.R` (+ `scripts/build_inputs_muscle.py`) on 2026-11-01.

| File | Status | Source and version | Method | Rows |
|---|---|---|---|---|
| `sample_sheet.csv` | auto | GEO GSE109178 sample records, retrieved 2026-11-01 | `scripts/01_download_labels.R` parses `!Sample_description`, `!Sample_characteristics_ch1` | 49 (34 included) |
| `mirna_families.tsv` | auto | TargetScanHuman 8.0 `miR_Family_Info.txt` (file dated 2018-04-06), filtered to species 9606 | seed = `Seed+m8` column (positions 2–8); muscle flag from the miRNA atlas below | 2,064 families |
| `mirna_targets.tsv` | auto | TargetScanHuman 8.0 `Predicted_Targets_Info.default_predictions.txt`, species 9606 | conserved-site default predictions, as CLAUDE.md §10 specifies | 115,308 rows, 219 families, 13,034 genes |
| `utr_lengths.tsv` | auto | TargetScanHuman 8.0 `UTR_Sequences.txt`, species 9606 rows only | alignment gap characters `-` stripped, then `nchar()`; longest transcript per gene kept | 19,429 genes |
| `fibrosis_sets.tsv` | auto | MSigDB 2026.1.Hs via `msigdbr` 26.1.1 | `NABA_CORE_MATRISOME` (C2:CGP, Naba et al. 2012, PMID 22159717) = 275 genes; `HALLMARK_TGF_BETA_SIGNALING` = 54 genes; overlap 2 | 329 rows |
| `fibroblast_markers.tsv` | auto | MSigDB 2026.1.Hs C8: `RUBENSTEIN_SKELETAL_MUSCLE_FAP_CELLS` + `RUBENSTEIN_SKELETAL_MUSCLE_FBN1_FAP_CELLS` (Rubenstein et al. 2020) | union of the two FAP (fibro-adipogenic progenitor) signatures; `in_fibrosis_set` TRUE for the 46 that are also in a fibrosis set | 372 markers, 326 usable for the covariate |
| `pig_mirnas.tsv` | **NOT BUILT** | Reschke et al. 2026, J Cachexia Sarcopenia Muscle, Table 1 | Table not obtained. CLAUDE.md §4 forbids building this file from anything else. **Step 15 is BLOCKED.** | — |

## Muscle-expression source (CLAUDE.md §4 asked for a proposal)

**Proposed and used:** the miRNA Tissue Atlas, 2025 release, Center for Bioinformatics,
Saarland University (`web.ccb.uni-saarland.de/mirnatissueatlas_2025`), human harmonised
RPMM table `all_human_miRNA_rpmm_harmonized_meta_data.parquet`, retrieved 2026-11-01.
Lineage of this resource: Ludwig et al. 2016 NAR (`10.1093/nar/gkw116`) → Keller et al. 2022
NAR (`10.1093/nar/gkab808`) → 2025 update (`10.1093/nar/gkae1036`).
**Aziz: the exact citation for the 2025 release still needs checking against the paper — see
`docs/fact_check.md`.**

- Skeletal muscle = rows with `Biotype = "tissue"` and `Organ_system = "limb_muscle"`
  (tissues `limb_muscle` and `fast_twitch`). Limb skeletal muscle, the same kind of muscle as
  the vastus lateralis biopsies in GSE109178. Cardiac and smooth muscle were excluded.
- Threshold: **mean RPMM ≥ 10** in either tissue → 353 of 2,656 miRNAs pass.
- Sanity check that the data are right: the highest-expressed miRNAs in this tissue come out as
  miR-486-5p, miR-10b-5p, miR-133a-3p, miR-1-3p, miR-22-3p, miR-378a-3p, miR-206. miR-1,
  miR-133a and miR-206 are the canonical muscle-specific "myomiRs", so the table is behaving.
- miR-29 family in muscle: miR-29a-3p 1,895 RPMM, miR-29b-3p 650, miR-29c-3p 439 → comfortably
  muscle-expressed.
