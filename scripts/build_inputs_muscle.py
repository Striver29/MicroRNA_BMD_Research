#!/usr/bin/env python3
"""Collapse the miRNA Tissue Atlas parquet to one max-RPMM value per miRNA in skeletal muscle.
Skeletal muscle = Biotype 'tissue', Organ_system 'limb_muscle' (Tissue: limb_muscle, fast_twitch).
Writes data/raw/mirna_atlas/muscle_rpmm.csv with columns mirna, rpmm_muscle_max."""
import sys, pyarrow.parquet as pq, pandas as pd
src = "data/raw/mirna_atlas/all_human_miRNA_rpmm_harmonized_meta_data.parquet"
d = pq.read_table(src).to_pandas()
m = d[(d.Biotype == "tissue") & (d.Organ_system == "limb_muscle") & (d.Metric == "mean")]
assert len(m) == 2, f"expected 2 skeletal-muscle rows, got {len(m)}: {list(m.Tissue)}"
cols = [c for c in m.columns if c.startswith("hsa-")]
out = m[cols].max(axis=0).rename_axis("mirna").reset_index(name="rpmm_muscle_max")
out.to_csv("data/raw/mirna_atlas/muscle_rpmm.csv", index=False)
print(f"tissues used: {list(m.Tissue)}")
print(f"wrote {len(out)} miRNAs; >=10 RPMM: {(out.rpmm_muscle_max>=10).sum()}")
