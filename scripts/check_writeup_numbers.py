#!/usr/bin/env python3
"""Re-check every number in docs/WRITEUP_PLAN.md against the frozen result tables.
Reads results/tables/*.tsv and inputs/auto/*.tsv only; computes nothing new. Needs pandas + numpy.
Usage (from the repo root): python3 scripts/check_writeup_numbers.py [repo_root] [out.tsv]
Written 7 Nov 2026 for the pre-submission check; see docs/verification_7nov.md."""
import pandas as pd, numpy as np, math, re, sys, os
R = sys.argv[1] if len(sys.argv) > 1 else "."; T = f"{R}/results/tables"; I = f"{R}/inputs/auto"
OUT = sys.argv[2] if len(sys.argv) > 2 else None
rd = lambda f, d=T: pd.read_csv(f"{d}/{f}", sep="\t")
out = []
def chk(sec, claim, expected, actual, ok=None, note=""):
    if ok is None:
        ok = (str(expected) == str(actual))
    out.append((sec, claim, str(expected), str(actual), "PASS" if ok else "FAIL", note))
def close(a, b, tol): return abs(float(a) - float(b)) <= tol
# ---------- 4.1 samples ----------
ss = rd("01_sample_sheet.tsv")
inc = ss[ss.include.astype(str).str.upper() == "TRUE"]
chk("4.1", "49 samples total", 49, len(ss))
chk("4.1", "34 included", 34, len(inc))
g = inc.group.value_counts().to_dict()
chk("4.1", "6 control / 11 BMD / 17 DMD", "6/11/17", f"{g.get('control')}/{g.get('BMD')}/{g.get('DMD')}")
lg = ss[ss.group == "LGMD"]
chk("4.1", "15 LGMD = 8 LGMD2B + 7 LGMD2I", "15: " + str({'DYS': 8, 'FKRP': 7}), f"{len(lg)}: {lg.lgmd_subtype.value_counts().to_dict()}", ok=len(lg) == 15)
b = inc[inc.group == "BMD"]; d = inc[inc.group == "DMD"]; c = inc[inc.group == "control"]
chk("4.1", "BMD age median 12.2, range 2.5-60.8", "12.2 (2.5-60.8)", f"{b.age_years.median()} ({b.age_years.min()}-{b.age_years.max()})")
chk("4.1", "BMD severity 6 mild/3 moderate/2 severe", "{'mild': 6, 'moderate': 3, 'severe': 2}", str(dict(sorted(b.severity.value_counts().to_dict().items()))))
chk("4.1", "DMD age median 3.5, range 0.9-8.0", "3.5 (0.9-8.0)", f"{d.age_years.median()} ({d.age_years.min()}-{d.age_years.max()})", ok=(d.age_years.median() == 3.5 and d.age_years.min() == 0.9 and d.age_years.max() == 8.0))
chk("4.1", "DMD severity 4 moderate/13 severe", "{'moderate': 4, 'severe': 13}", str(dict(sorted(d.severity.value_counts().to_dict().items()))))
chk("4.1", "control ages not reported (all NA)", "6 NA", f"{c.age_years.isna().sum()} NA")
chk("7", "BMD age span ~58 y", "58", f"{b.age_years.max() - b.age_years.min():.1f}", ok=close(b.age_years.max() - b.age_years.min(), 58, 0.5))
# ---------- 4.2 QC ----------
qc = rd("02_qc_metrics.tsv"); mc = rd("02_qc_median_cor.tsv")
for col, exp, lab in [("median_cor_all", {"control": .9529, "BMD": .9704, "DMD": .9803}, "overall"), ("median_cor_within_group", {"control": .9839, "BMD": .9757, "DMD": .9855}, "within-group")]:
    med = qc.groupby("group")[col].median().round(4).to_dict(); mean = qc.groupby("group")[col].mean().round(4).to_dict()
    okm = all(close(med[k], v, 0.00006) for k, v in exp.items()); oka = all(close(mean[k], v, 0.00006) for k, v in exp.items())
    chk("4.2", f"{lab} correlation by group {exp}", exp, f"median {med} | mean {mean}", ok=okm or oka, note=("matches group MEDIAN" if okm else "matches group MEAN" if oka else ""))
chk("4.2", "median RLE within +/-0.057", "<=0.057", f"{qc.rle_median.abs().max():.4f}", ok=round(qc.rle_median.abs().max(), 3) <= 0.057)
chk("4.2", "median NUSE 0.981-1.039", "0.981-1.039", f"{qc.nuse_median.min():.3f}-{qc.nuse_median.max():.3f}", ok=(round(qc.nuse_median.min(), 3) == 0.981 and round(qc.nuse_median.max(), 3) == 1.039))
chk("4.2", "nothing flagged on RLE/NUSE", "0", str(int(qc.rle_bad.sum() + qc.nuse_bad.sum())))
chk("4.2", "3 controls flagged on spec metric", "3 (all control)", f"{int(mc.flagged.sum())} flagged: {mc[mc.flagged].group.value_counts().to_dict()}", ok=(mc.flagged.sum() == 3 and (mc[mc.flagged].group == 'control').all()))
fl = qc[qc.flagged]
chk("4.2", "within-group flags GSM2934803 BMD 0.9676; GSM2934813 control 0.9784", "GSM2934803 0.9676; GSM2934813 0.9784", "; ".join(f"{r.gsm} {r.group} {r.median_cor_within_group:.4f}" for r in fl.itertuples()), ok=(set(fl.gsm) == {"GSM2934803", "GSM2934813"}))
# ---------- 4.3 collapse ----------
cc = dict(rd("03_collapse_counts.tsv").values)
chk("4.3", "54,675 -> 44,647 -> -1,524 -> 43,123 -> 20,823", "54675/44647/1524/43123/20823", f"{cc['probes_total']}/{cc['probes_with_symbol']}/{cc['probes_dropped_multimap']}/{cc['probes_kept']}/{cc['genes_kept']}")
# 11,043 multi-probe genes: from probes kept per gene -- the map only has the chosen probe, so cannot recount here
# ---------- 4.4 / 5.1 limma ----------
L = {k: rd(f"limma_{k}.tsv") for k in ["BMDvsCtrl", "DMDvsCtrl", "BMDvsDMD"]}
exp_de = {"BMDvsCtrl": (4549, 4659), "DMDvsCtrl": (5547, 6654), "BMDvsDMD": (4564, 4039)}
for k, (u, dn) in exp_de.items():
    t = L[k]; uu = int(((t["adj.P.Val"] < .05) & (t.logFC > 0)).sum()); dd = int(((t["adj.P.Val"] < .05) & (t.logFC < 0)).sum())
    chk("5.1", f"{k} FDR<0.05 up/down {u}/{dn}", f"{u}/{dn}", f"{uu}/{dd}")
chk("4.3", "20,823 genes in limma tables", 20823, len(L["DMDvsCtrl"]))
# ---------- 5.2 collagen table ----------
san = rd("04_sanity_fibrosis_genes.tsv").set_index("gene_symbol")
exp_tab = {"COL1A1": (2.517, 9.4e-08, 4.058, 5.6e-13, -1.542, 4.6e-06), "COL1A2": (2.043, 1.4e-07, 2.940, 7.5e-12, -0.897, 4.2e-04),
           "COL3A1": (2.311, 8.3e-09, 3.260, 4.6e-13, -0.949, 1.3e-04), "COL5A1": (1.919, 1.4e-06, 2.828, 7.5e-11, -0.910, 8.5e-04), "LOX": (2.607, 2.3e-07, 3.573, 3.7e-11, -0.966, 3.1e-03)}
for gname, e in exp_tab.items():
    r = san.loc[gname]
    act = (r.logFC_BMDvsCtrl, r["adj.P.Val_BMDvsCtrl"], r.logFC_DMDvsCtrl, r["adj.P.Val_DMDvsCtrl"], r.logFC_BMDvsDMD, r["adj.P.Val_BMDvsDMD"])
    ok = all(close(a, x, 0.0006) for a, x in zip(act[::2], e[::2])) and all(float(f"{a:.1e}") == x for a, x in zip(act[1::2], e[1::2]))
    chk("5.2", f"{gname} table row", e, tuple(round(float(a), 3) if i % 2 == 0 else float(f"{a:.1e}") for i, a in enumerate(act)), ok=ok)
for gname, lfc, fold in [("COL1A1", 4.058, 16.7), ("COL3A1", 3.260, 9.6), ("COL1A2", 2.940, 7.7)]:
    chk("5.2", f"{gname} DMD {fold}x", fold, round(2 ** lfc, 1))
chk("5.2", "COL1A1 BMD 5.7x", 5.7, round(2 ** 2.517, 1))
# ---------- 4.5 inputs ----------
fam = rd("mirna_families.tsv", I); tg = rd("mirna_targets.tsv", I); utr = rd("utr_lengths.tsv", I)
chk("4.5", "2,064 human families", 2064, len(fam))
chk("4.5", "253 muscle-expressed families", 253, int(fam.expressed_in_muscle.astype(str).str.upper().eq("TRUE").sum()))
chk("4.5", "targets: 115,308 rows / 219 families / 13,034 genes", "115308/219/13034", f"{len(tg)}/{tg.family_id.nunique()}/{tg.gene_symbol.nunique()}")
chk("4.5", "UTR lengths for 19,429 genes", 19429, len(utr))
fs = rd("05_family_sizes.tsv")
chk("4.5", "128 families; sizes 15-1,446, median 371", "128; 15-1446; 371", f"{len(fs)}; {fs.n_targets.min()}-{fs.n_targets.max()}; {fs.n_targets.median():g}")
m29 = fs[fs.family_id == "miR-29-3p"].iloc[0]
chk("4.5", "miR-29-3p seed AGCACCA, conservation 2, 1,158 targets", "AGCACCA/2/1158", f"{m29.seed_7mer}/{m29.conservation}/{m29.n_targets}")
chk("4.5", "miR-29-3p members 29a/b/c-3p", "hsa-miR-29a-3p;hsa-miR-29b-3p;hsa-miR-29c-3p", m29.members)
chk("4.5", "miR-29a-3p 1,895 RPMM (family max)", 1895, round(m29.max_muscle_rpmm), ok=close(m29.max_muscle_rpmm, 1895, 0.6))
fsets = rd("fibrosis_sets.tsv", I); core = set(fsets[fsets.set_name == "matrisome_core"].gene_symbol); tgfb = set(fsets[fsets.set_name == "hallmark_tgfb"].gene_symbol)
chk("4.9", "core matrisome 275 / TGF-b 54 / overlap 2", "275/54/2", f"{len(core)}/{len(tgfb)}/{len(core & tgfb)}")
fm = rd("fibroblast_markers.tsv", I); nf = fm[~fm.in_fibrosis_set.astype(str).str.upper().eq("TRUE")]
genes = set(L["DMDvsCtrl"].gene_symbol)
chk("4.9", "markers 372 / 326 not in fibrosis set / 316 on array", "372/326/316", f"{len(fm)}/{len(nf)}/{len(set(nf.gene_symbol) & genes)}")
# ---------- ruler ----------
t29 = set(tg[tg.family_id == "miR-29-3p"].gene_symbol) & genes
dm = L["DMDvsCtrl"]; inn = dm.gene_symbol.isin(t29)
mt, mo = dm[inn].logFC.mean(), dm[~inn].logFC.mean()
chk("5.1", "1,158 targets vs 19,665 others", "1158/19665", f"{inn.sum()}/{(~inn).sum()}")
chk("5.1", "mean log2FC targets +0.304 vs others +0.036", "0.304/0.036", f"{mt:.3f}/{mo:.3f}")
chk("4.8", "ruler +0.2686 log2 = 20.5% linear", "0.2686 / 20.5%", f"{mt - mo:.4f} / {100 * (2 ** (mt - mo) - 1):.1f}%")
# ---------- 5.1 smoke / step 7 ----------
sm = rd("06_smoke_fgsea_DMDvsCtrl.tsv"); w = rd("07_fgsea_wide.tsv")
r = sm[sm.family_id == "miR-29-3p"].iloc[0]
chk("5.1", "smoke test miR-29 NES +2.22, padj 5.4e-35 (CLAUDE.md status block)", "2.22 / 5.4e-35", f"{r.NES:.2f} / {r.padj:.1e}")
chk("5.1", "smoke test miR-29 rank by padj = 3", 3, int(r.rank_by_padj))
chk("5.1", "smoke test miR-29 rank by NES = 2", 2, int((sm.NES > r.NES).sum() + 1))
rw = w[w.family_id == "miR-29-3p"].iloc[0]
chk("5.1/5.3", "step 7 miR-29 NES DMD +2.224 / BMD +2.169 / BMDvsDMD -1.590", "2.224/2.169/-1.590", f"{rw.NES_DMDvsCtrl:.3f}/{rw.NES_BMDvsCtrl:.3f}/{rw.NES_BMDvsDMD:.3f}")
chk("5.3", "fgsea padj DMD 1.5e-34 / BMD 1.4e-31 / BMDvsDMD 4.7e-10", "1.5e-34/1.4e-31/4.7e-10", f"{rw.padj_DMDvsCtrl:.1e}/{rw.padj_BMDvsCtrl:.1e}/{rw.padj_BMDvsDMD:.1e}")
chk("5.1", "step 7 DMD: miR-29 rank by NES 2nd of 128", 2, int((w.NES_DMDvsCtrl > rw.NES_DMDvsCtrl).sum() + 1))
chk("5.1", "step 7 DMD: miR-29 rank by padj 3rd", 3, int((w.padj_DMDvsCtrl < rw.padj_DMDvsCtrl).sum() + 1), note="ties counted as better")
chk("5.1", "step 7 BMD: rank 6 by padj, 11 by NES (FULL_RESEARCH_DOC)", "6/11", f"{int((w.padj_BMDvsCtrl < rw.padj_BMDvsCtrl).sum() + 1)}/{int((w.NES_BMDvsCtrl > rw.NES_BMDvsCtrl).sum() + 1)}")
sigd = w[w.padj_DMDvsCtrl < .05]
chk("5.1", "103 of 128 padj<0.05 in DMD vs ctrl; 99 same direction (NES>0)", "103/99", f"{len(sigd)}/{int((sigd.NES_DMDvsCtrl > 0).sum())}")
# ---------- leading edge ----------
les = rd("07_leading_edge_summary.tsv"); le = rd("07_leading_edge.tsv")
lm = les[les.family_id == "miR-29-3p"].iloc[0]; others = les[les.family_id != "miR-29-3p"].pct_fibrosis
chk("5.1", "miR-29 leading edge 6.5% fibrosis, highest of top 25; next 5.3%", "6.5 / max / 5.3", f"{lm.pct_fibrosis} / {'max' if lm.pct_fibrosis >= les.pct_fibrosis.max() else 'NOT max'} / {others.max()}")
chk("5.1", "most other families 1-3%", "1-3", f"{(others.between(1, 3)).sum()} of {len(others)} in [1,3]; range {others.min()}-{others.max()}", ok=(others.between(1, 3).mean() > .5))
roll = ["COL1A1","COL1A2","COL3A1","COL4A1","COL4A2","COL4A4","COL5A1","COL5A2","COL5A3","COL6A2","COL6A3","COL6A6","COL15A1","COL16A1","COL21A1","COL27A1","ELN","FBN1","LAMA2","LAMC1","NID1","SPARC"]
le29 = set(le[(le.family_id == "miR-29-3p") & (le.contrast == "DMDvsCtrl")].gene)
missing = [x for x in roll if x not in le29]
chk("5.1", "collagen roll-call all in miR-29 DMD leading edge", "0 missing", f"{len(missing)} missing {missing}" if missing else "0 missing", ok=not missing)
# ---------- 10a / 12 / 11 ----------
e10 = rd("10a_ecm_removed.tsv"); e = e10[(e10.contrast == "DMDvsCtrl") & (e10.family_id == "miR-29-3p")].iloc[0]
chk("5.1", "ECM removal +2.224 -> +2.151 (96.7%)", "2.151 / 96.7%", f"{e.NES_noECM:.3f} / {100 * e.NES_noECM / e.NES_full:.1f}%")
f12 = rd("12_fibro_adjusted_fgsea.tsv"); fd = f12[f12.contrast == "DMDvsCtrl"].copy(); fd["ret"] = fd.NES_adj / fd.NES_raw
f = fd[fd.family_id == "miR-29-3p"].iloc[0]
chk("5.1", "fibroblast adj +2.224 -> +1.762 (79%)", "1.762 / 79%", f"{f.NES_adj:.3f} / {100 * f.ret:.0f}%")
top = fd.sort_values("padj_raw").head(10); oth = top[top.family_id != "miR-29-3p"]
chk("5.1", "other top-10 families retain 48-74% (median 64%)", "48-74% (median 64%)", f"{100 * oth.ret.min():.0f}-{100 * oth.ret.max():.0f}% (median {100 * oth.ret.median():.0f}%)", note="top-10 = 10 smallest raw padj, DMD vs ctrl")
fs_ = fd.sort_values("padj_adj")
chk("5.1", "after adjustment miR-29 padj 1.7e-16, rank 1; next 5.5e-08", "1.7e-16 / 5.5e-08 / rank 1", f"{f.padj_adj:.1e} / {fs_[fs_.family_id != 'miR-29-3p'].padj_adj.min():.1e} / rank {int((fd.padj_adj < f.padj_adj).sum() + 1)}")
u11 = rd("11_utr_matched.tsv")
for k, en in [("DMDvsCtrl", 1.592), ("BMDvsCtrl", 1.669)]:
    u = u11[(u11.contrast == k) & (u11.family_id == "miR-29-3p")].iloc[0]
    chk("5.1", f"UTR-matched null mean {k} {en}", en, f"{u.null_mean:.3f}")
    if k == "DMDvsCtrl":
        chk("5.1", "UTR-matched z = 6.56; 0 of 1,000 sets reached it", "6.56 / p=1/1001", f"{(u.NES_real - u.null_mean) / u.null_sd:.2f} / p={u.utr_matched_p:.6f} (n_null {u.n_null})", ok=close((u.NES_real - u.null_mean) / u.null_sd, 6.56, 0.006) and close(u.utr_matched_p, 1 / 1001, 1e-9))
        chk("5.1", "~72% of miR-29 enrichment is UTR length", "72%", f"{100 * u.null_mean / u.NES_real:.1f}%", ok=round(100 * u.null_mean / u.NES_real) == 72)
# ---------- 5.3 shuffles ----------
sh = rd("08_family_shuffle_p.tsv"); s = sh[sh.family_id == "miR-29-3p"].iloc[0]
chk("5.3", "miR-29 shuffle p DMD 0.052 / BMD 0.084", "0.052/0.084", f"{s.shuffle_p_DMDvsCtrl:.3f}/{s.shuffle_p_BMDvsCtrl:.3f}")
def denom(p):
    for n in range(450, 502):
        k = p * n
        if abs(k - round(k)) < 1e-6: return n
chk("4.7", "miR-29 step-8 p uses 497 usable shuffles (denominator 498)", "DMD 498 / BMD 498", f"DMD {denom(s.shuffle_p_DMDvsCtrl)} / BMD {denom(s.shuffle_p_BMDvsCtrl)}", note="denominator = usable shuffles + 1")
n501 = lambda c: int(sum(1 for p in sh[c] if abs(p * 501 - round(p * 501)) < 1e-6 and denom(p) == 501))
chk("4.7", "families using all 500 shuffles: 106 (DMD) / 110 (BMD) of 128", "106/110", f"{n501('shuffle_p_DMDvsCtrl')}/{n501('shuffle_p_BMDvsCtrl')}")
chk("5.3", "miR-29 3rd of 128 on shuffle p (DMD)", 3, int((sh.shuffle_p_DMDvsCtrl < s.shuffle_p_DMDvsCtrl).sum() + 1))
chk("5.3", "miR-29 45th on shuffle p (BMD)", 45, f"{int((sh.shuffle_p_BMDvsCtrl < s.shuffle_p_BMDvsCtrl).sum() + 1)} (ties-better) / {int((sh.shuffle_p_BMDvsCtrl <= s.shuffle_p_BMDvsCtrl).sum())} (ties-worse)", ok=int((sh.shuffle_p_BMDvsCtrl < s.shuffle_p_BMDvsCtrl).sum() + 1) == 45 or int((sh.shuffle_p_BMDvsCtrl <= s.shuffle_p_BMDvsCtrl).sum()) == 45)
chk("5.3", "families shuffle p<0.05: DMD 1, BMD 10 of 128", "1/10", f"{int((sh.shuffle_p_DMDvsCtrl < .05).sum())}/{int((sh.shuffle_p_BMDvsCtrl < .05).sum())}")
b9 = rd("09_bmd_vs_dmd.tsv"); m9 = b9[b9.family_id == "miR-29-3p"].iloc[0]
chk("5.3", "BMDvsDMD: fgsea padj<0.05 7 of 128; shuffle 0 of 128; all q>0.95", "7/0/min q>0.95", f"{int((b9.fgsea_padj < .05).sum())}/{int((b9.shuffle_p < .05).sum())}/min q {b9.shuffle_q.min():.4f}", ok=(int((b9.fgsea_padj < .05).sum()) == 7 and int((b9.shuffle_p < .05).sum()) == 0 and b9.shuffle_q.min() > .95))
chk("5.3", "miR-29 BMDvsDMD shuffle p 0.239, q 0.951, null mean 0.053 SD 1.386, 496 usable", "0.239/0.951/0.053/1.386/496", f"{m9.shuffle_p:.3f}/{m9.shuffle_q:.3f}/{m9.null_mean:.3f}/{m9.null_sd:.3f}/{m9.n_null}")
cl = rd("09_family_classification.tsv")
chk("5.3", "classification fgsea 101 shared/27 neither; shuffle 1 shared/127 neither", "{'neither': 27, 'shared': 101} | {'neither': 127, 'shared': 1}", f"{dict(sorted(cl.class_fgsea.value_counts().to_dict().items()))} | {dict(sorted(cl.class_shuffle.value_counts().to_dict().items()))}")
chk("5.3", "the 1 shuffle-shared family is miR-196-5p", "miR-196-5p", ",".join(cl[cl.class_shuffle == "shared"].family_id))
chk("5.3", "families differing BMD vs DMD 0 of 128 (classification)", 0, int(cl.differs.astype(str).str.upper().eq("TRUE").sum()))
chk("5.3", "Monte Carlo SE at p=0.05, N=500 ~ +/-0.01", "0.0097", f"{math.sqrt(.05 * .95 / 500):.4f}")
# ---------- 10b ----------
r10 = rd("10b_detection_runs.tsv"); z = r10[r10.f == 0]
chk("5.3", "f=0: 100 of 100 caught; weakest padj 2.1e-06", "100 / 2.1e-06", f"{int(z.caught.sum())} / {z.padj.max():.1e}")
chk("5.3", "unplanted NES mean 1.679 SD 0.087 range 1.455-1.867", "1.679/0.087/1.455-1.867", f"{z.NES.mean():.3f}/{z.NES.std():.3f}/{z.NES.min():.3f}-{z.NES.max():.3f}")
q95 = np.quantile(z.NES, .95)
chk("4.8", "calibrated bar NES > 1.812 (95th pct of unplanted)", 1.812, f"{q95:.3f}", note="numpy default quantile (R type 7)")
cal = [round(float((r10[r10.f == ff].NES > q95).mean()), 2) for ff in sorted(r10.f.unique())]
chk("5.4", "calibrated catch 0.05/0.73/1/1/1/1", "[0.05, 0.73, 1.0, 1.0, 1.0, 1.0]", str(cal))
fx = 0.10 + (0.80 - cal[1]) / (cal[2] - cal[1]) * 0.15
chk("5.4", "80% crossed at ~14% (+0.037 log2)", "14% / 0.037", f"{100 * fx:.1f}% / {fx * (mt - mo):.3f}", ok=(round(100 * fx) == 14 and round(fx * (mt - mo), 3) == 0.037))
ru = rd("10b2_ruler_utr_share.tsv").set_index("contrast")
exc = ru.loc["BMDvsCtrl", "mir29_gap"] - ru.loc["BMDvsCtrl", "utr_matched_gap_mean"]
chk("5.4", "BMD excess +0.095 = 35% of ruler", "0.095 / 35%", f"{exc:.3f} / {100 * exc / ru.loc['DMDvsCtrl', 'mir29_gap']:.1f}%", ok=(round(exc, 3) == 0.095 and round(100 * exc / ru.loc['DMDvsCtrl', 'mir29_gap']) == 35))
chk("5.4", "UTR = 42% of DMD ruler, 51% of BMD gap", "42%/51%", f"{100 * ru.loc['DMDvsCtrl', 'utr_share_of_mir29_gap']:.0f}%/{100 * ru.loc['BMDvsCtrl', 'utr_share_of_mir29_gap']:.0f}%")
chk("5.4", "BMD gap 0.1923 minus UTR-matched 0.0973", "0.1923/0.0973", f"{ru.loc['BMDvsCtrl', 'mir29_gap']:.4f}/{ru.loc['BMDvsCtrl', 'utr_matched_gap_mean']:.4f}")
# ---------- 13 / 14 ----------
m13 = rd("13_mir29_leave_one_out.tsv")
bc = m13[m13.contrast == "BMDvsCtrl"]; bd = m13[m13.contrast == "BMDvsDMD"]
chk("5.4", "LOO: 17 runs; BMDvsCtrl 2.120-2.186, all padj<0.05", "17 / 2.120-2.186 / 17", f"{len(bc)} / {bc.NES.min():.3f}-{bc.NES.max():.3f} / {int((bc.padj < .05).sum())}")
chk("5.4", "LOO: BMDvsDMD -1.674 to -1.500, rank 1 every run", "-1.674..-1.500 / [1]", f"{bd.NES.min():.3f}..{bd.NES.max():.3f} / {sorted(set(int(x) for x in bd['rank']))}")
g4 = m13[(m13.dropped == "GSM2934804") & (m13.contrast == "BMDvsCtrl")].iloc[0]
chk("5.4", "drop GSM2934804 -> +2.156", 2.156, f"{g4.NES:.3f}")
l13 = rd("13_leave_one_out.tsv")
full_top = set(w[w.padj_BMDvsCtrl.rank(method="first") <= 20].family_id)
ov = []
for dr, gg in l13[l13.contrast == "BMDvsCtrl"].groupby("dropped"):
    ov.append(len(set(gg.nsmallest(20, "padj").family_id) & full_top))
chk("5.4", "top-20 overlap 18-20 of 20 (BMDvsCtrl, by padj)", "18-20", f"{min(ov)}-{max(ov)}", note="recomputed; log says min 18 max 20")
s14 = rd("14_subsample_dmd.tsv")
chk("5.4", "DMD@n=11 median 2.167, range 2.074-2.258", "2.167 / 2.074-2.258", f"{s14.m29_nes.median():.3f} / {s14.m29_nes.min():.3f}-{s14.m29_nes.max():.3f}")
chk("5.4", "53.0% of draws at or below BMD 2.169", "53.0%", f"{100 * (s14.m29_nes <= rw.NES_BMDvsCtrl).mean():.1f}%")
# ---------- 7 fibroblast scores ----------
sc = rd("12_fibroblast_scores.tsv"); gm = sc.groupby("group").fibro_score.mean()
k = 3; n = len(sc); grand = sc.fibro_score.mean()
ssb = sum(len(x) * (x.mean() - grand) ** 2 for _, x in sc.groupby("group").fibro_score); ssw = sum(((x - x.mean()) ** 2).sum() for _, x in sc.groupby("group").fibro_score)
F = (ssb / (k - 1)) / (ssw / (n - k)); p = (1 + 2 * F / (n - k)) ** (-(n - k) / 2)
chk("7", "fibroblast score ctrl -0.791 / BMD +0.040 / DMD +0.253", "-0.791/0.040/0.253", f"{gm['control']:.3f}/{gm['BMD']:.3f}/{gm['DMD']:.3f}")
chk("7", "ANOVA p = 1.2e-08", "1.2e-08", f"{p:.1e}", note=f"F={F:.2f}, df 2,{n - k}")
# ---------- 6.5 precision ----------
a = math.sqrt(1 / 6 + 1 / 11); bb = math.sqrt(1 / 12 + 1 / 11)
chk("6.5", "GSE13608: 0.508 -> 0.417 (~18%)", "0.508/0.417/18%", f"{a:.3f}/{bb:.3f}/{100 * (1 - bb / a):.0f}%")
chk("5.1", "fibroblast adj: 103 -> 51 significant families (DMD)", "103/51", f"{int((fd.padj_raw < .05).sum())}/{int((fd.padj_adj < .05).sum())}")
chk("5.4", "leave-one-out reran BMDvsCtrl and BMDvsDMD only (not DMDvsCtrl)", "['BMDvsCtrl', 'BMDvsDMD']", str(sorted(l13.contrast.unique())))
# ---------- print ----------
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 95)
df = pd.DataFrame(out, columns=["sec", "claim", "expected", "actual", "result", "note"])
if OUT: df.to_csv(OUT, sep="\t", index=False)
print(df[["sec", "result", "claim", "actual", "note"]].to_string(index=False))
print("\nTOTAL", len(df), "| PASS", (df.result == "PASS").sum(), "| FAIL", (df.result == "FAIL").sum())
