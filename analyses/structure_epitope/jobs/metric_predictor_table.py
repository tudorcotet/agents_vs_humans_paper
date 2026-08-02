"""Per-predictor table of AUROC (binder vs non, P_expressed), Spearman ρ (vs pKd, P_kd), and BH-FDR
for every metric named in the section Summary.

The 'predictor' axis is intentionally ragged, because the metrics live on different models:
  - Interface size / chemistry / secondary structure  -> Protenix, Chai, ESMFold2, Boltz-2, AF2M
  - ipTM / pLDDT / ipSAE / PRODIGY ΔG                 -> Boltz-2, Protenix, Chai, AF2M (folding-model scores)
  - ESM-2 pll, SaProt                                 -> sequence-based, one value (no predictor axis)
  - contact motifs                                    -> population enrichment (OR/p), NOT a per-design score

FDR: BH within tier, per analysis (matches the draft's two-tier structure, extended across predictors):
  Tier "learned/sequence scores" = ipTM, pLDDT, ipSAE, ESM-2, SaProt
  Tier "structural/physical"      = interface size, interface chemistry, PRODIGY ΔG, secondary structure
ρ p-value is the Spearman asymptotic p (the draft uses a 10k permutation p; they agree closely at this n).

Run: uv run --with pandas --with numpy --with scipy --with scikit-learn --with pyarrow \
       python analyses/structure_epitope/jobs/metric_predictor_table.py
Out: results/metric_predictor_table.md  +  results/metric_predictor_table.csv
"""
from __future__ import annotations
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
from sklearn.metrics import roc_auc_score

ROOT = Path("/Users/amyhe/Desktop/trem2_2026/agents_vs_humans_paper")
RES = ROOT / "analyses/structure_epitope/results"

# ---- data sources ----
meta = pd.read_parquet(ROOT / "analyses/structure_epitope/data/designs_canonical.parquet")[
    ["design_id", "is_hit", "p_expressed", "p_kd", "pkd_arith_mean"]]
dp = pd.read_parquet(ROOT / "data/designs.parquet")
# px_mean_plddt and af2m_mean_plddt ship on a 0-100 scale while the other pLDDT columns are 0-1; rescale.
for _c in ["px_mean_plddt", "af2m_mean_plddt"]:
    if dp[_c].median() > 1.5:
        dp[_c] = dp[_c] / 100.0
lp = pd.read_csv(RES / "interface_descriptors_perpredictor.csv")
mop = pd.read_csv(RES / "motifs_L1_architecture_perpredictor.csv")  # per-predictor helix/strand/loop

# ---- metric registry: (tier, family, label, {predictor: (source, column)}) ----
IF = "iface"; SC = "scalar"; MOP = "motif_perpred"
STRUCT, LEARN = "structural/physical", "learned/sequence"
GEOM = ["Protenix", "Chai", "ESMFold2", "Boltz-2", "AF2M"]  # predictors we computed complex geometry on
PKEY = {"Protenix": "protenix", "Chai": "chai", "ESMFold2": "esmfold2", "Boltz-2": "boltz2", "AF2M": "af2m"}  # display->data key
REG = [
    (STRUCT, "Interface size", "buried area (BSA)",      {p: (IF, "bsa") for p in GEOM}),
    (STRUCT, "Interface size", "binder iface residues",  {p: (IF, "n_iface_binder") for p in GEOM}),
    (STRUCT, "Interface size", "target iface residues",  {p: (IF, "n_iface_target") for p in GEOM}),
    (STRUCT, "Interface chemistry", "hydrophobic frac",  {p: (IF, "iface_frac_hydrophobic") for p in GEOM}),
    (STRUCT, "Interface chemistry", "polar frac",        {p: (IF, "iface_frac_polar") for p in GEOM}),
    (STRUCT, "Interface chemistry", "charged frac",      {p: (IF, "iface_frac_charged") for p in GEOM}),
    (STRUCT, "Interface chemistry", "H-bond count",      {p: (IF, "hbonds") for p in GEOM}),
    (STRUCT, "Interface chemistry", "salt bridges",      {p: (IF, "salt_bridges") for p in GEOM}),
    (STRUCT, "Interface chemistry", "H-bond density",    {p: (IF, "hbond_density_per100") for p in GEOM}),
    (STRUCT, "Secondary structure", "interface helix frac",  {p: (MOP, "iface_helix") for p in GEOM}),
    (STRUCT, "Secondary structure", "interface strand frac", {p: (MOP, "iface_strand") for p in GEOM}),
    (STRUCT, "Contacts ΔG (PRODIGY)", "PRODIGY pKd",     {"Boltz-2": (SC, "prodigy_boltz2_pkd"), "Protenix": (SC, "prodigy_protenix_pkd"), "Chai": (SC, "prodigy_chai_pkd"), "AF2M": (SC, "prodigy_af2m_pkd")}),
    (LEARN, "Confidence: ipTM", "ipTM",                  {"Boltz-2": (SC, "boltz2_iptm"), "Protenix": (SC, "px_iptm"), "Chai": (SC, "chai_iptm"), "AF2M": (SC, "af2m_iptm")}),
    (LEARN, "Confidence: pLDDT", "pLDDT (mean)",         {"Boltz-2": (SC, "b2_mean_plddt"), "Protenix": (SC, "px_mean_plddt"), "Chai": (SC, "chai_mean_plddt"), "AF2M": (SC, "af2m_mean_plddt")}),
    (LEARN, "Confidence: ipSAE", "ipSAE",                {"Boltz-2": (SC, "submitted_ipsae"), "Protenix": (SC, "px_ipsae_d0chn_max"), "Chai": (SC, "chai_ipsae_d0chn_max"), "AF2M": (SC, "af2m_ipsae_d0chn_max")}),
    (LEARN, "Sequence PLL", "ESM-2 pll",                 {"(sequence)": (SC, "esm_pll_avg")}),
    (LEARN, "Sequence PLL", "SaProt pll",                {"(sequence)": (SC, "saprot_pll_norm")}),
]
PRED_ORDER = ["Boltz-2", "Protenix", "Chai", "ESMFold2", "Protenix+Chai", "(sequence)"]


def series_for(source, col, predictor):
    if source == IF:
        s = lp[lp.predictor == PKEY.get(predictor, predictor.lower())][["design_id", col]]  # lp keys are lowercase
    elif source == MOP:
        s = mop[mop.pred == PKEY.get(predictor, predictor.lower())][["design_id", col]]
    else:
        s = dp[["design_id", col]]
    return s.rename(columns={col: "val"}).merge(meta, on="design_id", how="right")


def auroc_disc(df):
    d = df[df.p_expressed == True].dropna(subset=["val"])            # noqa: E712
    b = d[d.is_hit == True]["val"]; n = d[d.is_hit == False]["val"]  # noqa: E712
    if len(b) < 3 or len(n) < 3:
        return np.nan, np.nan, np.nan
    p = stats.mannwhitneyu(b, n, alternative="two-sided").pvalue
    y = (d.is_hit == True).astype(int)                              # noqa: E712
    auc = max(roc_auc_score(y, d.val), roc_auc_score(y, -d.val))
    return round(auc, 3), p, len(d)


def rho_aff(df):
    d = df[df.p_kd == True].dropna(subset=["val", "pkd_arith_mean"])  # noqa: E712
    if len(d) < 6:
        return np.nan, np.nan, np.nan
    rho, p = stats.spearmanr(d.val, d.pkd_arith_mean)
    return round(rho, 3), p, len(d)


def bh(pvals):
    p = np.asarray(pvals, float); ok = ~np.isnan(p); q = np.full(len(p), np.nan)
    idx = np.where(ok)[0]; pp = p[ok]; n = len(pp)
    if n:
        o = np.argsort(pp); r = pp[o] * n / (np.arange(n) + 1)
        qq = np.minimum.accumulate(r[::-1])[::-1]
        out = np.empty(n); out[o] = np.clip(qq, 0, 1); q[idx] = out
    return q


rows = []
for tier, family, label, preds in REG:
    for predictor, (src, col) in preds.items():
        s = series_for(src, col, predictor)
        auc, pd_, nd = auroc_disc(s)
        rho, pa, na = rho_aff(s)
        rows.append(dict(tier=tier, family=family, metric=label, predictor=predictor, column=col,
                         auroc=auc, auroc_p=pd_, n_disc=nd, rho=rho, rho_p=pa, n_kd=na))
R = pd.DataFrame(rows)
# BH within tier, per analysis
for tier in R.tier.unique():
    m = R.tier == tier
    R.loc[m, "q_disc"] = bh(R.loc[m, "auroc_p"].values)
    R.loc[m, "q_aff"] = bh(R.loc[m, "rho_p"].values)
R.to_csv(RES / "metric_predictor_table.csv", index=False)


def cell(r, kind):
    if kind == "auroc":
        v, q = r.auroc, r.q_disc
    else:
        v, q = r.rho, r.q_aff
    if pd.isna(v):
        return "—"
    star = "*" if (not pd.isna(q) and q < 0.10) else ""
    return f"{v:+.2f}{star}".replace("+", "") if kind == "auroc" else f"{v:+.2f}{star}"


def qcell(r, kind):
    q = r.q_disc if kind == "auroc" else r.q_aff
    return "—" if pd.isna(q) else f"{q:.3f}"


L = ["# Per-predictor metric table — AUROC, ρ, and FDR\n",
     "For every metric named in the Summary. **AUROC** = binder-vs-non-binder discrimination on "
     "P_expressed (n=89; 37 binders/52 non); **ρ** = Spearman vs pKd on P_kd (n=36); **q** = BH-FDR "
     "within tier, per analysis (this is a broader family than the per-section panels, so these q's run "
     "slightly more conservative than the draft's headline values — e.g. Protenix BSA q=0.014 here vs "
     "0.004 in §2's 9-descriptor panel; the raw p-values and effect sizes are identical). `*` marks "
     "q<0.10. Cells are `—` where the metric does not exist for "
     "that predictor (the axis is ragged by construction — see header of `jobs/metric_predictor_table.py`).\n",
     "Predictors: **Protenix** is the section primary. Interface geometry + secondary structure are on "
     "Protenix/Chai/ESMFold2/Boltz-2/AF2M (Boltz-2 + AF2M = regenerated/PR#6 complexes; AF2M ships as PDB); "
     "learned confidence + PRODIGY on Boltz-2/Protenix/Chai/AF2M; ESM-2/SaProt are sequence-based (one "
     "value); **contact motifs are "
     "population enrichment (odds ratios), not a per-design score, so they are reported separately below, "
     "not as AUROC/ρ.** pLDDT is **mean-complex pLDDT** on 0–1 for all three (Protenix rescaled from "
     "0–100) — a slightly different quantity from §5's Boltz-2 *binder* pLDDT (ρ=0.14).\n"]

for analysis, kind, head in [("Binder-vs-non-binder discrimination", "auroc", "AUROC"),
                             ("Affinity (pKd) correlation", "rho", "ρ")]:
    L.append(f"## {analysis} — {head} (q)\n")
    L.append("| tier | metric | Boltz-2 | Protenix | Chai | ESMFold2 | AF2M | sequence |")
    L.append("|---|---|---|---|---|---|---|---|")
    for (tier, family, metric), g in R.groupby(["tier", "family", "metric"], sort=False):
        by = {r.predictor: r for r in g.itertuples()}
        def two(pcol):
            if pcol not in by:
                return "—"
            r = by[pcol]
            return f"{cell(r, kind)} ({qcell(r, kind)})"
        L.append(f"| {tier} | {metric} | {two('Boltz-2')} | {two('Protenix')} | {two('Chai')} | "
                 f"{two('ESMFold2')} | {two('AF2M')} | {two('(sequence)')} |")
    L.append("")

# ---- signed-difference table: Protenix (reference) − every other predictor, per design ----
def per_design(source, col, predictor):
    if source == IF:
        s = lp[lp.predictor == PKEY.get(predictor, predictor.lower())][["design_id", col]]
    elif source == MOP:
        s = mop[mop.pred == PKEY.get(predictor, predictor.lower())][["design_id", col]]
    else:
        s = dp[["design_id", col]]
    return s.rename(columns={col: "v"}).dropna()

sd = []
for tier, family, metric, preds in [(t, f, m, p) for (t, f, m, p) in REG if "Protenix" in p]:
    srcP, colP = preds["Protenix"]
    P = per_design(srcP, colP, "Protenix").rename(columns={"v": "vP"})
    for other in [o for o in ["Chai", "ESMFold2", "Boltz-2", "AF2M"] if o in preds]:
        srcO, colO = preds[other]
        O = per_design(srcO, colO, other).rename(columns={"v": "vO"})
        mrg = P.merge(O, on="design_id")
        if len(mrg) < 5:
            continue
        d = (mrg.vP - mrg.vO).values
        med = float(np.median(d)); pmed = float(np.median(mrg.vP))
        rel = 100 * med / pmed if pmed else np.nan
        wp = stats.wilcoxon(d).pvalue if np.any(d != 0) else np.nan
        sd.append(dict(family=family, metric=metric, vs=other, n=len(mrg),
                       median_diff=round(med, 3), rel_pct=round(rel, 1),
                       prot_gt_pct=round(100 * np.mean(d > 0)), wilcoxon_p=wp))
SD = pd.DataFrame(sd)
SD.to_csv(RES / "metric_predictor_signed_diff.csv", index=False)

L.append("## Signed difference from Protenix (reference) — median(Protenix − other), per design\n")
L.append("Positive = Protenix reads *higher* than the other predictor. `rel%` = median diff ÷ Protenix "
         "median. `Prot>other%` = fraction of designs where Protenix is higher. Wilcoxon signed-rank p "
         "(paired, all designs with both predictors).\n")
L.append("| metric | vs predictor | n | median Δ (Prot−other) | rel% | Prot>other% | Wilcoxon p |")
L.append("|---|---|---|---|---|---|---|")
for r in SD.itertuples():
    wp = "—" if pd.isna(r.wilcoxon_p) else f"{r.wilcoxon_p:.1e}"
    rel = "—" if pd.isna(r.rel_pct) else f"{r.rel_pct:+g}%"
    L.append(f"| {r.metric} | {r.vs} | {r.n} | {r.median_diff:+g} | {rel} | {r.prot_gt_pct:.0f}% | {wp} |")
L.append("")

L.append("## Contact motifs (not a per-design AUROC/ρ metric)\n")
L.append("Level-2/3 typed contact-pair analysis is population enrichment, not a per-design score: "
         "aromatic-to-apex OR 3.1 (p=0.025), Arg/Lys-to-tip OR 2.3 (p=0.085), **none surviving FDR**; "
         "no convergent residue-pair motif beyond universal apical engagement "
         "(`results/motifs_L2_contactpairs.csv`, `motifs_L3_convergent.csv`). Computed on Protenix+Chai.\n")

(RES / "metric_predictor_table.md").write_text("\n".join(L) + "\n")
print("\n".join(L))
print("wrote metric_predictor_table.md + .csv")
