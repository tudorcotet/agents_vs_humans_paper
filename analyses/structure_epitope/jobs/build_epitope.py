"""
Phase 3 core: consensus footprints, per-residue engagement frequency, and the pre-registered
H-Perera enrichment test. Reads results/footprints_long.parquet + the canonical table.

Consensus footprint of a design = TREM2 residues in the PRIMARY footprint (contact_45 & dASA>1)
of BOTH clean predictors (Protenix AND Chai-1). Predictor agreement reported.

Outputs:
  results/epitope_per_residue_freq.csv   (per construct residue x population: contact frequency)
  results/epitope_predictor_agreement.csv
  results/epitope_HPerera.csv            (aggregate + per-residue Fisher, BH-FDR)
  results/footprints_consensus.parquet   (design x residue consensus indicator)
"""
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats

def bh_fdr(pvals):
    """Benjamini-Hochberg adjusted p-values (q). Vectorized, no statsmodels dependency."""
    p = np.asarray(pvals, float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n); out[order] = np.clip(q, 0, 1)
    return out

ROOT = Path("/Users/amyhe/Desktop/trem2_2026/agents_vs_humans_paper")
RES = ROOT / "analyses/structure_epitope/results"
fp = pd.read_parquet(RES / "footprints_long.parquet")
meta = pd.read_parquet(ROOT / "analyses/structure_epitope/data/designs_canonical.parquet")[
    ["design_id", "cohort", "is_hit", "expressed", "p_screened", "p_expressed", "p_kd",
     "kd_arith_mean_nM_all", "pkd_arith_mean", "team"]]

# hydrophobic tip (UniProt -> construct = U-17)
TIP_U = [44, 69, 70, 71, 74, 89]
TIP_C = [u - 17 for u in TIP_U]
BASIC_U = [46, 47, 62, 76, 77, 87, 122, 123]
BASIC_C = [u - 17 for u in BASIC_U]

# ---- pivot to per-design per-residue primary, per predictor ----
prim = fp.pivot_table(index=["design_id", "construct_pos", "uniprot", "region"],
                      columns="predictor", values="primary", aggfunc="first").reset_index()
prim["consensus"] = prim.get("protenix", False).fillna(False) & prim.get("chai", False).fillna(False)
prim["union"] = prim.get("protenix", False).fillna(False) | prim.get("chai", False).fillna(False)

# ---- predictor agreement: per-design Jaccard of protenix vs chai primary footprints ----
agg = []
for did, g in prim.groupby("design_id"):
    a = set(g.loc[g["protenix"] == True, "construct_pos"])
    b = set(g.loc[g["chai"] == True, "construct_pos"])
    inter = len(a & b); uni = len(a | b)
    agg.append({"design_id": did, "n_protenix": len(a), "n_chai": len(b),
                "n_consensus": inter, "jaccard": (inter / uni if uni else np.nan)})
agree = pd.DataFrame(agg).merge(meta[["design_id", "cohort", "is_hit"]], on="design_id")
agree.to_csv(RES / "epitope_predictor_agreement.csv", index=False)
print("=== predictor agreement (protenix vs chai primary footprint) ===")
print(f"median Jaccard = {agree.jaccard.median():.3f}  (binders {agree[agree.is_hit==True].jaccard.median():.3f} / "
      f"non {agree[agree.is_hit==False].jaccard.median():.3f})")
print(f"median consensus footprint size = {agree.n_consensus.median():.0f} residues "
      f"(protenix {agree.n_protenix.median():.0f}, chai {agree.n_chai.median():.0f})")

# consensus footprint table
cons = prim[["design_id", "construct_pos", "uniprot", "region", "consensus"]].copy()
cons.to_parquet(RES / "footprints_consensus.parquet", index=False)

# ---- per-residue engagement frequency across populations (consensus footprint, IgSF+stalk) ----
cons = cons.merge(meta, on="design_id")
epi = cons[cons.region.isin(["IgSF", "stalk"])].copy()
def freq_table(sub, label):
    ids = sub.design_id.nunique()
    t = (sub.groupby(["construct_pos", "uniprot", "region"])["consensus"].mean()
         .rename(f"freq_{label}").reset_index())
    t[f"n_{label}"] = ids
    return t
pops = {
    "screened": epi[epi.p_screened == True],
    "expressed": epi[epi.p_expressed == True],
    "binder": epi[(epi.p_expressed == True) & (epi.is_hit == True)],
    "nonbinder": epi[(epi.p_expressed == True) & (epi.is_hit == False)],
    "human": epi[(epi.p_expressed == True) & (epi.cohort == "human")],
    "agent": epi[(epi.p_expressed == True) & (epi.cohort == "agent")],
}
base = None
for lab, sub in pops.items():
    t = freq_table(sub, lab)
    base = t if base is None else base.merge(t.drop(columns=[c for c in ["region"] if c in t]), on=["construct_pos", "uniprot"], how="outer")
base = base.sort_values("construct_pos")
base["diff_binder_minus_nonbinder"] = base["freq_binder"] - base["freq_nonbinder"]
base.to_csv(RES / "epitope_per_residue_freq.csv", index=False)
print("\n=== top 12 residues by binder contact frequency (expressed set) ===")
show = base.sort_values("freq_binder", ascending=False).head(12)
print(show[["uniprot", "region", "freq_binder", "freq_nonbinder", "diff_binder_minus_nonbinder"]].to_string(index=False))

# ---- H-Perera: aggregate tip contact, binder vs non-binder within P_expressed ----
exp = meta[meta.p_expressed == True].copy()
# per-design: does consensus footprint hit >=1 tip residue? basic? any IgSF?
def design_hits(residset):
    d = cons[cons.consensus == True].groupby("design_id")["construct_pos"].apply(set)
    return {did: bool(set(residset) & s) for did, s in d.items()}
tip_hit = design_hits(TIP_C); basic_hit = design_hits(BASIC_C)
exp["tip"] = exp.design_id.map(lambda d: tip_hit.get(d, False))
exp["basic"] = exp.design_id.map(lambda d: basic_hit.get(d, False))
exp["binder"] = exp.is_hit == True

def fisher2x2(flag):
    a = int(((exp[flag]) & (exp.binder)).sum())      # contact & binder
    b = int(((exp[flag]) & (~exp.binder)).sum())     # contact & non
    c = int(((~exp[flag]) & (exp.binder)).sum())     # no & binder
    d = int(((~exp[flag]) & (~exp.binder)).sum())    # no & non
    OR, p = stats.fisher_exact([[a, b], [c, d]])
    return a, b, c, d, OR, p
rows = []
for flag, name in [("tip", "hydrophobic_tip_aggregate"), ("basic", "basic_patch_aggregate")]:
    a, b, c, d, OR, p = fisher2x2(flag)
    rows.append({"test": name, "binder_contact": a, "nonbinder_contact": b,
                 "binder_nocontact": c, "nonbinder_nocontact": d, "odds_ratio": round(OR, 3),
                 "p": round(p, 4)})
agg_df = pd.DataFrame(rows)
print("\n=== H-Perera aggregate tests (within P_expressed, n=89) ===")
print(agg_df.to_string(index=False))

# ---- per-residue Fisher (contacted x binder) over IgSF, BH-FDR ----
per = []
exp_ids = set(exp.design_id)
cons_hit_by_res = (cons[(cons.consensus == True) & (cons.design_id.isin(exp_ids))]
                   .groupby("construct_pos")["design_id"].apply(set))
binder_ids = set(exp[exp.binder].design_id); nonb_ids = set(exp[~exp.binder].design_id)
igsf_res = sorted(cons[cons.region == "IgSF"].construct_pos.unique())
for cpos in igsf_res:
    hit = cons_hit_by_res.get(cpos, set())
    a = len(hit & binder_ids); c = len(binder_ids) - a
    b = len(hit & nonb_ids); d = len(nonb_ids) - b
    OR, p = stats.fisher_exact([[a, b], [c, d]])
    per.append({"construct_pos": cpos, "uniprot": cpos + 17,
                "binder_contact": a, "nonbinder_contact": b, "odds_ratio": OR, "p": p})
perdf = pd.DataFrame(per)
perdf["q_bh"] = bh_fdr(perdf.p.values)
perdf = perdf.sort_values("p")
perdf.to_csv(RES / "epitope_HPerera_perresidue.csv", index=False)
agg_df.to_csv(RES / "epitope_HPerera.csv", index=False)
print("\n=== per-residue enrichment: residues with q_BH < 0.05 ===")
sig = perdf[perdf.q_bh < 0.05]
if len(sig):
    print(sig[["uniprot", "binder_contact", "nonbinder_contact", "odds_ratio", "p", "q_bh"]].to_string(index=False))
else:
    print("NONE survive q<0.05.")
print("\ntop 10 by raw p:")
print(perdf.head(10)[["uniprot", "binder_contact", "nonbinder_contact", "odds_ratio", "p", "q_bh"]].to_string(index=False))
