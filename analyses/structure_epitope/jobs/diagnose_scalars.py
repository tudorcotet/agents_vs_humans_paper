"""
Diagnostic: are the b2_* (Boltz-2) / af2m_* (AF2M) SCALAR metric columns corrupt (computed from the
mislabeled complex files), or fine (computed upstream on correct structures)? Decisive if the
per-design scalars track the mislabeled-binder identity rather than the design's own binder.
"""
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
from sklearn.metrics import roc_auc_score
ROOT = Path("/Users/amyhe/Desktop/trem2_2026/agents_vs_humans_paper")
RES = ROOT / "analyses/structure_epitope/results"
df = pd.read_parquet(ROOT / "data/designs.parquet")
audit = pd.read_csv(RES / "structure_integrity_audit.csv")
meta = pd.read_parquet(ROOT / "analyses/structure_epitope/data/designs_canonical.parquet")

def matched(pred):
    a = audit[audit.predictor == pred].copy()
    def mid(r):
        if r.status == "self_match": return r.design_id
        s = str(r.matched_id)
        for tok in s.replace("[","").replace("]","").split(","):
            tok = tok.strip()
            if tok.isdigit(): return int(tok)
        return np.nan
    a["file_binder_id"] = a.apply(mid, axis=1)
    return a[["design_id","status","file_binder_id"]]

print("################ Is b2_* / af2m_* corrupt? ################")
for pred, ipcol, iptcol in [("boltz2","b2_ipsae_d0chn_max","b2_iptm"), ("af2m",None,"af2m_iptm")]:
    print(f"\n===== {pred} =====")
    m = matched(pred).merge(df[["design_id",iptcol]+([ipcol] if ipcol else [])], on="design_id")
    # (C) duplication test: designs whose FILE contains the same wrong binder -> identical scalar?
    grp = m[m.status=="match_OTHER_design"].groupby("file_binder_id")
    ident_groups = 0; varied_groups = 0; examples=[]
    for fid, g in grp:
        if len(g) < 2: continue
        vals = g[iptcol].round(4).nunique()
        if vals == 1: ident_groups += 1
        else: varied_groups += 1
        if len(examples)<3: examples.append((int(fid), len(g), list(g[iptcol].round(3))[:5]))
    print(f"  duplication test: groups of designs sharing one wrong binder — "
          f"IDENTICAL {iptcol}: {ident_groups} groups | VARIED: {varied_groups} groups")
    print(f"    (IDENTICAL within groups => scalars computed FROM the corrupt files = CORRUPT)")
    for fid,n,vals in examples: print(f"    file-binder {fid}: {n} designs, {iptcol}={vals}")
    # (A/B) agreement with competition + clean predictors
    sub = df[["design_id","submitted_ipsae","px_iptm","chai_iptm",iptcol]].dropna()
    r_comp = stats.spearmanr(sub[iptcol], sub.submitted_ipsae)[0]
    r_px = stats.spearmanr(sub[iptcol], sub.px_iptm)[0]
    r_chai = stats.spearmanr(sub[iptcol], sub.chai_iptm)[0]
    print(f"  {iptcol} vs submitted_ipsae ρ={r_comp:.2f} | vs px_iptm ρ={r_px:.2f} | vs chai_iptm ρ={r_chai:.2f}")
    # (D) binder classifier AUROC (expected among expressed) — df has is_hit/expressed
    e = df[(df.expressed==True) & df[iptcol].notna()]
    y=(e.is_hit==True).astype(int); auc=max(roc_auc_score(y,e[iptcol]),roc_auc_score(y,-e[iptcol]))
    print(f"  {iptcol} binder AUROC (expressed) = {auc:.2f}")

print("\n=== reference: clean predictors agree with each other + competition ===")
sub = df[["submitted_ipsae","px_iptm","chai_iptm"]].dropna()
print(f"  px_iptm vs chai_iptm ρ={stats.spearmanr(sub.px_iptm,sub.chai_iptm)[0]:.2f} | "
      f"px_iptm vs submitted_ipsae ρ={stats.spearmanr(sub.px_iptm,sub.submitted_ipsae)[0]:.2f}")
# pb_boltz2 (ProteinBase Boltz-2 rerun, screened) vs b2
if "pb_boltz2_iptm" in df:
    s2=df[["pb_boltz2_iptm","b2_iptm","submitted_ipsae"]].dropna()
    print(f"  pb_boltz2_iptm vs submitted_ipsae ρ={stats.spearmanr(s2.pb_boltz2_iptm,s2.submitted_ipsae)[0]:.2f} "
          f"(ProteinBase Boltz-2 rerun, screened) | pb_boltz2 vs b2 ρ={stats.spearmanr(s2.pb_boltz2_iptm,s2.b2_iptm)[0]:.2f}")
