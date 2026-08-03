"""
Phase 3: stalk-violation analysis + curated functional-set overlap.
Reads results/footprints_consensus.parquet (design x construct residue, consensus indicator) and
the canonical table. Consensus footprint = both Protenix & Chai agree (primary contact).

Stalk (competition rule: target IgSF 19-132, NOT stalk 133-172) = construct 116-157 = UniProt
133-174. His-tag/linker (construct 158-175) contacts are prediction artifacts (binder docking the
tag), reported as a QC flag, not an epitope.
"""
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
ROOT = Path("/Users/amyhe/Desktop/trem2_2026/agents_vs_humans_paper")
RES = ROOT / "analyses/structure_epitope/results"
cons = pd.read_parquet(RES / "footprints_consensus.parquet")
meta = pd.read_parquet(ROOT / "analyses/structure_epitope/data/designs_canonical.parquet")[
    ["design_id", "cohort", "is_hit", "expressed", "p_screened", "p_expressed", "team",
     "kd_arith_mean_nM_all"]]

c = cons[cons.consensus == True].copy()
# per-design region counts
reg = c.groupby(["design_id", "region"]).size().unstack(fill_value=0)
for col in ["IgSF", "stalk", "linker", "his_tag", "met"]:
    if col not in reg: reg[col] = 0
reg = reg.reset_index().merge(meta, on="design_id", how="right").fillna({"IgSF":0,"stalk":0,"linker":0,"his_tag":0,"met":0})
reg["any_stalk"] = reg["stalk"] > 0
reg["stalk_dominant"] = reg["stalk"] > reg["IgSF"]
reg["any_histag"] = reg["his_tag"] > 0
reg.to_csv(RES / "epitope_region_counts.csv", index=False)

def rate(sub, col):
    n = len(sub); k = int(sub[col].sum())
    return k, n, (100*k/n if n else 0)

print("################ STALK-VIOLATION ANALYSIS ################")
for popname, mask in [("P_screened", reg.p_screened==True), ("P_expressed", reg.p_expressed==True)]:
    pop = reg[mask]
    print(f"\n--- {popname} (n={len(pop)}) ---")
    for col in ["any_stalk", "stalk_dominant", "any_histag"]:
        k, n, pct = rate(pop, col)
        print(f"  {col:15s}: {k}/{n} = {pct:.1f}%")
    # by cohort
    for col in ["any_stalk", "stalk_dominant"]:
        hu = pop[pop.cohort=="human"]; ag = pop[pop.cohort=="agent"]
        hk,hn,_ = rate(hu,col); ak,an,_ = rate(ag,col)
        OR,p = stats.fisher_exact([[hk,hn-hk],[ak,an-ak]])
        print(f"  {col}: human {hk}/{hn} ({100*hk/hn:.1f}%) vs agent {ak}/{an} ({100*ak/an:.1f}%)  Fisher p={p:.3f} OR={OR:.2f}")

print("\n=== designs contacting the stalk (any_stalk), with detail ===")
sv = reg[reg.any_stalk & (reg.p_screened==True)].sort_values("stalk", ascending=False)
print(sv[["design_id","team","cohort","is_hit","IgSF","stalk","stalk_dominant","any_histag"]].to_string(index=False))

print("\n=== His-tag artifact contacts (screened) ===")
ht = reg[reg.any_histag & (reg.p_screened==True)]
print(f"  {len(ht)} screened designs have binder-His-tag contact (prediction artifact / QC flag)")
if len(ht):
    print(ht[["design_id","team","cohort","is_hit","IgSF","stalk","his_tag"]].to_string(index=False))

# ---- curated functional-set overlap ----
print("\n################ CURATED FUNCTIONAL-SET CONTACT (expressed) ################")
TIP_U=[44,69,70,71,74,89]; BASIC_U=[46,47,62,76,77,87,122,123]
CDR_U=list(range(40,43))+list(range(69,73))+list(range(88,92))  # Perera union CDR1/2/3
sets={"tip":TIP_U,"basic":BASIC_U,"cdr_perera":CDR_U}
byd = c[c.region=="IgSF"].groupby("design_id")["uniprot"].apply(set).to_dict()
rows=[]
for did in meta.design_id:
    s=byd.get(did,set())
    row={"design_id":did}
    for name,us in sets.items():
        row[f"n_{name}"]=len(s & set(us)); row[f"frac_{name}"]=len(s & set(us))/len(us)
    rows.append(row)
ov=pd.DataFrame(rows).merge(meta,on="design_id")
ov.to_csv(RES/"epitope_reference_overlap.csv",index=False)
exp=ov[ov.p_expressed==True]
for name in sets:
    b=exp[exp.is_hit==True][f"frac_{name}"]; nb=exp[exp.is_hit==False][f"frac_{name}"]
    U,p=stats.mannwhitneyu(b,nb,alternative="two-sided")
    print(f"  {name:12s}: binder frac median={b.median():.2f} vs nonbinder={nb.median():.2f}  MWU p={p:.3f}")
print("\nwrote epitope_region_counts.csv, epitope_reference_overlap.csv")
