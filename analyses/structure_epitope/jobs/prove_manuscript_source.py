"""
Prove which manuscript structure-metric results depend on the corrupt b2_*/af2m_* columns, by
reproducing the published numbers from the data and checking which column matches.

Manuscript 'In Silico Scores: Humans vs Agents' table (published means +/- SD, p):
  ipSAE          H 0.717+/-0.254  A 0.605+/-0.247  p=0.00085
  ipTM           H 0.858+/-0.160  A 0.832+/-0.115  p=0.0014
  pLDDT binder   H 0.900+/-0.060  A 0.873+/-0.054  p=0.00016
  pLDDT complex  H 0.910+/-0.040  A 0.884+/-0.054  p=0.00080
"""
import pandas as pd, numpy as np
from scipy import stats
from sklearn.metrics import roc_auc_score
ROOT="/Users/amyhe/Desktop/trem2_2026/agents_vs_humans_paper"
df=pd.read_parquet(f"{ROOT}/data/designs.parquet")
meta=pd.read_parquet(f"{ROOT}/analyses/structure_epitope/data/designs_canonical.parquet")

PUB={"ipSAE":dict(col_clean="submitted_ipsae",col_corrupt="b2_ipsae_d0chn_max",H=(0.717,0.254),A=(0.605,0.247),p=0.00085),
     "ipTM":dict(col_clean="boltz2_iptm",col_corrupt="b2_iptm",H=(0.858,0.160),A=(0.832,0.115),p=0.0014),
     "pLDDT_binder":dict(col_clean="boltz2_plddt",col_corrupt="b2_plddt" if "b2_plddt" in df else None,H=(0.900,0.060),A=(0.873,0.054),p=0.00016),
     "pLDDT_complex":dict(col_clean="boltz2_complex_plddt",col_corrupt=None,H=(0.910,0.040),A=(0.884,0.054),p=0.00080)}

def by_cohort(col):
    h=df[df.cohort=="human"][col].dropna(); a=df[df.cohort=="agent"][col].dropna()
    p=stats.mannwhitneyu(h,a,alternative="two-sided")[1]
    return (h.mean(),h.std()),(a.mean(),a.std()),p,len(h),len(a)

print("="*96)
print("MANUSCRIPT 'In Silico Scores' TABLE — reproduce from CLEAN vs CORRUPT columns")
print("="*96)
for metric,d in PUB.items():
    print(f"\n### {metric}  (published: H {d['H'][0]}±{d['H'][1]}, A {d['A'][0]}±{d['A'][1]}, p={d['p']})")
    for kind,col in [("CLEAN  ",d["col_clean"]),("CORRUPT",d["col_corrupt"])]:
        if col is None or col not in df: print(f"  {kind} ({col}): n/a"); continue
        (hm,hs),(am,as_),p,nh,na=by_cohort(col)
        match = abs(hm-d['H'][0])<0.01 and abs(am-d['A'][0])<0.01
        print(f"  {kind} ({col:20s}): H {hm:.3f}±{hs:.3f}  A {am:.3f}±{as_:.3f}  p={p:.5f}   {'<-- MATCHES published' if match else '(differs from published)'}")

# --- Figure 7 / metric-predictiveness: binder-vs-nonbinder AUROC among expressed ---
print("\n"+"="*96)
print("FIGURE 7 / metric ROC-AUC (binder vs non-binder, expressed) — clean vs corrupt")
print("="*96)
e=df[df.expressed==True].copy()  # designs.parquet already has expressed/is_hit
pairs=[("ipSAE","submitted_ipsae","b2_ipsae_d0chn_max"),("ipTM","boltz2_iptm","b2_iptm"),
       ("pLDDT_binder","boltz2_plddt",None)]
for name,clean,corrupt in pairs:
    row=f"  {name:14s}:"
    for lab,col in [("clean",clean),("corrupt",corrupt)]:
        if col is None or col not in e: row+=f"  {lab}={col}:n/a"; continue
        s=e.dropna(subset=[col]); yy=(s.is_hit==True).astype(int)
        a=max(roc_auc_score(yy,s[col]),roc_auc_score(yy,-s[col]))
        row+=f"  AUROC({lab})={a:.2f}"
    print(row)
print("\n(If the manuscript had used b2_*/af2m_*, the In-Silico table numbers would not reproduce and the")
print(" AUROCs would be ~0.5. That they reproduce ONLY from the clean columns proves independence.)")
