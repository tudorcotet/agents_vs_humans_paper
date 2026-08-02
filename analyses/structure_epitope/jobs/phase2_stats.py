"""Phase 2: active-vs-inactive descriptor statistics + binary classifier metrics + pKd correlation.
Reads interface_descriptors.csv + canonical table."""
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
from sklearn.metrics import average_precision_score, roc_auc_score
ROOT = Path("/Users/amyhe/Desktop/trem2_2026/agents_vs_humans_paper")
RES = ROOT / "analyses/structure_epitope/results"
d = pd.read_csv(RES / "interface_descriptors.csv")
meta = pd.read_parquet(ROOT / "analyses/structure_epitope/data/designs_canonical.parquet")
df = d.merge(meta[["design_id","cohort","is_hit","expressed","p_expressed","p_kd",
                   "pkd_arith_mean","binding_strength"]], on="design_id")
DESC = ["bsa","n_iface_binder","n_iface_target","hbonds","salt_bridges",
        "iface_frac_hydrophobic","iface_frac_polar","iface_frac_charged","hbond_density_per100"]

def bh(p):
    p=np.asarray(p,float); n=len(p); o=np.argsort(p)
    r=p[o]*n/(np.arange(n)+1); q=np.minimum.accumulate(r[::-1])[::-1]
    out=np.empty(n); out[o]=np.clip(q,0,1); return out
def cliffs(a,b):
    a=np.asarray(a); b=np.asarray(b); gt=sum((x>b).sum() for x in a); lt=sum((x<b).sum() for x in a)
    return (gt-lt)/(len(a)*len(b))

exp = df[df.p_expressed==True]
b = exp[exp.is_hit==True]; nb = exp[exp.is_hit==False]
print(f"=== ACTIVE vs INACTIVE within P_expressed (binders n={len(b)}, non n={len(nb)}) ===")
rows=[]
for c in DESC:
    x=b[c].dropna(); y=nb[c].dropna()
    U,p=stats.mannwhitneyu(x,y,alternative="two-sided")
    delta=cliffs(x.values,y.values)
    # classifier: does higher descriptor => binder? AP/AUROC on expressed
    sub=exp.dropna(subset=[c]); yv=(sub.is_hit==True).astype(int); sc=sub[c].values
    ap=average_precision_score(yv,sc); auc=roc_auc_score(yv,sc)
    ap2=average_precision_score(yv,-sc); auc2=roc_auc_score(yv,-sc)  # in case inverse
    rows.append({"descriptor":c,"binder_med":round(x.median(),2),"nonbinder_med":round(y.median(),2),
                 "cliffs_delta":round(delta,3),"MWU_p":p,
                 "AP":round(max(ap,ap2),3),"AUROC":round(max(auc,auc2),3),"AUROC_dir":"+" if auc>=auc2 else "-"})
r=pd.DataFrame(rows); r["q_bh"]=bh(r.MWU_p.values)
r=r.sort_values("MWU_p")
print(r.to_string(index=False))
r.to_csv(RES/"phase2_active_vs_inactive.csv",index=False)

# baseline: hit rate among expressed
base=(exp.is_hit==True).mean()
print(f"\nbaseline positive rate (AP null) = {base:.3f}")

# ---- pKd correlation (Tier-2 descriptors), P_kd ----
print(f"\n=== descriptor vs pKd (Spearman), P_kd n={df.p_kd.sum()} ===")
kd=df[df.p_kd==True]
rows=[]
for c in DESC:
    sub=kd.dropna(subset=[c,"pkd_arith_mean"])
    rho,p=stats.spearmanr(sub[c],sub.pkd_arith_mean)
    rows.append({"descriptor":c,"spearman_rho":round(rho,3),"p":round(p,3),"n":len(sub)})
rr=pd.DataFrame(rows).sort_values("spearman_rho",key=abs,ascending=False)
print(rr.to_string(index=False))
rr.to_csv(RES/"phase2_desc_vs_pkd.csv",index=False)

# ---- strong vs weak among binders (Kd median split) ----
kd2=kd.copy(); med=kd2.kd_arith_mean_nM_all.median() if "kd_arith_mean_nM_all" in kd2 else None
print(f"\n(descriptor medians: binders vs non-binders shown above; BSA binder {b.bsa.median():.0f} vs non {nb.bsa.median():.0f})")
