"""
Part D / manuscript line 295: does Boltz-2 ipSAE reward helicity beyond its coupling to binding?
H3: compare standardized association of binder helix fraction with ipSAE vs. with binding outcome.
If helix predicts ipSAE more strongly than it predicts binding, the ranking metric rewards a
structural property only loosely coupled to binding. Helix fraction from DSSP (pydssp).
"""
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
from sklearn.linear_model import LogisticRegression, LinearRegression
ROOT = Path("/Users/amyhe/Desktop/trem2_2026/agents_vs_humans_paper")
RES = ROOT / "analyses/structure_epitope/results"
fold = pd.read_csv(RES / "fold_topology.csv")  # has cohort,is_hit,p_screened,p_expressed,method_family
df = pd.read_parquet(ROOT / "data/designs.parquet")[["design_id","submitted_ipsae"]]
d = fold.merge(df, on="design_id")

print("=== D.1 feasibility recap ===")
print(f"designs with strand_frac>0.20: {int((d.strand_frac>0.20).sum())}/141;  "
      f"predominantly alpha (helix_frac>0.6): {int((d.helix_frac>0.6).sum())}/141")

def z(x):
    x = np.asarray(x, float); return (x - np.nanmean(x)) / np.nanstd(x)

# H3a: helix_frac -> ipSAE (all with ipsae)
a = d.dropna(subset=["submitted_ipsae","helix_frac"])
rho_i, p_i = stats.spearmanr(a.helix_frac, a.submitted_ipsae)
lin = LinearRegression().fit(z(a.helix_frac).reshape(-1,1), z(a.submitted_ipsae))
print(f"\nH3a  helix_frac vs ipSAE (n={len(a)}): Spearman rho={rho_i:.3f} (p={p_i:.4f}); "
      f"standardized beta={lin.coef_[0]:.3f}")

# H3b: helix_frac -> binding (expressed)
e = d[d.p_expressed==True].dropna(subset=["helix_frac"]).copy()
e["y"] = (e.is_hit==True).astype(int)
rho_b, p_b = stats.spearmanr(e.helix_frac, e.y)
log = LogisticRegression().fit(z(e.helix_frac).reshape(-1,1), e.y)
print(f"H3b  helix_frac vs binding (n={len(e)}): point-biserial/Spearman rho={rho_b:.3f} (p={p_b:.4f}); "
      f"standardized logit beta={log.coef_[0][0]:.3f}")

print("\n=== H3 verdict ===")
print(f"|assoc with ipSAE| = {abs(rho_i):.3f}  vs  |assoc with binding| = {abs(rho_b):.3f}")
if abs(rho_i) > abs(rho_b) + 0.1:
    print("-> ipSAE tracks helicity MORE than binding does: supports 'metric rewards helicity beyond "
          "predictive value' (line 295), on our data.")
elif abs(rho_b) > abs(rho_i) + 0.1:
    print("-> binding tracks helicity more than ipSAE: line-295 claim NOT supported.")
else:
    print("-> comparable/weak: line-295 claim underpowered; report descriptively, do not assert.")

# within-PXDesign stratum (largest, spans both cohorts) to de-confound tool
px = d[d.method_family=="PXDesign"].dropna(subset=["submitted_ipsae","helix_frac"])
if len(px)>=10:
    r,pp = stats.spearmanr(px.helix_frac, px.submitted_ipsae)
    print(f"\nwithin PXDesign (n={len(px)}): helix_frac vs ipSAE rho={r:.3f} (p={pp:.3f})")

out = pd.DataFrame([{"test":"helix_vs_ipsae","n":len(a),"rho":round(rho_i,3),"p":round(p_i,4),"std_beta":round(lin.coef_[0],3)},
                    {"test":"helix_vs_binding","n":len(e),"rho":round(rho_b,3),"p":round(p_b,4),"std_beta":round(log.coef_[0][0],3)}])
out.to_csv(RES/"helical_bias_H3.csv", index=False)
print("\nwrote helical_bias_H3.csv")
