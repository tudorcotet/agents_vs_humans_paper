"""
Phase 6 (pre-registered): metrics vs affinity (pKd) on P_kd. Full statistical treatment per prereg:
Spearman rho, BCa bootstrap 95% CI (5000), permutation p (10000), BH-FDR within tier,
partial correlation controlling for binder length, per-cohort, literature-copy sensitivity.
Seeds fixed for reproducibility.
"""
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
ROOT = Path("/Users/amyhe/Desktop/trem2_2026/agents_vs_humans_paper")
RES = ROOT / "analyses/structure_epitope/results"
RNG = np.random.default_rng(20260725)
NBOOT, NPERM = 5000, 10000

df = pd.read_parquet(ROOT / "data/designs.parquet")
# §5 uses PROTENIX interface descriptors (the section's primary structural model), not the
# Protenix+Chai average, for section-wide consistency with §2 (Amy's call).
desc = pd.read_csv(RES / "interface_descriptors_perpredictor.csv")
desc = desc[desc.predictor == "protenix"].drop(columns=["predictor"]).reset_index(drop=True)
meta = pd.read_parquet(ROOT / "analyses/structure_epitope/data/designs_canonical.parquet")
d = (meta[["design_id","cohort","p_kd","pkd_arith_mean","is_literature_copy"]]
     .merge(df[["design_id","submitted_ipsae","boltz2_iptm","boltz2_plddt","esm_pll_avg",
                "saprot_pll_norm","prodigy_protenix_pkd","prodigy_chai_pkd","sequence_length"]],
            on="design_id")
     .merge(desc[["design_id","bsa","n_iface_binder","n_iface_target","iface_frac_hydrophobic",
                  "salt_bridges","hbonds"]], on="design_id"))
kd = d[d.p_kd == True].copy()
print(f"P_kd n = {len(kd)}")

TIER1 = ["submitted_ipsae","boltz2_iptm","boltz2_plddt","esm_pll_avg","saprot_pll_norm"]
TIER2 = ["bsa","n_iface_binder","n_iface_target","iface_frac_hydrophobic",
         "prodigy_protenix_pkd","prodigy_chai_pkd","salt_bridges","hbonds"]
TIER3 = ["sequence_length"]  # exploratory: length is a strong raw confound (verified rho~0.65)

def bca_ci(x, y, alpha=0.05):
    obs = stats.spearmanr(x, y)[0]
    n = len(x); boots = np.empty(NBOOT)
    for i in range(NBOOT):
        idx = RNG.integers(0, n, n)
        boots[i] = stats.spearmanr(x[idx], y[idx])[0]
    boots = boots[~np.isnan(boots)]
    # bias correction
    z0 = stats.norm.ppf((boots < obs).mean()) if 0 < (boots < obs).mean() < 1 else 0.0
    # acceleration via jackknife
    jack = np.empty(n)
    for i in range(n):
        m = np.ones(n, bool); m[i] = False
        jack[i] = stats.spearmanr(x[m], y[m])[0]
    jbar = jack.mean(); num = ((jbar - jack)**3).sum(); den = 6*(((jbar - jack)**2).sum()**1.5)
    a = num/den if den != 0 else 0.0
    zl, zu = stats.norm.ppf(alpha/2), stats.norm.ppf(1-alpha/2)
    def adj(z): return stats.norm.cdf(z0 + (z0+z)/(1-a*(z0+z)))
    lo, hi = np.quantile(boots, [adj(zl), adj(zu)])
    return obs, lo, hi

def perm_p(x, y):
    obs = abs(stats.spearmanr(x, y)[0]); cnt = 0
    yy = y.copy()
    for _ in range(NPERM):
        RNG.shuffle(yy)
        if abs(stats.spearmanr(x, yy)[0]) >= obs: cnt += 1
    return (cnt+1)/(NPERM+1)

def partial_spearman(x, y, z):  # control for z
    rxy = stats.spearmanr(x, y)[0]; rxz = stats.spearmanr(x, z)[0]; ryz = stats.spearmanr(y, z)[0]
    return (rxy - rxz*ryz)/np.sqrt((1-rxz**2)*(1-ryz**2))

def bh(p):
    p=np.asarray(p,float); n=len(p); o=np.argsort(p)
    r=p[o]*n/(np.arange(n)+1); q=np.minimum.accumulate(r[::-1])[::-1]
    out=np.empty(n); out[o]=np.clip(q,0,1); return out

def run_tier(cols, name):
    rows=[]
    for c in cols:
        sub = kd.dropna(subset=[c,"pkd_arith_mean"])
        x = sub[c].values.astype(float); y = sub.pkd_arith_mean.values.astype(float)
        rho, lo, hi = bca_ci(x, y)
        pp = perm_p(x, y)
        pr = partial_spearman(x, y, sub.sequence_length.values.astype(float))
        # per cohort
        rh = {}
        for coh in ["human","agent"]:
            s2 = sub[sub.cohort==coh]
            rh[coh] = stats.spearmanr(s2[c], s2.pkd_arith_mean)[0] if len(s2)>=6 else np.nan
        # sensitivity: exclude literature copies
        s3 = sub[sub.is_literature_copy != True]
        rho_nolit = stats.spearmanr(s3[c], s3.pkd_arith_mean)[0]
        rows.append({"tier":name,"metric":c,"n":len(sub),"rho":round(rho,3),
                     "ci_lo":round(lo,3),"ci_hi":round(hi,3),"perm_p":round(pp,4),
                     "partial_rho_len":round(pr,3),"rho_human":round(rh["human"],3),
                     "rho_agent":round(rh["agent"],3),"rho_no_litcopy":round(rho_nolit,3)})
    return pd.DataFrame(rows)

r1 = run_tier(TIER1, "Tier1_confidence")
r1["q_bh"] = bh(r1.perm_p.values)
r2 = run_tier(TIER2, "Tier2_structural")
r2["q_bh"] = bh(r2.perm_p.values)
r3 = run_tier(TIER3, "Tier3_exploratory")
r3["q_bh"] = r3.perm_p  # single test
out = pd.concat([r1, r2, r3], ignore_index=True)
out.to_csv(RES / "metric_correlations.csv", index=False)
pd.set_option("display.width", 240); pd.set_option("display.max_columns", 20)
print("\n=== TIER 1 (confidence scores) vs pKd ===")
print(r1.to_string(index=False))
print("\n=== TIER 2 (structural / energetic) vs pKd ===")
print(r2.to_string(index=False))
print("\nCI excludes 0 (|rho| CI both same sign) AND survives BH q<0.10:")
sig = out[((out.ci_lo>0)&(out.ci_hi>0))|((out.ci_lo<0)&(out.ci_hi<0))]
sig = sig[sig.q_bh<0.10]
print(sig[["tier","metric","rho","ci_lo","ci_hi","perm_p","q_bh","partial_rho_len","rho_no_litcopy"]].to_string(index=False))
