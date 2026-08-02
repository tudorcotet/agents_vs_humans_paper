"""
Phase 3: epitope binning. Jaccard distance on consensus IgSF footprints -> hierarchical clustering,
silhouette-selected cut. Supersedes the manuscript's '5 human / 3 agent patches (Jaccard 0.7)'.
Report bins, cohort membership, hit rate (Wilson CI), Kd per bin.
"""
import numpy as np, pandas as pd
from pathlib import Path
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from sklearn.metrics import silhouette_score
ROOT = Path("/Users/amyhe/Desktop/trem2_2026/agents_vs_humans_paper")
RES = ROOT / "analyses/structure_epitope/results"
cons = pd.read_parquet(RES / "footprints_consensus.parquet")
meta = pd.read_parquet(ROOT / "analyses/structure_epitope/data/designs_canonical.parquet")

def wilson(k, n, z=1.96):
    if n == 0: return (np.nan, np.nan)
    p = k/n; d = 1+z*z/n
    c = (p+z*z/(2*n))/d; h = z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return (round(100*(c-h),1), round(100*(c+h),1))

# binary matrix: design x IgSF residue (consensus)
igsf = cons[(cons.region == "IgSF")]
mat = igsf.pivot_table(index="design_id", columns="uniprot", values="consensus", aggfunc="first").fillna(False).astype(int)
# restrict to SCREENED designs (compare with the manuscript's screened-set patch count)
scr_ids = set(meta[meta.p_screened == True].design_id)
mat = mat.loc[mat.index.isin(scr_ids)]
# drop designs with empty footprint (can't cluster a zero vector meaningfully)
nonzero = mat.sum(axis=1) > 0
mat = mat[nonzero]
print(f"clustering {len(mat)} screened designs with non-empty footprint (of {len(scr_ids)} screened)")

# jaccard distance
from scipy.spatial.distance import pdist
D = pdist(mat.values, metric="jaccard")
Z = linkage(D, method="average")

# silhouette sweep
best_k, best_s = None, -1
sweep = []
Dsq = squareform(D)
for k in range(2, 9):
    labels = fcluster(Z, k, criterion="maxclust")
    if len(set(labels)) < 2: continue
    s = silhouette_score(Dsq, labels, metric="precomputed")
    sweep.append((k, round(s, 3)))
    if s > best_s: best_s, best_k = s, k
print("silhouette by k:", sweep)
print(f"chosen k = {best_k} (silhouette {best_s:.3f})")

labels = fcluster(Z, best_k, criterion="maxclust")
mat_meta = meta.set_index("design_id").loc[mat.index]
res = pd.DataFrame({"design_id": mat.index, "bin": labels,
                    "cohort": mat_meta.cohort.values, "is_hit": mat_meta.is_hit.values,
                    "kd": mat_meta.kd_arith_mean_nM_all.values})
res.to_csv(RES / "epitope_bins.csv", index=False)

print("\n=== epitope bins (screened, silhouette-selected) ===")
for b in sorted(set(labels)):
    g = res[res.bin == b]
    nh = int((g.cohort == "human").sum()); na = int((g.cohort == "agent").sum())
    hits = int((g.is_hit == True).sum()); n = len(g)
    lo, hi = wilson(hits, n)
    kd = g.kd.dropna()
    # centroid: residues contacted by >50% of bin
    sub = mat.loc[g.design_id]; cent = sub.mean(axis=0); top = list(cent[cent > 0.5].index)
    print(f"  bin {b}: n={n} (human {nh}, agent {na}) | hit {hits}/{n} = {100*hits/n:.0f}% (Wilson {lo}-{hi}%) | "
          f"Kd med {kd.median():.0f} nM (n={len(kd)})")
    print(f"         core residues (>50% of bin): {top}")

# cohort distinct-patch count analog (how many bins does each cohort occupy with >=2 designs?)
print("\n=== bins occupied per cohort (>=2 designs) — supersedes '5 human / 3 agent patches' ===")
for coh in ["human", "agent"]:
    occ = res[res.cohort == coh].groupby("bin").size()
    print(f"  {coh}: occupies {int((occ>=2).sum())} bins (>=2 designs); distribution {dict(occ)}")
print("\nwrote epitope_bins.csv")
