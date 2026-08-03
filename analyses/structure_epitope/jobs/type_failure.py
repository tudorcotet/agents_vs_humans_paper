"""
Bennett Type 0/I/II failure decomposition.
Type 0 = not expressed. Among expressed non-binders: Type I (folding) = binder conformation in the
complex differs from its free-monomer prediction (high CA-RMSD, ESMFold monomer vs binder-in-complex);
Type II (interface) = folds as predicted but no binding (low RMSD). Uses gemmi (ESMFold CIFs lack
occupancy -> Biopython fails). RMSD by Kabsch on 1:1 CA correspondence (same sequence, same numbering).
"""
import numpy as np, pandas as pd, gemmi
from pathlib import Path
ROOT = Path("/Users/amyhe/Desktop/trem2_2026/agents_vs_humans_paper")
RES = ROOT / "analyses/structure_epitope/results"
TARGET_LEN = 175

def ca_by_resnum(chain):
    out = {}
    for r in chain:
        a = r.find_atom("CA", "*")
        if a is not None: out[r.seqid.num] = np.array([a.pos.x, a.pos.y, a.pos.z])
    return out

def binder_ca_complex(design_id, pred="protenix"):
    st = gemmi.read_structure(str(ROOT / f"data/structures/{pred}/design_{design_id:03d}.cif"))
    m = st[0]
    bnd = next((c for c in m if len([r for r in c]) != TARGET_LEN), None)
    return ca_by_resnum(bnd) if bnd else {}

def monomer_ca(design_id):
    for sub in ["esmfold", "proteintyper"]:
        p = ROOT / f"data/structures/{sub}/design_{design_id:03d}.cif"
        if p.exists():
            st = gemmi.read_structure(str(p)); m = st[0]
            ch = max(m, key=lambda c: len([r for r in c]))
            return ca_by_resnum(ch)
    return {}

def kabsch_rmsd(P, Q):
    Pc = P - P.mean(0); Qc = Q - Q.mean(0)
    V, S, Wt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(V @ Wt))
    R = V @ np.diag([1, 1, d]) @ Wt
    Pr = Pc @ R
    return np.sqrt(((Pr - Qc) ** 2).sum() / len(P))

meta = pd.read_parquet(ROOT / "analyses/structure_epitope/data/designs_canonical.parquet")
rows = []
for did in meta.design_id:
    cc = binder_ca_complex(did); mc = monomer_ca(did)
    common = sorted(set(cc) & set(mc))
    rmsd = np.nan
    if len(common) >= 5:
        P = np.array([cc[i] for i in common]); Q = np.array([mc[i] for i in common])
        rmsd = round(kabsch_rmsd(P, Q), 2)
    rows.append({"design_id": did, "monomer_rmsd": rmsd, "n_common": len(common)})
rr = pd.DataFrame(rows).merge(meta[["design_id","population_tier","expressed","is_hit","cohort"]], on="design_id")

# classify
RMSD_CUT = 3.0  # CA-RMSD threshold for "folds differently in complex"
def failtype(r):
    if r.population_tier == "screened_not_expressed": return "Type0_no_expression"
    if r.is_hit == True: return "Active"
    if r.expressed == True:
        if pd.isna(r.monomer_rmsd): return "expressed_nonbinder_unknown"
        return "TypeI_folding" if r.monomer_rmsd > RMSD_CUT else "TypeII_interface"
    return "non_screened"
rr["failure_class"] = rr.apply(failtype, axis=1)
rr.to_csv(RES / "type_failure.csv", index=False)

print("=== monomer RMSD availability ===")
print("designs with RMSD:", int(rr.monomer_rmsd.notna().sum()), "/141")
print("\n=== failure decomposition (screened) ===")
scr = rr[rr.population_tier != "non_screened"]
print(scr.failure_class.value_counts().to_string())
print("\n=== RMSD: binders vs expressed non-binders ===")
b = rr[rr.is_hit==True].monomer_rmsd.dropna(); nb = rr[(rr.expressed==True)&(rr.is_hit==False)].monomer_rmsd.dropna()
print(f"binder median RMSD {b.median():.2f} Å (n={len(b)}); non-binder {nb.median():.2f} Å (n={len(nb)})")
from scipy import stats
print(f"MWU p={stats.mannwhitneyu(b,nb)[1]:.3f}")
print(f"\nType I (folding, RMSD>{RMSD_CUT}): {int((scr.failure_class=='TypeI_folding').sum())}; "
      f"Type II (interface): {int((scr.failure_class=='TypeII_interface').sum())}")
print("wrote type_failure.csv")
