"""
Phase 5 + Part D: binder fold / topology from DSSP (pydssp, MIT), fold x cohort / outcome / epitope.
Uses protenix binder chains (clean). Interpretable topology labels, not cluster numbers.
"""
import warnings, re, numpy as np, pandas as pd
from pathlib import Path
from Bio.PDB import MMCIFParser
import pydssp, torch
from scipy import stats
warnings.filterwarnings("ignore")
ROOT = Path("/Users/amyhe/Desktop/trem2_2026/agents_vs_humans_paper")
RES = ROOT / "analyses/structure_epitope/results"
TARGET_LEN = 175

def binder_coords(design_id):
    p = ROOT / f"data/structures/protenix/design_{design_id:03d}.cif"
    if not p.exists(): return None
    model = MMCIFParser(QUIET=True).get_structure("m", str(p))[0]
    chs = list(model.get_chains())
    tgt = next((c for c in chs if len([r for r in c if r.id[0]==" "])==TARGET_LEN), None)
    bnd = next(c for c in chs if c is not tgt)
    coords = []
    for r in bnd:
        if r.id[0] != " ": continue
        try: coords.append([r["N"].coord, r["CA"].coord, r["C"].coord, r["O"].coord])
        except KeyError: continue
    return np.array(coords, dtype=np.float32) if coords else None

def ss_string(coords):
    x = torch.tensor(coords)  # (L,4,3)
    ss = pydssp.assign(x, out_type="c3")  # '-'(loop),'H','E'
    return "".join(ss)

def segments(ss, char, minlen):
    return [m.group() for m in re.finditer(f"{char}+", ss) if len(m.group()) >= minlen]

def topology(ss):
    nH = len(segments(ss, "H", 4)); nE = len(segments(ss, "E", 3))
    hf = ss.count("H")/len(ss); ef = ss.count("E")/len(ss)
    if nE >= 2 and ef > 0.15: lab = "beta/mixed"
    elif nE >= 1 and ef > 0.10: lab = "alpha+beta"
    elif nH >= 4: lab = ">=4-helix"
    elif nH == 3: lab = "3-helix bundle"
    elif nH == 2: lab = "2-helix/hairpin"
    elif nH == 1: lab = "single-helix"
    else: lab = "loop/other"
    return lab, nH, nE, round(hf, 3), round(ef, 3)

meta = pd.read_parquet(ROOT / "analyses/structure_epitope/data/designs_canonical.parquet")
rows = []
for did in meta.design_id:
    c = binder_coords(did)
    if c is None or len(c) < 5:
        rows.append({"design_id": did, "topology": "n/a", "helix_frac": np.nan, "strand_frac": np.nan,
                     "n_helix": 0, "n_strand": 0, "binder_len": (0 if c is None else len(c))}); continue
    ss = ss_string(c); lab, nH, nE, hf, ef = topology(ss)
    rows.append({"design_id": did, "topology": lab, "helix_frac": hf, "strand_frac": ef,
                 "n_helix": nH, "n_strand": nE, "binder_len": len(c), "ss": ss})
fold = pd.DataFrame(rows).merge(meta[["design_id","cohort","is_hit","p_screened","p_expressed",
                                      "kd_arith_mean_nM_all","method_family"]], on="design_id")
fold.to_csv(RES / "fold_topology.csv", index=False)

print("=== topology distribution (all 141) ===")
print(fold.topology.value_counts().to_string())
print(f"\nmedian helix_frac = {fold.helix_frac.median():.2f}, median strand_frac = {fold.strand_frac.median():.2f}")

# Part D.1 feasibility: designs with strand fraction > 0.20
beta = fold[fold.strand_frac > 0.20]
print(f"\n=== Part D.1: designs with strand_frac > 0.20 : n = {len(beta)} (feasibility gate ~10) ===")
print(f"  predominantly alpha (strand_frac<0.10): {int((fold.strand_frac<0.10).sum())}/141")

print("\n=== fold x cohort (screened) — agent monoculture? ===")
scr = fold[fold.p_screened==True]
ct = pd.crosstab(scr.topology, scr.cohort)
print(ct.to_string())
# diversity: distinct topologies per cohort (>=2 designs)
for coh in ["human","agent"]:
    occ = scr[scr.cohort==coh].topology.value_counts()
    print(f"  {coh}: {int((occ>=2).sum())} topologies with >=2 designs; {dict(occ)}")

print("\n=== fold x outcome (expressed): hit rate + median Kd per topology ===")
exp = fold[fold.p_expressed==True]
for topo, g in exp.groupby("topology"):
    hits=int((g.is_hit==True).sum()); n=len(g); kd=g.kd_arith_mean_nM_all.dropna()
    print(f"  {topo:18s}: hit {hits}/{n} = {100*hits/n:3.0f}%  Kd med {kd.median() if len(kd) else float('nan'):.0f} nM (n={len(kd)})")

# fold x epitope contingency (using epitope_bins)
try:
    bins = pd.read_csv(RES / "epitope_bins.csv")
    fe = fold.merge(bins[["design_id","bin"]], on="design_id")
    print("\n=== fold x epitope-bin contingency (screened) ===")
    ct2 = pd.crosstab(fe.topology, fe.bin)
    print(ct2.to_string())
    if ct2.shape[0]>1 and ct2.shape[1]>1:
        chi2,p,dof,_=stats.chi2_contingency(ct2)
        v=np.sqrt(chi2/(ct2.values.sum()*(min(ct2.shape)-1)))
        print(f"  chi2={chi2:.1f} p={p:.3f} Cramer's V={v:.3f}  (epitope bins are degenerate: one dominant apical bin)")
except FileNotFoundError:
    pass
print("\nwrote fold_topology.csv")
