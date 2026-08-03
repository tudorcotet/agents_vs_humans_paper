"""Robustness of the H-Perera / basic-patch enrichment across footprint definition and cutoff.
Also labels residue identities. Reads footprints_long.parquet + canonical table."""
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
ROOT = Path("/Users/amyhe/Desktop/trem2_2026/agents_vs_humans_paper")
RES = ROOT / "analyses/structure_epitope/results"
fp = pd.read_parquet(RES / "footprints_long.parquet")
meta = pd.read_parquet(ROOT / "analyses/structure_epitope/data/designs_canonical.parquet")

CONSTRUCT = "MHNTTVFQGVAGQSLQVSCPYDSMKHWGRRKAWCRQLGEKGPCQRVVSTHNLWLLSFLRRWNGSTAITDDTLGGTLTITLRNLQPHDAGLYQCQSLHGSEADTLRKVLVEVLADPLDHRDAGDLWFPGESESFEDAHVEHSISRSLLEGEIPFPPTSGGGSGGGSHHHHHHHHHH"
def aa(uniprot):  # residue letter at a UniProt position
    cpos = uniprot - 17
    return CONSTRUCT[cpos - 1] if 1 <= cpos <= len(CONSTRUCT) else "?"

def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    r = p[o] * n / (np.arange(n) + 1); q = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty(n); out[o] = np.clip(q, 0, 1); return out

TIP_U = [44, 69, 70, 71, 74, 89]; BASIC_U = [46, 47, 62, 76, 77, 87, 122, 123]

exp = meta[meta.p_expressed == True][["design_id", "is_hit"]].copy()
exp["binder"] = exp.is_hit == True
exp_ids = set(exp.design_id); binder_ids = set(exp[exp.binder].design_id)
nonb_ids = set(exp[~exp.binder].design_id)

print("=== residue identities of key enriched positions ===")
for u in [67, 75, 72, 76, 74, 69, 71]:
    print(f"  UniProt {u} = {aa(u)}  (construct {u-17})")

def footprint(defn, cutoff):
    """return dict design_id -> set(construct_pos) under a footprint definition."""
    col = f"contact_{int(cutoff*10)}"
    d = fp.copy()
    d["hit"] = d[col] & (d["dasa"] > 1.0)  # keep dASA gate
    piv = d.pivot_table(index=["design_id", "construct_pos"], columns="predictor",
                        values="hit", aggfunc="first").reset_index()
    for p in ["protenix", "chai"]:
        if p not in piv: piv[p] = False
    piv[["protenix", "chai"]] = piv[["protenix", "chai"]].fillna(False)
    if defn == "consensus": piv["f"] = piv.protenix & piv.chai
    elif defn == "union": piv["f"] = piv.protenix | piv.chai
    else: piv["f"] = piv[defn]
    return piv[piv.f].groupby("design_id")["construct_pos"].apply(set).to_dict()

def agg_test(fpr, resid_u):
    resid_c = set(u - 17 for u in resid_u)
    flag = {d: bool(resid_c & fpr.get(d, set())) for d in exp_ids}
    a = sum(flag[d] and d in binder_ids for d in exp_ids)
    b = sum(flag[d] and d in nonb_ids for d in exp_ids)
    c = len(binder_ids) - a; e = len(nonb_ids) - b
    OR, p = stats.fisher_exact([[a, b], [c, e]])
    return a, b, OR, p

print("\n=== aggregate tests across footprint defn x cutoff (within P_expressed) ===")
print(f"{'defn':10s} {'cut':4s} | {'TIP a/b OR p':28s} | {'BASIC a/b OR p'}")
for defn in ["consensus", "union", "protenix", "chai"]:
    for cut in [4.0, 4.5, 5.0]:
        fpr = footprint(defn, cut)
        ta, tb, tOR, tp = agg_test(fpr, TIP_U)
        ba, bb, bOR, bp = agg_test(fpr, BASIC_U)
        print(f"{defn:10s} {cut:<4} | tip {ta:2d}/{tb:2d} OR={tOR:5.2f} p={tp:.3f}     | "
              f"basic {ba:2d}/{bb:2d} OR={bOR:5.2f} p={bp:.4f}")

print("\n=== per-residue FDR survivors (q<0.05) across footprint defn (cutoff 4.5) ===")
igsf_c = list(range(2, 116))
for defn in ["consensus", "union", "protenix", "chai"]:
    fpr = footprint(defn, 4.5)
    hit_by_res = {}
    for d, s in fpr.items():
        for c in s: hit_by_res.setdefault(c, set()).add(d)
    rows = []
    for c in igsf_c:
        h = hit_by_res.get(c, set())
        a = len(h & binder_ids); bb = len(binder_ids) - a
        b = len(h & nonb_ids); dd = len(nonb_ids) - b
        OR, p = stats.fisher_exact([[a, b], [bb, dd]])
        rows.append((c + 17, aa(c + 17), a, b, OR, p))
    R = pd.DataFrame(rows, columns=["uniprot", "aa", "bind", "nonb", "OR", "p"])
    R["q"] = bh(R.p.values)
    sig = R[R.q < 0.05].sort_values("p")
    labels = ", ".join(f"{r.aa}{r.uniprot}(q={r.q:.3f})" for _, r in sig.iterrows())
    print(f"  {defn:10s}: {labels if labels else 'none'}")
