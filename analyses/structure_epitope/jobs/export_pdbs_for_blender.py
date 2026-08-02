"""
Export clean, chain-labeled, commonly-oriented PDBs for Blender + Molecular Nodes (Fig 10 gallery).

For every screened design: superpose its TREM2 chain onto a common reference frame (design 17's
target, IgSF CAs), apply the transform to the whole complex, trim the target to the IgSF domain
(construct 1-115 = Met + UniProt 19-132; drops the floppy stalk/linker/His), relabel chains
A = TREM2 target, B = binder, and write one PDB. A single operator camera then gives a consistent
TREM2 orientation across every panel. Output: figures/blender/pdbs/design_NNN.pdb + manifest.csv.

Run: uv run python analyses/structure_epitope/jobs/export_pdbs_for_blender.py
"""
import numpy as np, pandas as pd, gemmi
from pathlib import Path
ROOT = Path("/Users/amyhe/Desktop/trem2_2026/agents_vs_humans_paper")
OUT = ROOT / "analyses/structure_epitope/figures/blender/pdbs"; OUT.mkdir(parents=True, exist_ok=True)
TARGET_LEN = 175
IGSF_LAST = 115        # construct pos: keep target residues 1..115 (Met + IgSF)
REF_ID = 17            # common orientation reference (best binder)
SUP_RANGE = range(2, 116)   # superpose on IgSF CAs (construct 2..115)

def chains_of(st):
    m = st[0]
    tgt = next((c for c in m if len([r for r in c]) == TARGET_LEN), None)
    if tgt is None: tgt = max(m, key=lambda c: len([r for r in c]))
    bnd = next(c for c in m if c.name != tgt.name)
    return m, tgt, bnd

def ca_map(chain):
    out = {}
    for r in chain:
        a = r.find_atom("CA", "*")
        if a is not None: out[r.seqid.num] = np.array([a.pos.x, a.pos.y, a.pos.z])
    return out

def kabsch(P, Q):  # rotate/translate P onto Q
    Pc, Qc = P.mean(0), Q.mean(0)
    H = (P - Pc).T @ (Q - Qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    t = Qc - R @ Pc
    return R, t

# reference frame
refst = gemmi.read_structure(str(ROOT / f"data/structures/protenix/design_{REF_ID:03d}.cif"))
_, rtgt, _ = chains_of(refst)
refca = ca_map(rtgt)
refP = np.array([refca[i] for i in SUP_RANGE if i in refca])

meta = pd.read_parquet(ROOT / "analyses/structure_epitope/data/designs_canonical.parquet")
scr = meta[meta.p_screened == True]
rows = []
for did in scr.design_id:
    p = ROOT / f"data/structures/protenix/design_{did:03d}.cif"
    if not p.exists(): continue
    st = gemmi.read_structure(str(p)); m, tgt, bnd = chains_of(st)
    tca = ca_map(tgt)
    idx = [i for i in SUP_RANGE if i in tca and i in refca]
    P = np.array([tca[i] for i in idx]); Q = np.array([refca[i] for i in idx])
    R, t = kabsch(P, Q)
    # apply transform to every atom
    for ch in m:
        for r in ch:
            for a in r:
                v = np.array([a.pos.x, a.pos.y, a.pos.z]); w = R @ v + t
                a.pos = gemmi.Position(float(w[0]), float(w[1]), float(w[2]))
    # trim target to IgSF (construct 1..115)
    to_del = [r.seqid.num for r in tgt if r.seqid.num > IGSF_LAST]
    for n in to_del:
        for i, r in enumerate(tgt):
            if r.seqid.num == n: del tgt[i]; break
    # relabel chains: target -> A, binder -> B
    tgt.name = "A"; bnd.name = "B"
    st.setup_entities()
    out = OUT / f"design_{did:03d}.pdb"
    st.write_pdb(str(out))
    row = scr[scr.design_id == did].iloc[0]
    rows.append({"design_id": did, "pdb": f"pdbs/design_{did:03d}.pdb", "name": row["name"],
                 "leaderboard_rank": row.get("selection_rank"), "submitted_ipsae": row.get("submitted_ipsae"),
                 "cohort": row.cohort, "team": row.team, "is_hit": bool(row.is_hit == True),
                 "kd_nM": row.kd_arith_mean_nM_all, "binder_chain_len": len([r for r in bnd])})
man = pd.DataFrame(rows)
# add topology if available
try:
    fold = pd.read_csv(ROOT / "analyses/structure_epitope/results/fold_topology.csv")[["design_id","topology"]]
    man = man.merge(fold, on="design_id", how="left")
except FileNotFoundError:
    pass
man = man.sort_values("leaderboard_rank", na_position="last")  # competition leaderboard order (Fig 1)
man.to_csv(OUT.parent / "manifest.csv", index=False)
print(f"wrote {len(man)} PDBs -> {OUT}")
print("target chain A trimmed to IgSF (construct 1-115 = UniProt 19-132); binder = chain B; common frame = design 17")
print("\nsuggested gallery picks (top binders):")
print(man[man.is_hit].head(6)[["design_id","name","cohort","kd_nM","topology"]].to_string(index=False))
