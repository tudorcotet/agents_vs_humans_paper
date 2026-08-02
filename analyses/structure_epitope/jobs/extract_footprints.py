"""
Contact-footprint extractor (Phase 3 foundation).

For each design and each CLEAN predictor (Protenix, Chai-1), compute, for every TREM2 target
residue, whether it is contacted by the binder:
  - heavy-atom min distance <= {4.0, 4.5, 5.0} A  (Biopython NeighborSearch / KDTree)
  - dASA > 1 A^2 on complex formation  (Bio.PDB.SASA Shrake-Rupley, BSD; SASA_isolated - SASA_complex)
Primary footprint (prereg): contacted at 4.5 A AND dASA > 1.

Target chain A = 175-aa BLI construct; UniProt = construct_pos + 17. IgSF = 2..115,
stalk = 116..157, linker+His = 158..175 (excluded from epitope, kept only for the min-dist scan
so stalk contacts are detectable).

Output (long): results/footprints_long.parquet  (design_id, predictor, construct_pos, uniprot,
region, min_dist, contact_40/45/50, dasa, iface_dasa, primary).

Run:  uv run python analyses/structure_epitope/jobs/extract_footprints.py [--limit N]
"""
import argparse, os, sys, time, warnings
from pathlib import Path
import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, NeighborSearch, Selection
from Bio.PDB.SASA import ShrakeRupley

warnings.filterwarnings("ignore")
ROOT = Path("/Users/amyhe/Desktop/trem2_2026/agents_vs_humans_paper")
OUT = ROOT / "analyses/structure_epitope/results"
PREDICTORS = ["protenix", "chai"]
TARGET_LEN = 175
CUTOFFS = [4.0, 4.5, 5.0]

def region(cpos):
    if 2 <= cpos <= 115: return "IgSF"
    if 116 <= cpos <= 157: return "stalk"
    if 158 <= cpos <= 165: return "linker"
    return "his_tag" if cpos >= 166 else "met"

def target_binder_chains(model):
    """chain A = 175-aa TREM2 target; the other polymer chain = binder."""
    chains = list(model.get_chains())
    tgt = None; bnd = None
    for ch in chains:
        res = [r for r in ch if r.id[0] == " "]
        if len(res) == TARGET_LEN:
            tgt = ch
        else:
            bnd = ch
    if tgt is None:  # fallback: longest is target
        tgt = max(chains, key=lambda c: len([r for r in c if r.id[0]==" "]))
        bnd = [c for c in chains if c is not tgt][0]
    return tgt, bnd

def heavy_atoms(chain):
    return [a for a in chain.get_atoms() if a.element != "H"]

def process(design_id, predictor, sr):
    p = ROOT / f"data/structures/{predictor}/design_{design_id:03d}.cif"
    if not p.exists(): return None
    st = MMCIFParser(QUIET=True).get_structure("m", str(p))
    model = st[0]
    tgt, bnd = target_binder_chains(model)
    tgt_res = [r for r in tgt if r.id[0] == " "]
    # --- distance contacts ---
    b_heavy = heavy_atoms(bnd)
    ns = NeighborSearch(b_heavy)
    rows = []
    for r in tgt_res:
        cpos = r.id[1]
        mind = np.inf
        for a in r.get_atoms():
            if a.element == "H": continue
            near = ns.search(a.coord, CUTOFFS[-1], level="A")
            if near:
                d = min(np.linalg.norm(a.coord - x.coord) for x in near)
                mind = min(mind, d)
        rows.append({"design_id": design_id, "predictor": predictor,
                     "construct_pos": cpos, "uniprot": cpos + 17, "region": region(cpos),
                     "min_dist": (None if np.isinf(mind) else round(float(mind), 3))})
    dfd = pd.DataFrame(rows)
    for c in CUTOFFS:
        dfd[f"contact_{int(c*10)}"] = dfd["min_dist"].notna() & (dfd["min_dist"] <= c)
    # --- dASA (target side): SASA of target-alone minus SASA of target-in-complex ---
    # complex SASA
    sr.compute(model, level="R")
    complex_sasa = {r.id[1]: r.sasa for r in tgt_res}
    # target alone: detach binder, recompute, then reattach
    bnd_id = bnd.id
    model.detach_child(bnd_id)
    sr.compute(model, level="R")
    alone_sasa = {r.id[1]: r.sasa for r in tgt if r.id[0] == " "}
    model.add(bnd)  # reattach
    dfd["dasa"] = dfd["construct_pos"].map(lambda c: round(alone_sasa.get(c, 0) - complex_sasa.get(c, 0), 3))
    dfd["iface_dasa"] = dfd["dasa"] > 1.0
    dfd["primary"] = dfd["contact_45"] & dfd["iface_dasa"]
    return dfd

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--bench", action="store_true")
    args = ap.parse_args()
    df = pd.read_parquet(ROOT / "data/designs.parquet")
    ids = sorted(df.design_id.tolist())
    if args.limit: ids = ids[:args.limit]
    sr = ShrakeRupley()  # default probe 1.40, n_points 100
    OUT.mkdir(parents=True, exist_ok=True)
    all_rows = []
    t0 = time.time()
    for i, did in enumerate(ids):
        for pred in PREDICTORS:
            t1 = time.time()
            res = process(did, pred, sr)
            if res is not None: all_rows.append(res)
            if args.bench:
                print(f"  design {did} {pred}: {time.time()-t1:.2f}s "
                      f"primary={int(res['primary'].sum())} c45={int(res['contact_45'].sum())}")
        if (i+1) % 20 == 0:
            print(f"...{i+1}/{len(ids)} designs, {time.time()-t0:.1f}s elapsed", file=sys.stderr)
    out = pd.concat(all_rows, ignore_index=True)
    if not args.bench:
        out.to_parquet(OUT / "footprints_long.parquet", index=False)
        print(f"wrote {OUT/'footprints_long.parquet'}  rows={len(out)}  "
              f"designs={out.design_id.nunique()}  time={time.time()-t0:.1f}s")
    else:
        print(f"\nbench total {time.time()-t0:.1f}s for {len(ids)} designs x{len(PREDICTORS)}")

if __name__ == "__main__":
    main()
