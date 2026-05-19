"""
Epitope footprints on TREM2 from the Boltz-2 complex predictions.

For each screened design, the TREM2 (chain A) residues within 5 A of any
binder (chain B) heavy atom define that design's epitope footprint. Footprints
are clustered (Jaccard distance, average linkage) into discrete epitope
patches; the figure reports how many distinct patches each cohort engages.

Inputs : data/proteinbase/boltz2/design_NNN.cif
Output : analyses/helen/epitope_cache.npz
    design_id   : (N,) int
    footprints  : (N, R) bool   — TREM2 residue contacted (R = TREM2 length)
    trem2_resid : (R,) int      — TREM2 author residue numbers (column order)

Deterministic, CPU-only, no network.
"""
from __future__ import annotations

import re
from pathlib import Path

import gemmi
import numpy as np

from scripts.utils import load_designs, repo_root

ROOT = repo_root()
CIF_DIR = ROOT / "data" / "proteinbase" / "boltz2"
CACHE = ROOT / "analyses" / "helen" / "epitope_cache.npz"
CONTACT_A = 5.0  # angstrom, heavy-atom


def _design_id(path: Path) -> int:
    return int(re.search(r"design_(\d+)", path.stem).group(1))


def _trem2_and_binder(model: gemmi.Model) -> tuple[gemmi.Chain, gemmi.Chain]:
    """Chain A is TREM2 (longer, fixed across designs); B is the binder."""
    trem2 = next(c for c in model if c.name == "A")
    binder = next(c for c in model if c.name != "A")
    return trem2, binder


def main() -> None:
    scr = load_designs(only_screened=True)[["design_id"]]
    cifs = sorted(CIF_DIR.glob("design_*.cif"), key=_design_id)
    cifs = [p for p in cifs if _design_id(p) in set(scr.design_id)]
    if not cifs:
        raise FileNotFoundError(f"No Boltz-2 CIFs in {CIF_DIR} (git lfs pull?)")

    ids = np.array([_design_id(p) for p in cifs], dtype=int)
    if CACHE.exists():
        c = np.load(CACHE)
        if set(c["design_id"].tolist()) == set(ids.tolist()):
            print(f"[epitope] cache hit ({CACHE.name}) — skipping")
            return

    # Reference TREM2 residue numbering from the first complex.
    ref = gemmi.read_structure(str(cifs[0]))[0]
    trem2_ref, _ = _trem2_and_binder(ref)
    resids = [r.seqid.num for r in trem2_ref]
    rindex = {n: k for k, n in enumerate(resids)}
    foot = np.zeros((len(cifs), len(resids)), dtype=bool)

    for di, path in enumerate(cifs):
        m = gemmi.read_structure(str(path))[0]
        trem2, binder = _trem2_and_binder(m)
        ns = gemmi.NeighborSearch(m, gemmi.UnitCell(), CONTACT_A).populate()
        for res in trem2:
            k = rindex.get(res.seqid.num)
            if k is None:
                continue
            for atom in res:
                marks = ns.find_atoms(atom.pos, "\0", radius=CONTACT_A)
                for mk in marks:
                    cra = mk.to_cra(m)
                    if cra.chain.name == binder.name:
                        foot[di, k] = True
                        break
                if foot[di, k]:
                    break

    np.savez(CACHE, design_id=ids, footprints=foot,
             trem2_resid=np.array(resids, dtype=int))
    print(f"[epitope] wrote {CACHE} — {foot.shape}, "
          f"mean footprint size {foot.sum(1).mean():.1f} residues")


if __name__ == "__main__":
    main()
