"""
Epitope footprints on TREM2 from the Boltz-2 complex predictions.

For each screened design, the TREM2 (chain A) residues within 5 A of any
binder (chain B) heavy atom define that design's epitope footprint. Footprints
are clustered (Jaccard distance, average linkage) into discrete epitope
patches; the figure reports how many distinct patches each cohort engages.

The TREM2 chain (A) is NOT numbered identically across designs (8 of the
100 carry a longer, frame-shifted construct), so footprints are aligned to a
common reference sequence (the modal TREM2 chain) by pairwise alignment
before indexing — never by raw author residue number.

Inputs : data/structures/boltz2/design_NNN.cif
Output : analyses/helen/epitope_cache.npz
    design_id   : (N,) int
    footprints  : (N, R) bool   — contacted, indexed on the reference TREM2
    trem2_resid : (R,) int      — 1..R position along the reference TREM2

Deterministic, CPU-only, no network.
"""
from __future__ import annotations

import re
from collections import Counter
from pathlib import Path

import gemmi
import numpy as np
from Bio import Align

from scripts.utils import load_designs, repo_root

ROOT = repo_root()
CIF_DIR = ROOT / "data" / "structures" / "boltz2"
CACHE = ROOT / "analyses" / "helen" / "epitope_cache.npz"
CONTACT_A = 5.0  # angstrom, heavy-atom


def _design_id(path: Path) -> int:
    return int(re.search(r"design_(\d+)", path.stem).group(1))


def _trem2_and_binder(model: gemmi.Model) -> tuple[gemmi.Chain, gemmi.Chain]:
    """Chain A is TREM2; the other chain is the binder."""
    trem2 = next(c for c in model if c.name == "A")
    binder = next(c for c in model if c.name != "A")
    return trem2, binder


_ALIGNER = Align.PairwiseAligner(mode="global", match_score=2,
                                 mismatch_score=-1, open_gap_score=-6,
                                 extend_gap_score=-1)


def _design_to_ref(design_seq: str, ref_seq: str) -> dict[int, int]:
    """Map design TREM2 residue order-index -> reference position (matched
    columns only) via global alignment. Robust to the frame-shifted construct."""
    aln = _ALIGNER.align(design_seq, ref_seq)[0]
    da, ra = aln.aligned
    out: dict[int, int] = {}
    for (d0, d1), (r0, _r1) in zip(da, ra, strict=True):
        for k in range(int(d1 - d0)):
            out[int(d0) + k] = int(r0) + k
    return out


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

    # Read each complex once; cache (model, trem2 seq, chains).
    parsed = []
    seqs = []
    for path in cifs:
        m = gemmi.read_structure(str(path))[0]
        trem2, binder = _trem2_and_binder(m)
        seq = gemmi.one_letter_code([r.name for r in trem2]).upper()
        parsed.append((m, trem2, binder, seq))
        seqs.append(seq)

    # Reference = the modal TREM2 sequence (92/100 share one construct);
    # every design is aligned to it so footprints share a common frame.
    ref_seq = Counter(seqs).most_common(1)[0][0]
    foot = np.zeros((len(cifs), len(ref_seq)), dtype=bool)

    for di, (m, trem2, binder, seq) in enumerate(parsed):
        idx_map = _design_to_ref(seq, ref_seq)
        ns = gemmi.NeighborSearch(m, gemmi.UnitCell(), CONTACT_A).populate()
        for j, res in enumerate(trem2):
            rp = idx_map.get(j)
            if rp is None:
                continue
            for atom in res:
                hit = False
                for mk in ns.find_atoms(atom.pos, "\0", radius=CONTACT_A):
                    if mk.to_cra(m).chain.name == binder.name:
                        hit = True
                        break
                if hit:
                    foot[di, rp] = True
                    break

    np.savez(CACHE, design_id=ids, footprints=foot,
             trem2_resid=np.arange(1, len(ref_seq) + 1, dtype=int))
    print(f"[epitope] wrote {CACHE} — {foot.shape}, "
          f"mean footprint size {foot.sum(1).mean():.1f} residues")


if __name__ == "__main__":
    main()
