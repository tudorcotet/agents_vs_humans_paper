"""
All-vs-all structural clustering of the 100 screened binders via Foldseek
TM-score.

Inputs : data/structures/esmfold/design_NNN.cif (single-chain binder folds)
Binary : .tools/foldseek/bin/foldseek (local static binary)
Output : analyses/helen/foldseek_cache.npz
    design_id : (N,) int       — order of the TM matrix
    tm        : (N, N) float32 — symmetric similarity: Foldseek `alntmscore`
                (TM normalised by ALIGNMENT length, not query/target length),
                symmetrised as max of the two directions, diagonal forced 1.0.

Deterministic, CPU-only, no network — consistent with the analyses contract.
Re-run is a no-op if the cache covers the same design set.
"""
from __future__ import annotations

import re
import subprocess
import tempfile
from pathlib import Path

import gemmi
import numpy as np

from scripts.utils import repo_root

ROOT = repo_root()
ESM_DIR = ROOT / "data" / "structures" / "esmfold"
FOLDSEEK = ROOT / ".tools" / "foldseek" / "bin" / "foldseek"
CACHE = ROOT / "analyses" / "helen" / "foldseek_cache.npz"


def _design_id(path: Path) -> int:
    return int(re.search(r"design_(\d+)", path.stem).group(1))


def main() -> None:
    cifs = sorted(ESM_DIR.glob("design_*.cif"), key=_design_id)
    ids = np.array([_design_id(p) for p in cifs], dtype=int)
    if not cifs:
        raise FileNotFoundError(f"No ESMFold CIFs in {ESM_DIR} (git lfs pull?)")
    if not FOLDSEEK.exists():
        raise FileNotFoundError(f"Foldseek binary missing at {FOLDSEEK}")
    if CACHE.exists():
        c = np.load(CACHE)
        if set(c["design_id"].tolist()) == set(ids.tolist()):
            print(f"[foldseek] cache hit ({CACHE.name}) — skipping")
            return

    idx = {int(i): k for k, i in enumerate(ids)}
    tm = np.zeros((len(ids), len(ids)), dtype=np.float32)
    np.fill_diagonal(tm, 1.0)

    with tempfile.TemporaryDirectory() as td:
        tdp = Path(td)
        # ESMFold CIFs are atom-only mmCIF that Foldseek's parser rejects;
        # round-trip through gemmi to PDB first.
        pdb_dir = tdp / "pdb"
        pdb_dir.mkdir()
        for p in cifs:
            st = gemmi.read_structure(str(p))
            st.setup_entities()
            st.write_pdb(str(pdb_dir / f"{p.stem}.pdb"))
        aln = tdp / "aln.tsv"
        # TMalign mode (--alignment-type 1), exhaustive so every pair is scored.
        subprocess.run(
            [str(FOLDSEEK), "easy-search", str(pdb_dir), str(pdb_dir),
             str(aln), str(tdp / "tmp"),
             "--alignment-type", "1", "--exhaustive-search", "1",
             "--tmscore-threshold", "0.0", "-e", "1000000",
             "--format-output", "query,target,alntmscore"],
            check=True, capture_output=True, text=True,
        )
        for line in aln.read_text().splitlines():
            q, t, score = line.split("\t")
            qi, ti = idx.get(_design_id(Path(q))), idx.get(_design_id(Path(t)))
            if qi is None or ti is None:
                continue
            s = float(score)
            tm[qi, ti] = max(tm[qi, ti], s)
            tm[ti, qi] = max(tm[ti, qi], s)

    # Foldseek self-hits report alntmscore slightly >1 (alignment-length
    # normalisation); force an exact unit diagonal.
    np.fill_diagonal(tm, 1.0)
    np.savez(CACHE, design_id=ids, tm=tm)
    print(f"[foldseek] wrote {CACHE} — TM matrix {tm.shape}, "
          f"median off-diag {np.median(tm[~np.eye(len(ids), dtype=bool)]):.3f}")


if __name__ == "__main__":
    main()
