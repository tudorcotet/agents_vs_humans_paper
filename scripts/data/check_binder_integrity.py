"""Assert every shipped complex actually contains the binder it is named after.

Catches the failure mode from issue #5: Modal reuses warm containers across
``.map()`` inputs, so a predictor that writes into a fixed ``/tmp`` dir and then
globs it will happily return the *previous* design's structure. Chain A (the
target) stays correct, so nothing downstream looks wrong until you compare
sequences.

    uv run python scripts/data/check_binder_integrity.py
    uv run python scripts/data/check_binder_integrity.py --root ../rbx1_gem_paper --key pb_id
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import gemmi
import pandas as pd


def read_target(root: Path) -> str:
    """The target chain, so we can tell it apart from a same-length binder."""
    fasta = next((root / "data" / "target").glob("*.fasta"))
    return "".join(
        line.strip() for line in fasta.read_text().splitlines() if not line.startswith(">")
    )


def binder_chains(path: Path, target: str) -> list[str]:
    """One-letter sequence of every chain that isn't the target."""
    model = gemmi.read_structure(str(path))[0]
    seqs = [gemmi.one_letter_code([r.name for r in c]).upper().replace("X", "") for c in model]
    return [s for s in seqs if s != target] or seqs


def audit(root: Path, key: str) -> dict[str, tuple[int, int]]:
    target = read_target(root)
    designs = pd.read_parquet(root / "data" / "designs.parquet")
    # Filenames are either the key verbatim (rbx1: "amber-bat-maple") or the
    # zero-padded design_NNN form, so index both.
    expected: dict[str, str] = {}
    for k, seq in zip(designs[key], designs["sequence"], strict=True):
        expected[str(k)] = seq
        if isinstance(k, (int, float)) or str(k).isdigit():
            expected[f"design_{int(k):03d}"] = seq

    report: dict[str, tuple[int, int]] = {}
    for model_dir in sorted((root / "data" / "structures").iterdir()):
        if not model_dir.is_dir():
            continue
        files = sorted(list(model_dir.glob("*.cif")) + list(model_dir.glob("*.pdb")))
        if not files:
            continue
        ok = bad = 0
        for path in files:
            want = expected.get(path.stem)
            if want is None:
                continue
            # ESMFold and friends drop the initiator Met; that is not corruption.
            got = binder_chains(path, target)
            if any(s == want or want == "M" + s for s in got):
                ok += 1
            else:
                bad += 1
                if bad <= 5:
                    print(f"  {model_dir.name}/{path.name}: binder != {path.stem}")
        report[model_dir.name] = (ok, bad)
    return report


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    ap.add_argument("--key", default="design_id", help="designs.parquet column matching the filename stem")
    args = ap.parse_args()

    report = audit(args.root, args.key)
    failed = False
    for model, (ok, bad) in report.items():
        status = "ok" if bad == 0 else "CORRUPT"
        print(f"{model:12s} {ok:5d} correct  {bad:5d} wrong   {status}")
        failed |= bad > 0
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
