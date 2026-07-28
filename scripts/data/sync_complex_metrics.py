"""Rewrite data/metrics/{boltz2,af2m}/*.json from the corrected scalars.

PR #6 replaced the corrupt Boltz-2 / AF2-Multimer complexes and the ``b2_*`` /
``af2m_*`` columns in ``designs.{csv,parquet}``, but left the per-design metric
JSONs alone — they still hold the values scored off the mislabeled structures.
``build_grand_metrics.py`` reads those JSONs, so any rebuild of
``grand_metrics.csv`` silently reintroduces the corrupt numbers.

This backfills them from the corrected table PR #6 shipped
(``analyses/structure_epitope/results/recomputed_scalars.csv``) so the canonical
build path produces correct values by construction.

    uv run python scripts/data/sync_complex_metrics.py          # report only
    uv run python scripts/data/sync_complex_metrics.py --apply  # rewrite JSONs
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RECOMPUTED = ROOT / "analyses" / "structure_epitope" / "results" / "recomputed_scalars.csv"

# JSON key -> column suffix. Mirrors the payload the Modal predictors write; the
# build script derives `*_struct_exists` from disk, so it is not stored here.
FIELDS = (
    "iptm", "ptm", "mean_plddt",
    "ipsae_d0res_min", "ipsae_d0res_max",
    "ipsae_d0chn_min", "ipsae_d0chn_max",
    "ipsae_d0dom_min", "ipsae_d0dom_max",
    "iptm_d0chn_min", "iptm_d0chn_max",
    "iptm_af_min", "iptm_af_max",
    "pdockq", "pdockq2", "lis", "n_interface",
)
MODELS = {"boltz2": "b2", "af2m": "af2m"}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true", help="write the JSONs (default: dry run)")
    args = ap.parse_args()

    df = pd.read_csv(RECOMPUTED).set_index("design_id")
    changed = written = 0

    for model, prefix in MODELS.items():
        out_dir = ROOT / "data" / "metrics" / model
        for design_id, row in df.iterrows():
            slug = f"design_{int(design_id):03d}"
            path = out_dir / f"{slug}.json"
            payload = {
                "pb_id": slug,
                "predictor": model,
                "status": row[f"{prefix}_status"],
            }
            for f in FIELDS:
                v = row[f"{prefix}_{f}"]
                payload[f] = int(v) if f == "n_interface" and pd.notna(v) else v
            # n_interface is the only integer field; the rest stay float.

            if path.exists():
                before = json.loads(path.read_text())
                if any(before.get(f) != payload[f] for f in FIELDS):
                    changed += 1
            if args.apply:
                path.write_text(json.dumps(payload))
                written += 1

    print(f"{changed} JSONs differ from the corrected scalars")
    print(f"{'wrote ' + str(written) if args.apply else 'dry run — pass --apply to write'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
