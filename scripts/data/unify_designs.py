"""Merge every grand_metrics column into data/designs.csv so the repo has
ONE wide canonical table instead of two split datasets.

Idempotent: re-running just re-applies the merge.

Outputs:
  data/designs.csv       — 141 rows × ~280 cols (was 141 × 123)
  data/designs.parquet   — typed sibling, regenerated from the CSV.

Per-design JSONs under data/metrics/<scorer>/ stay on disk (in LFS) but
become derivative — their scalar values now live on each design's row.
The PAE matrices under data/metrics/pae/ are NOT scalarised (residue ×
residue payloads) — those remain the source of truth for per-residue
analyses.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from scripts.utils import repo_root


def main() -> None:
    root = repo_root()
    designs_path = root / "data" / "designs.csv"
    grand_path = root / "data" / "grand_metrics.csv"

    designs = pd.read_csv(designs_path)
    grand = pd.read_csv(grand_path)

    # Drop grand-only redundant identity / index columns that designs.csv
    # already owns canonically. `slug` we re-derive from design_id below.
    redundant = [
        "slug",  # = f"design_{design_id:03d}", reconstructable
        "pb_id", "name", "team", "is_human", "sequence_length",
        "is_hit", "submitted_to_lab", "binding_label", "kd_arith_mean_nM_all",
    ]
    overlap = [c for c in redundant if c in grand.columns]
    grand_merge = grand.drop(columns=overlap)

    # Anything in grand still overlapping with designs.csv (e.g. pb_*) — keep
    # the designs.csv version, drop the grand duplicate.
    dup_cols = [c for c in grand_merge.columns if c in designs.columns and c != "design_id"]
    grand_merge = grand_merge.drop(columns=dup_cols)

    print(f"designs.csv:           {designs.shape}")
    print(f"grand_metrics.csv:     {grand.shape}")
    print(f"redundant identity dropped: {len(overlap)}")
    print(f"duplicate column dropped:   {len(dup_cols)}")
    print(f"net new columns to merge:   {len(grand_merge.columns) - 1}")

    unified = designs.merge(grand_merge, on="design_id", how="left")
    print(f"unified:               {unified.shape}")

    # Sanity: no row count drift.
    assert len(unified) == len(designs) == 141, f"row count drift: {len(unified)}"

    # Write back. CSV first (human-readable), parquet second (typed sibling).
    unified.to_csv(designs_path, index=False)
    unified.to_parquet(root / "data" / "designs.parquet", index=False)
    print(f"wrote {designs_path}")
    print(f"wrote {root / 'data' / 'designs.parquet'}")

    # Per-group preview.
    groups: dict[str, int] = {}
    for c in unified.columns:
        for prefix in ("esm_pll_", "netsolp_", "saprot_",
                       "prodigy_boltz2_", "prodigy_protenix_", "prodigy_chai_", "prodigy_af2m_",
                       "destress_boltz2_", "destress_protenix_", "destress_chai_", "destress_af2m_",
                       "b2_", "px_", "chai_", "af2m_", "tp_", "pb_"):
            if c.startswith(prefix):
                key = prefix.rstrip("_")
                groups[key] = groups.get(key, 0) + 1
                break
        else:
            groups.setdefault("designs.csv (curated)", 0)
            groups["designs.csv (curated)"] += 1
    print("\nColumn groups in unified designs.csv:")
    for g, n in sorted(groups.items()):
        print(f"  {g:30s} {n:>3} cols")


if __name__ == "__main__":
    main()
