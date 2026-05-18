"""Canonical data loaders. Every analysis script imports from here.

Why this exists: every analysis used to hard-code a `Path(__file__).parents[N]`
prefix. When the repo layout drifts (someone moves a file, renames a folder),
every script breaks at once. Funnel all loads through this module instead.

Usage:
    from scripts.utils import load_designs, load_replicates
    df = load_designs()                       # the canonical 141-row table
    df = load_designs(only_screened=True)     # the 100 designs sent to the BLI assay
    df = load_designs(only_hits=True)         # the binder subset only
    reps = load_replicates()                  # per-replicate (long) BLI rows
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import pandas as pd


def repo_root() -> Path:
    """Resolve the repo root by walking up to the directory that owns `pyproject.toml`."""
    here = Path(__file__).resolve()
    for parent in [here, *here.parents]:
        if (parent / "pyproject.toml").exists():
            return parent
    raise RuntimeError(f"Could not locate repo root from {here}")


def _designs_path() -> Path:
    root = repo_root()
    parquet = root / "data" / "designs.parquet"
    if parquet.exists():
        return parquet
    return root / "data" / "designs.csv"


# Columns the CSV stores as `True`/`False`/`""` strings — cast on load so callers
# can write `df.is_hit & df.expressed` without pandas screaming about object dtype.
_BOOL_COLUMNS = (
    "is_human",
    "submitted_to_lab",
    "is_hit",
    "expressed",
    "weird_replicates_flag",
    "assay_methods_mixed",
    "modality_blast_overridden",
    "is_literature_copy",
    "is_framework_borrowed",
    "has_g4s_linker",
    "methionine_leader",
)


def _coerce_bools(df: pd.DataFrame) -> pd.DataFrame:
    for col in _BOOL_COLUMNS:
        if col not in df.columns:
            continue
        s = df[col]
        if s.dtype == bool:
            continue
        df[col] = s.map({True: True, False: False, "True": True, "False": False,
                          "true": True, "false": False, 1: True, 0: False}).astype("boolean")
    return df


@lru_cache(maxsize=8)
def load_designs(
    *,
    only_screened: bool = False,
    only_hits: bool = False,
    only_human: bool | None = None,
) -> pd.DataFrame:
    """Load the canonical 141-row design table.

    Bool-typed columns (`is_hit`, `expressed`, `is_human`, …) are cast to pandas
    nullable `boolean` on the way out so you can write `df.is_hit & df.expressed`
    without pandas screaming about object dtype.

    Args:
        only_screened: keep only the 100 designs that went to the BLI assay
                       (drops the 41 that ranked below the ipSAE cutoff).
        only_hits:     keep only designs whose `is_hit=True`. Implies `only_screened`.
        only_human:    True → human cohort only; False → agent cohort only;
                       None (default) → both.
    """
    path = _designs_path()
    if path.suffix == ".parquet":
        df = pd.read_parquet(path)
    else:
        df = pd.read_csv(path)

    df = _coerce_bools(df)

    if only_hits:
        df = df[df["is_hit"].fillna(False)]
    elif only_screened:
        df = df[df["submitted_to_lab"].fillna(False)]

    if only_human is True:
        df = df[df["is_human"].fillna(False)]
    elif only_human is False:
        df = df[~df["is_human"].fillna(False)]

    return df.reset_index(drop=True)


@lru_cache(maxsize=4)
def load_replicates() -> pd.DataFrame:
    """Per-replicate BLI rows. Multiple per design.

    Filter to the curated subset with::

        reps = load_replicates()
        reps = reps[reps.selected & ~reps.excluded & reps.fixed]
    """
    root = repo_root()
    return pd.read_csv(root / "data" / "raw_lab" / "bli_replicates.csv")


@lru_cache(maxsize=4)
def load_bli_results() -> pd.DataFrame:
    """Per-design BLI rows (one per design after replicate selection)."""
    root = repo_root()
    return pd.read_csv(root / "data" / "raw_lab" / "bli_results.csv")


def cohort_palette() -> dict[str, str]:
    """Two-colour palette for human / agent. Loaded from theme/palettes.json so
    every figure uses the same hues."""
    import json

    pal_file = repo_root() / "theme" / "palettes.json"
    if pal_file.exists():
        pal = json.loads(pal_file.read_text())
        if "agent_vs_human" in pal:
            return pal["agent_vs_human"]
    # Fallback (these are the brand defaults; don't change ad hoc).
    return {"human": "#30C5F5", "agent": "#FF6B35"}
