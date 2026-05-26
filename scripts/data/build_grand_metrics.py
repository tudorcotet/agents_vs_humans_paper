"""Pool per-model metric JSONs into a single wide ``grand_metrics.csv``.

Phase 3 of the rerun. After ``run_complex.py`` has filled
``data/metrics/{boltz2,protenix,chai,af2m}/`` and
``data/structures/{…}/`` with one file per design, this script joins
everything together and writes ``data/grand_metrics.csv``.

Column-prefix convention (also documented in ``docs/DATA.md``):

* ``pb_*``    — **ProteinBase mirror**, sourced from ``designs.csv``.
                100/141 coverage (screened designs only). Superset of
                monomer ProteinTyper fields PLUS ProteinBase's own
                Boltz-2 complex run, interface residue counts, and
                wet-lab roll-ups (``pb_n_bli_curves``, ``pb_kd_M_*``).
                These are the columns ProteinBase ships verbatim.
* ``tp_*``    — **Unified ProteinTyper monomer panel.** Same upstream
                tool as the typer subset of ``pb_*``, reshaped to one
                column per monomer metric covering all 141 designs:
                  - 41 non-screened → local rerun JSONs at
                    ``data/metrics/proteintyper/design_NNN.json``
                    (new ``sequence.metrics[]`` schema, 5 fields).
                  - 100 screened → fall back to ``pb_<field>``.
                Use ``tp_*`` when you want a single column per metric
                across all 141 designs; use ``pb_*`` when provenance
                matters or you need the wet-lab / complex / interface
                extras that aren't in ``tp_*``.
* ``b2_*``    — Boltz-2 complex (chain A=target, B=binder), our rerun.
* ``px_*``    — Protenix-v2 complex, our rerun.
* ``chai_*``  — Chai-1 complex, our rerun.
* ``af2m_*``  — AlphaFold2-Multimer (ColabFold) complex, our rerun.

Derived consensus columns at the end:

* ``ipsae_pass_4folders``  — # of {b2, px, chai, af2m} with d0chn_max >= 0.4
* ``iptm_pass_4folders``   — # of {b2, px, chai, af2m} with iptm >= 0.7

The script is idempotent. Missing JSONs leave nulls; missing CIFs leave
``<prefix>_struct_exists`` = False. Re-run after pulling new data.

Usage::

    mise run build:grand
    # equivalently:
    uv run python scripts/data/build_grand_metrics.py
"""

from __future__ import annotations

import argparse
import json
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import pandas as pd

from scripts.utils.load_data import repo_root

# Models we expect under data/metrics/<model>/ and data/structures/<model>/.
# `proteintyper` produces its own folder but those columns are already on
# designs.csv via build_designs.py — we re-export them with `tp_` here so
# `grand_metrics.csv` is fully self-contained.
COMPLEX_MODELS: tuple[str, ...] = ("boltz2", "protenix", "chai", "af2m")
PREFIX: dict[str, str] = {
    "boltz2": "b2",
    "protenix": "px",
    "chai": "chai",
    "af2m": "af2m",
    "proteintyper": "tp",
}
# AF2-M writes PDB; the other three write CIF.
STRUCTURE_EXT: dict[str, str] = {
    "boltz2": "cif",
    "protenix": "cif",
    "chai": "cif",
    "af2m": "pdb",
}

# Fields we lift off the typer monomer panel onto the grand row. These are
# the canonical ProteinBase column names (matching `pb_*` on designs.csv).
TYPER_FIELDS: tuple[str, ...] = (
    "esmfold_plddt",
    "proteinmpnn_score",
    "proteinmpnn_seq_recovery",
    "redesigned_proteinmpnn_score",
    "molecular_weight",
    "isoelectric_point",
    "novelty",
    "seqidentity",
    "seqidentity_afdb50",
    "evalue_afdb50",
    "tm_score_afdb50",
    "ted_confidence",
    "design_class",
    "classification",
    "foldstring",
)

# Fields we lift off each complex predictor's JSON. Names mirror the per-
# model output of `compute_ipsae` plus the native confidence scalars.
COMPLEX_FIELDS_BASE: tuple[str, ...] = (
    "iptm",
    "ptm",
    "mean_plddt",
    "ipsae_d0res_min",
    "ipsae_d0res_max",
    "ipsae_d0chn_min",
    "ipsae_d0chn_max",
    "ipsae_d0dom_min",
    "ipsae_d0dom_max",
    "iptm_d0chn_min",
    "iptm_d0chn_max",
    "iptm_af_min",
    "iptm_af_max",
    "pdockq",
    "pdockq2",
    "lis",
    "n_interface",
)

# Predictor-specific extra fields.
EXTRA_FIELDS: dict[str, tuple[str, ...]] = {
    "boltz2": (),
    "protenix": ("ranking_score", "model_name"),
    "chai": ("aggregate_score",),
    "af2m": (),
}

# ---------------------------------------------------------------------------
# Post-folding scorers — sequence-only + per-complex
# ---------------------------------------------------------------------------

# Sequence-only scorers write `data/metrics/<scorer>/<slug>.json` with the
# fields below. Column names in grand_metrics.csv are prefixed with the
# scorer name (e.g. `esm_pll_total`, `netsolp_solubility`, `saprot_pll_norm`).
SEQUENCE_SCORERS: dict[str, tuple[str, ...]] = {
    "esm_pll": ("esm_pll_total", "esm_pll_avg", "length"),
    "netsolp": ("netsolp_solubility", "netsolp_usability"),
    "saprot": ("saprot_pll", "saprot_pll_norm", "length"),
}

# Per-complex-model scorers write `data/metrics/<scorer>_<model>/<slug>.json`.
# Column names are emitted as `<scorer>_<model>_<field>` for every (scorer,
# model, field) triple — e.g. `prodigy_boltz2_pkd`, `destress_chai_gravy`.
PER_MODEL_SCORERS: dict[str, tuple[str, ...]] = {
    "prodigy": (
        "prodigy_kd",
        "prodigy_pkd",
        "prodigy_dg",
        "prodigy_temperature",
    ),
    "destress": (
        "rosetta_total_per_aa",
        "rosetta_fa_atr_per_aa",
        "rosetta_fa_rep_per_aa",
        "rosetta_fa_sol_per_aa",
        "rosetta_fa_elec_per_aa",
        "rosetta_hbond_sr_bb_per_aa",
        "rosetta_hbond_lr_bb_per_aa",
        "rosetta_hbond_bb_sc_per_aa",
        "rosetta_hbond_sc_per_aa",
        "rosetta_rama_prepro_per_aa",
        "rosetta_fa_dun_per_aa",
        "evoef2_total_per_aa",
        "budeff_total_per_aa",
        "budeff_steric_per_aa",
        "budeff_desolvation_per_aa",
        "budeff_charge_per_aa",
        "budeff_interaction_total",
        "budeff_interaction_steric",
        "budeff_interaction_desolvation",
        "budeff_interaction_charge",
        "isoelectric_point",
        "molecular_weight",
        "num_residues",
        "instability_index",
        "gravy",
    ),
}


def _read_json(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text())
    except Exception:
        return None


def _flatten(
    slug: str, model: str, data: dict[str, Any] | None, fields: Iterable[str]
) -> dict[str, Any]:
    """Return ``{prefix_field: value}`` for one (slug, model) row."""
    prefix = PREFIX[model]
    out: dict[str, Any] = {}
    out[f"{prefix}_status"] = (data or {}).get("status")
    for f in fields:
        out[f"{prefix}_{f}"] = (data or {}).get(f)
    return out


def _extract_typer_metric(data: dict[str, Any], name: str) -> Any:
    """Pull a metric named ``name`` out of the new-schema ProteinTyper JSON.

    Schema::

        {"sequence": {"metrics": [{"metric_type": {"name": "design_class"},
                                    "value": {"slug": "peptide",
                                              "human_readable": "Peptide"}},
                                   {"metric_type": {"name": "molecular_weight"},
                                    "value": {"value": 8439.32}},
                                   ...]}}

    Returns the unwrapped scalar (the slug, or the float, or whatever the
    inner ``value`` field exposes). Falls back to the legacy flat shape
    where ``data[name]`` already holds the scalar.
    """
    flat = data.get(name)
    if flat is not None:
        return flat
    metrics = (data.get("sequence") or {}).get("metrics") or []
    for m in metrics:
        mt = m.get("metric_type") or {}
        if mt.get("name") == name:
            v = m.get("value")
            if isinstance(v, dict):
                # design_class → {"slug": ..., "human_readable": ...}
                if "slug" in v:
                    return v.get("slug")
                # numeric metrics → {"value": ...}
                if "value" in v:
                    return v.get("value")
                return v
            return v
    return None


# Map our internal column suffix → upstream metric name (new ProteinTyper
# schema). The five names below are the only fields the current API still
# exposes; the other nine legacy TYPER_FIELDS entries (esmfold_plddt,
# proteinmpnn_*, ted_confidence, foldstring, ...) live on designs.csv as
# pb_* for the 100 mirrored designs but are simply not produced by the new
# typer endpoint for the 41 rerun designs.
_TYPER_NEW_SCHEMA_NAMES: dict[str, str] = {
    "design_class": "design_class",
    "molecular_weight": "molecular_weight",
    "isoelectric_point": "isoelectric_point",
    "seqidentity": "seqIdentity",  # upstream uses camelCase
    "novelty": "novelty",
}


def _collect_typer(slug: str, designs_row: pd.Series | None = None) -> dict[str, Any]:
    """Unified typer monomer panel for one design.

    Sources, in order of preference:
      1. Local JSON at ``data/metrics/proteintyper/<slug>.json``. This is
         the new-schema payload from the 41-design rerun (sequence.metrics[]
         shape). 5 fields populated: design_class, molecular_weight,
         isoelectric_point, seqIdentity (camelCase upstream), novelty.
      2. The matching ``pb_*`` column on ``designs.csv`` for whichever
         tp_* fields the new schema doesn't expose. 100/141 coverage —
         this is where esmfold_plddt, proteinmpnn_*, ted_confidence,
         foldstring, classification, seqidentity_afdb50, tm_score_afdb50
         come from.

    Net effect: every tp_* column is populated for any design that has
    EITHER the local rerun JSON OR a ProteinBase mirror entry. Only the
    handful of fields that exist in neither remain null.
    """
    typer_path = repo_root() / "data" / "metrics" / "proteintyper" / f"{slug}.json"
    data = _read_json(typer_path)

    flat: dict[str, Any] = {}
    for f in TYPER_FIELDS:
        flat[f"tp_{f}"] = None

    # 1. New-schema rerun JSON.
    if data:
        for col, upstream in _TYPER_NEW_SCHEMA_NAMES.items():
            flat[f"tp_{col}"] = _extract_typer_metric(data, upstream)
        # Legacy flat fields, if the JSON happens to have them.
        for f in TYPER_FIELDS:
            if flat.get(f"tp_{f}") is None and f in data:
                flat[f"tp_{f}"] = data.get(f)

    # 2. ProteinBase mirror via pb_* on designs.csv (covers the 100
    #    screened designs that don't have a local rerun JSON).
    if designs_row is not None:
        for f in TYPER_FIELDS:
            if flat[f"tp_{f}"] is None:
                pb_key = f"pb_{f}"
                if pb_key in designs_row.index:
                    v = designs_row[pb_key]
                    if pd.notna(v):
                        flat[f"tp_{f}"] = v

    # Status reflects "did we get anything at all" for this design.
    flat["tp_status"] = "ok" if any(v is not None for v in flat.values()) else None
    return flat


def _collect_complex(slug: str, model: str) -> dict[str, Any]:
    json_path = repo_root() / "data" / "metrics" / model / f"{slug}.json"
    ext = STRUCTURE_EXT.get(model, "cif")
    struct_path = repo_root() / "data" / "structures" / model / f"{slug}.{ext}"
    data = _read_json(json_path)
    fields = COMPLEX_FIELDS_BASE + EXTRA_FIELDS.get(model, ())
    row = _flatten(slug, model, data, fields)
    row[f"{PREFIX[model]}_struct_exists"] = struct_path.exists()
    return row


def _collect_sequence_scorer(slug: str, scorer: str) -> dict[str, Any]:
    """Pull a sequence-only scorer's per-design JSON into a flat row.

    Columns: ``<scorer>_status`` + every entry of ``SEQUENCE_SCORERS[scorer]``
    (already namespaced — we trust the upstream column names).
    """
    json_path = repo_root() / "data" / "metrics" / scorer / f"{slug}.json"
    data = _read_json(json_path)
    out: dict[str, Any] = {f"{scorer}_status": (data or {}).get("status")}
    for f in SEQUENCE_SCORERS[scorer]:
        # Disambiguate the shared ``length`` column between esm_pll and saprot.
        col = f if f.startswith(scorer) else f"{scorer}_{f}"
        out[col] = (data or {}).get(f)
    return out


def _collect_per_model_scorer(slug: str, scorer: str, model: str) -> dict[str, Any]:
    """Pull a (scorer, complex-model) JSON into a flat row.

    Columns: ``<scorer>_<model>_status`` + ``<scorer>_<model>_<field>`` for
    every entry of ``PER_MODEL_SCORERS[scorer]``.
    """
    json_path = (
        repo_root() / "data" / "metrics" / f"{scorer}_{model}" / f"{slug}.json"
    )
    data = _read_json(json_path)
    prefix = f"{scorer}_{model}"
    out: dict[str, Any] = {f"{prefix}_status": (data or {}).get("status")}
    for f in PER_MODEL_SCORERS[scorer]:
        # Strip the ``prodigy_`` / leading scorer namespace from the upstream
        # field name if present so the final column reads `prodigy_boltz2_kd`
        # instead of `prodigy_boltz2_prodigy_kd`.
        upstream = f
        clean = f
        if clean.startswith(f"{scorer}_"):
            clean = clean[len(scorer) + 1 :]
        out[f"{prefix}_{clean}"] = (data or {}).get(upstream)
    return out


def _add_consensus(df: pd.DataFrame) -> pd.DataFrame:
    """Cross-folder consensus columns: soft thresholds (ipSAE >= 0.4, iPTM >= 0.7)
    counted across every complex model we re-ran (boltz2, protenix, chai, af2m)."""
    n = len(COMPLEX_MODELS)
    ipsae_cols = [f"{PREFIX[m]}_ipsae_d0chn_max" for m in COMPLEX_MODELS]
    iptm_cols = [f"{PREFIX[m]}_iptm" for m in COMPLEX_MODELS]
    for c in ipsae_cols + iptm_cols:
        if c not in df.columns:
            df[c] = pd.NA
    df[f"ipsae_pass_{n}folders"] = (
        df[ipsae_cols].apply(pd.to_numeric, errors="coerce") >= 0.4
    ).sum(axis=1)
    df[f"iptm_pass_{n}folders"] = (
        df[iptm_cols].apply(pd.to_numeric, errors="coerce") >= 0.7
    ).sum(axis=1)
    return df


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out",
        default="data/grand_metrics.csv",
        help="Output CSV path, relative to repo root.",
    )
    args = parser.parse_args(argv)

    root = repo_root()
    designs = pd.read_csv(root / "data" / "designs.csv")

    designs_by_id = designs.set_index("design_id")
    rows: list[dict[str, Any]] = []
    for did in designs["design_id"].astype(int):
        slug = f"design_{did:03d}"
        des_row = designs_by_id.loc[did] if did in designs_by_id.index else None
        row: dict[str, Any] = {"slug": slug, "design_id": int(did)}
        row.update(_collect_typer(slug, des_row))
        for model in COMPLEX_MODELS:
            row.update(_collect_complex(slug, model))
        for scorer in SEQUENCE_SCORERS:
            row.update(_collect_sequence_scorer(slug, scorer))
        for scorer in PER_MODEL_SCORERS:
            for model in COMPLEX_MODELS:
                row.update(_collect_per_model_scorer(slug, scorer, model))
        rows.append(row)

    df = pd.DataFrame(rows)

    # Meta + label columns we always want exposed.
    meta = [c for c in ["design_id", "pb_id", "name", "team", "is_human",
                        "sequence_length", "is_hit", "submitted_to_lab",
                        "binding_label", "kd_arith_mean_nM_all"]
            if c in designs.columns]

    # Every pb_* column from designs.csv except the file-path stubs (those
    # are addressable via `data/structures/.../`, no point in carrying them
    # as wide-table columns). 100/141 coverage; the 41 non-screened designs
    # null across pb_*.
    pb_path_cols = {"pb_boltz2_cif", "pb_esmfold_cif", "pb_pae_json",
                    "pb_stylized_png"}
    pb_cols = [c for c in designs.columns
               if c.startswith("pb_") and c not in pb_path_cols]

    keep = meta + pb_cols
    df = df.merge(designs[keep], on="design_id", how="left", suffixes=("", "_dup"))
    # Drop any *_dup columns that arise when meta names collide (e.g. pb_id).
    df = df.drop(columns=[c for c in df.columns if c.endswith("_dup")])
    df = _add_consensus(df)

    leading = ["slug"] + meta
    cols = leading + [c for c in df.columns if c not in leading]
    df = df[cols]

    out_path = root / args.out
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)

    by_model = {m: int((df[f"{PREFIX[m]}_status"] == "ok").sum()) for m in COMPLEX_MODELS}
    n = len(COMPLEX_MODELS)
    half = (n + 1) // 2
    by_sequence_scorer = {
        s: int((df[f"{s}_status"] == "ok").sum())
        for s in SEQUENCE_SCORERS
        if f"{s}_status" in df.columns
    }
    by_per_model_scorer = {
        f"{s}_{m}": int((df[f"{s}_{m}_status"] == "ok").sum())
        for s in PER_MODEL_SCORERS
        for m in COMPLEX_MODELS
        if f"{s}_{m}_status" in df.columns
    }
    print(
        f"Wrote {out_path}  ({len(df)} rows x {len(df.columns)} cols)\n"
        f"  ok per complex model: {by_model}\n"
        f"  ok per sequence scorer: {by_sequence_scorer}\n"
        f"  ok per (scorer x model): {by_per_model_scorer}\n"
        f"  ipsae_pass_{n}folders >= {half}: {(df[f'ipsae_pass_{n}folders'] >= half).sum()}\n"
        f"  iptm_pass_{n}folders  >= {half}: {(df[f'iptm_pass_{n}folders'] >= half).sum()}"
    )


if __name__ == "__main__":
    main()
