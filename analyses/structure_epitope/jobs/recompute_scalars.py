"""Recompute the corrupt b2_*/af2m_* scalar columns from the REGENERATED complexes.

The shipped ``data/structures/{boltz2,af2m}/`` had the wrong binder chain for most
designs (issue #5); the ``b2_*``/``af2m_*`` columns in ``data/designs.{csv,parquet}``
were computed from those mislabeled files and are corrupt with them. This script
recomputes all 38 columns (plus the two ``*_pass_4folders`` consensus columns) from
the regenerated, integrity-checked Boltz-2 / AF2-Multimer structures + their PAE.

It uses the repo's OWN upstream scorer — ``scripts/modal/modal_boltz2_avh.compute_ipsae``
(PAE_CUTOFF=15, DIST_CUTOFF=8) — imported directly (modal stubbed) so there is zero
drift from the code that produced the clean ``px_*``/``chai_*`` columns. The corrected
values are therefore drop-in consistent with the other two predictors.

Context + validation: analyses/structure_epitope/STRUCTURE_REGENERATION.md

Usage (from repo root):
    uv run --with numpy --with gemmi python analyses/structure_epitope/jobs/recompute_scalars.py            # compute + validate, no write
    uv run --with numpy --with gemmi python analyses/structure_epitope/jobs/recompute_scalars.py --apply    # + patch designs.{csv,parquet}
"""

from __future__ import annotations

import argparse
import glob
import importlib.util
import io
import json
import os
import sys
from unittest import mock

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
REGEN = os.path.join(ROOT, "analyses", "structure_epitope", "data", "structures_regen")
BOLTZ_PRED = os.path.join(REGEN, "boltz_results_prod_yamls", "predictions")
AF2M_DIR = os.path.join(REGEN, "af2m")
TARGET_LEN = 175

# COMPLEX_FIELDS_BASE order from scripts/data/build_grand_metrics.py (native scalars
# first, then the compute_ipsae interface metrics). status + struct_exists bracket them.
FIELDS = (
    "iptm", "ptm", "mean_plddt",
    "ipsae_d0res_min", "ipsae_d0res_max",
    "ipsae_d0chn_min", "ipsae_d0chn_max",
    "ipsae_d0dom_min", "ipsae_d0dom_max",
    "iptm_d0chn_min", "iptm_d0chn_max",
    "iptm_af_min", "iptm_af_max",
    "pdockq", "pdockq2", "lis", "n_interface",
)


def _load_repo_scorer():
    """Import the exact upstream compute_ipsae with modal stubbed out (no drift)."""
    sys.modules["modal"] = mock.MagicMock()
    path = os.path.join(ROOT, "scripts", "modal", "modal_boltz2_avh.py")
    spec = importlib.util.spec_from_file_location("_mb_scorer", path)
    mb = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mb)
    return mb.compute_ipsae


def _boltz_native(conf: dict) -> tuple[float, float, float]:
    iptm = float(conf.get("iptm", 0.0))
    ptm = float(conf.get("ptm", 0.0))
    if "plddt" in conf:
        p = conf["plddt"]
        plddt = float(np.mean(p)) if isinstance(p, list) else float(p)
    elif "complex_plddt" in conf:
        plddt = float(conf["complex_plddt"])
    else:
        plddt = 0.0
    return iptm, ptm, plddt


def _af2m_native(scores: dict) -> tuple[float, float, float]:
    iptm = float(scores.get("iptm", 0.0))
    ptm = float(scores.get("ptm", 0.0))
    plddt = float(np.mean(np.array(scores["plddt"]))) if "plddt" in scores else 0.0
    return iptm, ptm, plddt


def recompute(compute_ipsae, seq_by_id: dict[int, str]) -> dict[str, dict[int, dict]]:
    """Return {'b2': {design_id: row}, 'af2m': {design_id: row}} with FIELDS + status."""
    out = {"b2": {}, "af2m": {}}
    for n in range(1, 142):
        binder_len = len(seq_by_id[n])
        # --- Boltz-2 ---
        d = f"{BOLTZ_PRED}/design_{n:03d}"
        cif = f"{d}/design_{n:03d}_model_0.cif"
        npz = f"{d}/pae_design_{n:03d}_model_0.npz"
        conf = json.load(open(f"{d}/confidence_design_{n:03d}_model_0.json"))
        pae = np.load(npz)["pae"]
        if pae.ndim == 3:
            pae = pae[0]
        sc = compute_ipsae(pae, cif, TARGET_LEN, binder_len, "A", "B")
        sc["iptm"], sc["ptm"], sc["mean_plddt"] = _boltz_native(conf)
        sc["status"], sc["struct_exists"] = "ok", True
        out["b2"][n] = sc
        # --- AF2-Multimer ---
        pdb = glob.glob(f"{AF2M_DIR}/design_{n:03d}_unrelaxed_rank_001_*.pdb")[0]
        scores = json.load(open(glob.glob(f"{AF2M_DIR}/design_{n:03d}_scores_rank_001_*.json")[0]))
        pae_a = np.array(scores["pae"])
        sa = compute_ipsae(pae_a, pdb, TARGET_LEN, binder_len, "A", "B")
        sa["iptm"], sa["ptm"], sa["mean_plddt"] = _af2m_native(scores)
        sa["status"], sa["struct_exists"] = "ok", True
        out["af2m"][n] = sa
    return out


def build_frames(new: dict, ids: list[int]) -> dict[str, pd.DataFrame]:
    """Numeric DataFrames of the new b2_*/af2m_* columns, indexed by design_id."""
    frames = {}
    for pred, prefix in [("b2", "b2"), ("af2m", "af2m")]:
        rows = {}
        for n in ids:
            r = new[pred][n]
            row = {f"{prefix}_status": r["status"]}
            for f in FIELDS:
                row[f"{prefix}_{f}"] = r[f]
            row[f"{prefix}_n_interface"] = int(r["n_interface"])
            row[f"{prefix}_struct_exists"] = bool(r["struct_exists"])
            rows[n] = row
        frames[prefix] = pd.DataFrame.from_dict(rows, orient="index")
    return frames


def add_consensus(df: pd.DataFrame) -> pd.DataFrame:
    """Reproduce build_grand_metrics._add_consensus with the (now corrected) columns."""
    models = ["b2", "px", "chai", "af2m"]
    ipsae_cols = [f"{m}_ipsae_d0chn_max" for m in models]
    iptm_cols = [f"{m}_iptm" for m in models]
    df["ipsae_pass_4folders"] = (
        df[ipsae_cols].apply(pd.to_numeric, errors="coerce") >= 0.4
    ).sum(axis=1).astype(int)
    df["iptm_pass_4folders"] = (
        df[iptm_cols].apply(pd.to_numeric, errors="coerce") >= 0.7
    ).sum(axis=1).astype(int)
    return df


def _fmt(v) -> str:
    if isinstance(v, bool):
        return "True" if v else "False"
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if isinstance(v, float) or isinstance(v, np.floating):
        return repr(float(v))
    return str(v)


def validate(df_old: pd.DataFrame, df_new: pd.DataFrame, frames: dict) -> None:
    print("\n" + "=" * 78)
    print("VALIDATION")
    print("=" * 78)
    ids = df_new["design_id"].tolist()

    # (1) iptm cross-check vs the independently-extracted regen scalars
    for pred, csv, col in [
        ("b2", f"{REGEN}/boltz2_regen_scalars.csv", "b2new_iptm"),
        ("af2m", f"{REGEN}/af2m_regen_scalars.csv", "af2new_iptm"),
    ]:
        chk = pd.read_csv(csv).set_index("design_id")[col]
        mine = df_new.set_index("design_id")[f"{pred}_iptm"]
        d = (chk.reindex(ids).values - mine.reindex(ids).values)
        print(f"[iptm cross-check] {pred}_iptm vs {col}: max|Δ|={np.nanmax(np.abs(d)):.2e}  (141 designs)")

    def sp(a, b):
        m = df_new[a].notna() & df_new[b].notna()
        return stats.spearmanr(df_new.loc[m, a], df_new.loc[m, b]).statistic if m.sum() > 3 else float("nan")

    # (2) new corrected columns agree with the clean references
    print("\n[agreement with clean refs — Spearman ρ]")
    for pred in ["b2", "af2m"]:
        print(f"  {pred}_iptm vs   submitted_ipsae={sp(f'{pred}_iptm','submitted_ipsae'):.2f}"
              f"  boltz2_iptm={sp(f'{pred}_iptm','boltz2_iptm'):.2f}"
              f"  px_iptm={sp(f'{pred}_iptm','px_iptm'):.2f}"
              f"  chai_iptm={sp(f'{pred}_iptm','chai_iptm'):.2f}")
        print(f"  {pred}_ipsae_d0chn_max vs px={sp(f'{pred}_ipsae_d0chn_max','px_ipsae_d0chn_max'):.2f}"
              f"  chai={sp(f'{pred}_ipsae_d0chn_max','chai_ipsae_d0chn_max'):.2f}")

    # (3) new vs OLD corrupt — should be ~0 (independent), confirming the columns changed
    print("\n[new vs OLD corrupt — Spearman ρ, expect ≈0]")
    for pred in ["b2", "af2m"]:
        for suf in ["iptm", "ipsae_d0chn_max"]:
            o = df_old.set_index("design_id")[f"{pred}_{suf}"]
            n = df_new.set_index("design_id")[f"{pred}_{suf}"]
            m = o.notna() & n.notna()
            rho = stats.spearmanr(o[m], n[m]).statistic
            print(f"  {pred}_{suf}: ρ(old,new)={rho:.2f}")

    # (4) binder discrimination AUROC among expressed
    try:
        from sklearn.metrics import roc_auc_score
        e = df_new[df_new["is_hit"].notna()].copy()
        if "expressed" in df_new.columns:
            e = df_new[df_new["expressed"] == True].copy()  # noqa: E712
        for pred in ["b2", "af2m"]:
            for suf in ["iptm", "ipsae_d0chn_max"]:
                s = e.dropna(subset=[f"{pred}_{suf}"])
                y = (s["is_hit"] == True).astype(int)  # noqa: E712
                if y.nunique() == 2:
                    a = max(roc_auc_score(y, s[f"{pred}_{suf}"]), roc_auc_score(y, -s[f"{pred}_{suf}"]))
                    print(f"[AUROC] {pred}_{suf}: {a:.2f}")
    except Exception as ex:
        print("AUROC skipped:", ex)

    # (5) consensus distributions old vs new
    print("\n[consensus counts, value: #designs]")
    for col in ["ipsae_pass_4folders", "iptm_pass_4folders"]:
        print(f"  OLD {col}: {df_old[col].value_counts().sort_index().to_dict()}")
        print(f"  NEW {col}: {df_new[col].value_counts().sort_index().to_dict()}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--apply", action="store_true", help="write designs.{csv,parquet}")
    args = ap.parse_args()

    compute_ipsae = _load_repo_scorer()
    pq_path = os.path.join(ROOT, "data", "designs.parquet")
    csv_path = os.path.join(ROOT, "data", "designs.csv")
    df_old = pd.read_parquet(pq_path)
    seq_by_id = dict(zip(df_old["design_id"].astype(int), df_old["sequence"]))
    ids = df_old["design_id"].astype(int).tolist()

    print(f"Recomputing b2_*/af2m_* for {len(ids)} designs "
          f"(scorer=modal_boltz2_avh.compute_ipsae, pae_cutoff=15, dist_cutoff=8)...")
    new = recompute(compute_ipsae, seq_by_id)
    frames = build_frames(new, ids)

    # assemble the corrected numeric frame
    df_new = df_old.copy()
    for prefix in ["b2", "af2m"]:
        f = frames[prefix]
        for col in f.columns:
            df_new[col] = df_new["design_id"].astype(int).map(f[col])
    df_new = add_consensus(df_new)

    validate(df_old, df_new, frames)

    changed = list(frames["b2"].columns) + list(frames["af2m"].columns) + ["ipsae_pass_4folders", "iptm_pass_4folders"]
    print(f"\nColumns to patch: {len(changed)} "
          f"({len(frames['b2'].columns)} b2_ + {len(frames['af2m'].columns)} af2m_ + 2 consensus)")

    # full corrected-scalars record
    rec_path = os.path.join(ROOT, "analyses", "structure_epitope", "results", "recomputed_scalars.csv")
    df_new[["design_id"] + changed].to_csv(rec_path, index=False)
    print(f"Wrote {os.path.relpath(rec_path, ROOT)}")

    if not args.apply:
        print("\n[dry run] designs.{csv,parquet} NOT modified. Re-run with --apply to write.")
        return

    # 1) parquet — numeric, exact
    df_new.to_parquet(pq_path, index=False)
    print(f"Patched {os.path.relpath(pq_path, ROOT)}")

    # 2) csv — read as strings so untouched columns round-trip verbatim; overwrite only target cells
    s = pd.read_csv(csv_path, dtype=str, keep_default_na=False)
    assert (s["design_id"].astype(int).tolist() == ids), "designs.csv row order != parquet"
    idmap = {int(v): i for i, v in enumerate(s["design_id"])}
    for col in changed:
        vals = df_new.set_index("design_id")[col]
        for did, i in idmap.items():
            s.at[i, col] = _fmt(vals.loc[did])
    buf = io.StringIO()
    s.to_csv(buf, index=False)
    open(csv_path, "w").write(buf.getvalue())
    print(f"Patched {os.path.relpath(csv_path, ROOT)}")
    print("\nDone. b2_*/af2m_* + consensus now reflect the regenerated structures.")


if __name__ == "__main__":
    main()
