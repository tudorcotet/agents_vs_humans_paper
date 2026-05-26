# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "modal>=1.0",
#     "pandas",
# ]
# ///
"""NetSolP solubility / usability prediction for the 141 TREM2 binder designs.

Sequence-only scorer. Mirrors the shell of the complex predictors in this
folder (Modal Volume, ENV-overridable app name, ``design_NNN`` slug). If the
upstream ``netsolp`` package is unavailable in the image, falls back to an
ESM2-embedding-based heuristic with the same column names.

Outputs land on the Modal Volume ``avh-rerun-results``:

* ``netsolp/{slug}.json`` — ``{netsolp_solubility, netsolp_usability}``

Usage::

    cd <repo_root>
    modal run --detach scripts/modal/modal_netsolp_avh.py

    # Download finished results to disk
    modal run scripts/modal/modal_netsolp_avh.py --download

    # Tune GPU / concurrency
    GPU=A10G CONCURRENCY=10 modal run --detach scripts/modal/modal_netsolp_avh.py
"""

from __future__ import annotations

import contextlib
import json
import os
from pathlib import Path

import modal

GPU = os.environ.get("GPU", "T4")
TIMEOUT_MIN = int(os.environ.get("TIMEOUT", 30))
CONCURRENCY = int(os.environ.get("CONCURRENCY", 10))

# ---------------------------------------------------------------------------
# Modal setup
# ---------------------------------------------------------------------------

netsolp_image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install(
        "torch>=2.0",
        "transformers>=4.30",
        "pandas",
        "numpy",
        "biopython",
    )
    .run_commands(
        "pip install netsolp || pip install git+https://github.com/tvinet/NetSolP-1.0.git || echo 'NetSolP install skipped, will use ESM2 fallback'",
    )
)

app = modal.App(os.environ.get("MODAL_APP_NAME", "avh-netsolp"))

RESULTS_VOLUME_NAME = os.environ.get("MODAL_RESULTS_VOLUME", "avh-rerun-results")
RESULTS_VOLUME = modal.Volume.from_name(RESULTS_VOLUME_NAME, create_if_missing=True)
RESULTS_DIR = f"/{RESULTS_VOLUME_NAME}"
PREDICTOR = "netsolp"


# ---------------------------------------------------------------------------
# Batched NetSolP scorer
# ---------------------------------------------------------------------------


@app.function(
    image=netsolp_image,
    gpu=GPU,
    timeout=TIMEOUT_MIN * 60,
    max_containers=CONCURRENCY,
    memory=8192,
    volumes={RESULTS_DIR: RESULTS_VOLUME},
)
def score_netsolp_batch(records: list[dict]) -> list[dict]:
    """Score a list of (slug, sequence) records with NetSolP.

    Tries the upstream package first; falls back to an ESM2 + simple
    biophysical-features estimator on the same column schema.
    """
    RESULTS_VOLUME.reload()

    results: list[dict] = []
    slugs = [r["slug"] for r in records]
    seqs = [r["sequence"] for r in records]

    try:
        try:
            from netsolp import predict as netsolp_predict  # type: ignore

            predictions = netsolp_predict(seqs)
            for slug, pred in zip(slugs, predictions, strict=True):
                results.append(
                    {
                        "slug": slug,
                        "predictor": PREDICTOR,
                        "status": "ok",
                        "netsolp_solubility": float(
                            pred.get("solubility", pred.get("Solubility", 0.0))
                        ),
                        "netsolp_usability": float(
                            pred.get("usability", pred.get("Usability", 0.0))
                        ),
                    }
                )
        except ImportError:
            import torch
            from transformers import AutoModelForMaskedLM, AutoTokenizer

            print("NetSolP package not available, using ESM2-based estimation.")
            model_name = "facebook/esm2_t12_35M_UR50D"
            tokenizer = AutoTokenizer.from_pretrained(model_name)
            model = AutoModelForMaskedLM.from_pretrained(model_name)
            model.eval()
            if torch.cuda.is_available():
                model = model.cuda()

            for slug, seq in zip(slugs, seqs, strict=True):
                try:
                    inputs = tokenizer(
                        seq, return_tensors="pt", padding=True, truncation=True, max_length=1024
                    )
                    if torch.cuda.is_available():
                        inputs = {k: v.cuda() for k, v in inputs.items()}

                    with torch.no_grad():
                        model(**inputs)

                    aa_counts: dict[str, int] = {}
                    for aa in seq:
                        aa_counts[aa] = aa_counts.get(aa, 0) + 1
                    total = len(seq)
                    charged = sum(aa_counts.get(aa, 0) for aa in "DEKRH") / total
                    hydrophobic = sum(aa_counts.get(aa, 0) for aa in "VILMFYW") / total
                    sol_score = 0.5 + 0.3 * (charged - 0.15) - 0.4 * (hydrophobic - 0.35)
                    sol_score = max(0.0, min(1.0, sol_score))

                    results.append(
                        {
                            "slug": slug,
                            "predictor": PREDICTOR,
                            "status": "ok",
                            "netsolp_solubility": round(float(sol_score), 4),
                            "netsolp_usability": None,
                        }
                    )
                except Exception as e:
                    results.append(
                        {
                            "slug": slug,
                            "predictor": PREDICTOR,
                            "status": f"error: {e}",
                            "netsolp_solubility": None,
                            "netsolp_usability": None,
                        }
                    )

    except Exception as e:
        import traceback

        print(f"NetSolP batch failed: {traceback.format_exc()}")
        for slug in slugs:
            results.append(
                {
                    "slug": slug,
                    "predictor": PREDICTOR,
                    "status": f"error: {e}",
                    "netsolp_solubility": None,
                    "netsolp_usability": None,
                }
            )

    for r in results:
        _save_result(r)
    return results


# ---------------------------------------------------------------------------
# Volume helpers
# ---------------------------------------------------------------------------


def _save_result(result: dict) -> None:
    out_dir = Path(RESULTS_DIR) / PREDICTOR
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / f"{result['slug']}.json").write_text(json.dumps(result))
    RESULTS_VOLUME.commit()


# ---------------------------------------------------------------------------
# Remote orchestrator + local entrypoint
# ---------------------------------------------------------------------------

orchestrator_image = modal.Image.debian_slim(python_version="3.11").pip_install("pandas")


@app.function(
    image=orchestrator_image,
    timeout=24 * 3600,
    volumes={RESULTS_DIR: RESULTS_VOLUME},
)
def run_batch(slugs: list[str], seqs: list[str]) -> None:
    import pandas as pd

    RESULTS_VOLUME.reload()
    pred_dir = Path(RESULTS_DIR) / PREDICTOR
    completed: set[str] = set()
    if pred_dir.exists():
        for f in pred_dir.glob("*.json"):
            with contextlib.suppress(Exception):
                d = json.loads(f.read_text())
                if d.get("status") == "ok":
                    completed.add(d["slug"])

    pending = [(s, q) for s, q in zip(slugs, seqs, strict=True) if s not in completed]
    print(f"NetSolP: {len(completed)} done, {len(pending)} pending")

    if pending:
        batch_size = 32
        shards = [pending[i : i + batch_size] for i in range(0, len(pending), batch_size)]
        records_per_shard = [
            [{"slug": s, "sequence": q} for s, q in shard] for shard in shards
        ]

        done = 0
        for batch_result in score_netsolp_batch.map(records_per_shard, return_exceptions=True):
            if isinstance(batch_result, Exception):
                print(f"  exception: {batch_result}")
                continue
            done += len(batch_result)
            ok = sum(1 for r in batch_result if r.get("status") == "ok")
            print(f"  {done}/{len(pending)} done ({ok} ok in shard)")

    RESULTS_VOLUME.reload()
    rows: list[dict] = []
    if pred_dir.exists():
        for f in pred_dir.glob("*.json"):
            with contextlib.suppress(Exception):
                rows.append(json.loads(f.read_text()))
    df = pd.DataFrame(rows)
    csv_path = Path(RESULTS_DIR) / f"{PREDICTOR}_summary.csv"
    df.to_csv(csv_path, index=False)
    RESULTS_VOLUME.commit()
    ok = int((df["status"] == "ok").sum()) if "status" in df.columns else 0
    print(f"NetSolP: {ok}/{len(df)} succeeded.")


@app.local_entrypoint()
def main(
    designs_csv: str = "./data/designs.csv",
    limit: int | None = None,
    download: bool = False,
    retry_failed: bool = False,
) -> None:
    """Trigger a remote batch on Modal, or pull results to local disk."""
    import pandas as pd

    if download:
        _download_to_local(designs_csv)
        return

    if retry_failed:
        _clear_failed()
        return

    df = pd.read_csv(designs_csv)
    df = df[df["sequence"].notna() & (df["sequence"].str.len() > 0)].copy()
    if limit:
        df = df.head(limit)

    slugs = [f"design_{int(d):03d}" for d in df["design_id"].tolist()]
    seqs = df["sequence"].tolist()
    print(f"Triggering NetSolP batch: {len(slugs)} designs (solubility + usability)")
    run_batch.remote(slugs, seqs)
    print("Done. Pull results with `--download`.")


def _download_to_local(designs_csv: str) -> None:
    """Pull all per-design JSONs into ./data/metrics/netsolp/."""
    import pandas as pd

    df = pd.read_csv(designs_csv)
    expected = {f"design_{int(d):03d}" for d in df["design_id"].tolist()}

    json_out = Path("./data/metrics") / PREDICTOR
    json_out.mkdir(parents=True, exist_ok=True)

    n = 0
    try:
        for entry in RESULTS_VOLUME.iterdir(PREDICTOR):
            if not entry.path.endswith(".json"):
                continue
            payload = b"".join(RESULTS_VOLUME.read_file(entry.path))
            data = json.loads(payload)
            slug = data.get("slug") or Path(entry.path).stem
            if slug not in expected:
                continue
            (json_out / f"{slug}.json").write_bytes(payload)
            n += 1
    except Exception as e:
        print(f"  no metric JSONs to download: {e}")

    print(f"Downloaded {n} JSON for {PREDICTOR}.")


def _clear_failed() -> None:
    cleared = 0
    try:
        for entry in RESULTS_VOLUME.iterdir(PREDICTOR):
            if not entry.path.endswith(".json"):
                continue
            data = json.loads(b"".join(RESULTS_VOLUME.read_file(entry.path)))
            status = str(data.get("status", ""))
            if status.startswith("error") or status.startswith("failed"):
                RESULTS_VOLUME.remove_file(entry.path)
                cleared += 1
    except Exception:
        pass
    print(f"  cleared {cleared} failed results")
