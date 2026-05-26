# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "modal>=1.0",
#     "pandas",
# ]
# ///
"""ESM2 t33 650M masked-marginal PLL on the 141 TREM2 binder designs.

Sequence-only scorer. Mirrors the shell of the complex predictors in this
folder (Modal Volume, ENV-overridable app name, ``design_NNN`` slug), but
runs no folding — we mask every position in turn and read back the log-prob
of the wild-type token to compute the masked-marginal pseudo-log-likelihood.

Outputs land on the Modal Volume ``avh-rerun-results``:

* ``esm_pll/{slug}.json`` — per-design ``{esm_pll_total, esm_pll_avg, length}``

Usage::

    cd <repo_root>
    modal run --detach scripts/modal/modal_esm_pll_avh.py

    # Download finished results to disk (idempotent, skip already-on-disk)
    modal run scripts/modal/modal_esm_pll_avh.py --download

    # Tune GPU / concurrency
    GPU=A100 CONCURRENCY=20 modal run --detach scripts/modal/modal_esm_pll_avh.py
"""

from __future__ import annotations

import contextlib
import json
import math
import os
from pathlib import Path

import modal

GPU = os.environ.get("GPU", "A10G")
TIMEOUT_MIN = int(os.environ.get("TIMEOUT", 30))
CONCURRENCY = int(os.environ.get("CONCURRENCY", 20))

# ---------------------------------------------------------------------------
# Modal setup
# ---------------------------------------------------------------------------

esm_image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install("fair-esm", "torch>=2.0", "pandas", "tqdm")
)

app = modal.App(os.environ.get("MODAL_APP_NAME", "avh-esm-pll"))

RESULTS_VOLUME_NAME = os.environ.get("MODAL_RESULTS_VOLUME", "avh-rerun-results")
RESULTS_VOLUME = modal.Volume.from_name(RESULTS_VOLUME_NAME, create_if_missing=True)
RESULTS_DIR = f"/{RESULTS_VOLUME_NAME}"
PREDICTOR = "esm_pll"


# ---------------------------------------------------------------------------
# Batched ESM2 PLL scorer
# ---------------------------------------------------------------------------


@app.function(
    image=esm_image,
    gpu=GPU,
    timeout=TIMEOUT_MIN * 60,
    max_containers=CONCURRENCY,
    volumes={RESULTS_DIR: RESULTS_VOLUME},
)
def compute_pll_batch(records: list[dict]) -> list[dict]:
    """Score a list of (slug, sequence) records with ESM2 t33 650M masked marginals.

    One masked forward pass per residue. Returns per-record dict with
    ``esm_pll_total``, ``esm_pll_avg``, ``length``.
    """
    import esm
    import torch

    RESULTS_VOLUME.reload()

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Loading ESM2 t33 650M on {device}")
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    model = model.eval().to(device)
    batch_converter = alphabet.get_batch_converter()

    def pll(seq: str) -> tuple[float, float]:
        data = [("p", seq)]
        *_, batch_tokens = batch_converter(data)
        log_probs: list[float] = []
        for i in range(len(seq)):
            masked = batch_tokens.clone()
            masked[0, i + 1] = alphabet.mask_idx
            with torch.no_grad():
                tp = torch.log_softmax(model(masked.to(device))["logits"], dim=-1)
            log_probs.append(tp[0, i + 1, alphabet.get_idx(seq[i])].item())
        s = math.fsum(log_probs)
        return s, s / len(seq)

    out: list[dict] = []
    for r in records:
        slug = r["slug"]
        seq = r["sequence"]
        try:
            total, avg = pll(seq)
            result = {
                "slug": slug,
                "predictor": PREDICTOR,
                "status": "ok",
                "length": len(seq),
                "esm_pll_total": total,
                "esm_pll_avg": avg,
            }
        except Exception as e:
            import traceback

            print(f"  {slug}: {traceback.format_exc()}")
            result = {
                "slug": slug,
                "predictor": PREDICTOR,
                "status": f"error: {e}",
                "length": len(seq) if isinstance(seq, str) else None,
                "esm_pll_total": None,
                "esm_pll_avg": None,
            }
        _save_result(result)
        out.append(result)
        print(f"  {slug}: pll_avg={result.get('esm_pll_avg')}")

    return out


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
    print(f"ESM PLL: {len(completed)} done, {len(pending)} pending")

    if pending:
        # One worker per CONCURRENCY shard — sequences are short and each
        # masked-marginal pass is cheap on an A10G.
        n_shards = max(1, min(CONCURRENCY, len(pending)))
        shards: list[list[dict]] = [[] for _ in range(n_shards)]
        for i, (slug, seq) in enumerate(pending):
            shards[i % n_shards].append({"slug": slug, "sequence": seq})

        done = 0
        for batch_result in compute_pll_batch.map(shards, return_exceptions=True):
            if isinstance(batch_result, Exception):
                print(f"  exception: {batch_result}")
                continue
            done += len(batch_result)
            print(f"  {done}/{len(pending)} done")

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
    print(f"ESM PLL: {ok}/{len(df)} succeeded.")


@app.local_entrypoint()
def main(
    designs_csv: str = "./data/designs.csv",
    limit: int | None = None,
    download: bool = False,
    retry_failed: bool = False,
) -> None:
    """Trigger a remote batch on Modal, or pull results to local disk.

    With ``--detach`` the batch survives client disconnects (Modal best
    practice for multi-hour jobs).
    """
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
    print(f"Triggering ESM PLL batch: {len(slugs)} designs (ESM2 t33 650M, masked marginals)")
    run_batch.remote(slugs, seqs)
    print("Done. Pull results with `--download`.")


def _download_to_local(designs_csv: str) -> None:
    """Pull all per-design JSONs into ./data/metrics/esm_pll/."""
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
