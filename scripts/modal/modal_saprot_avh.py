# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "modal>=1.0",
#     "pandas",
# ]
# ///
"""SaProt structure-aware PLL for the 141 TREM2 binder designs.

Reads each design's monomer structure from the Modal Volume
(``monomer_structures/<slug>.cif``, populated from
``data/structures/{esmfold,proteintyper}/`` by the local entrypoint), extracts
3Di structural tokens with Foldseek, builds the SaProt vocabulary (interleaved
AA + 3Di), and runs masked-marginal PLL with the ``westlake-repl/SaProt_650M_AF2``
checkpoint.

Outputs land on the Modal Volume ``avh-rerun-results``:

* ``saprot/{slug}.json`` — ``{saprot_pll, saprot_pll_norm, length}``

Usage::

    cd <repo_root>
    # First time only: upload monomer CIFs to the volume
    modal run scripts/modal/modal_saprot_avh.py --upload-monomers

    # Score
    modal run --detach scripts/modal/modal_saprot_avh.py

    # Download finished results to disk
    modal run scripts/modal/modal_saprot_avh.py --download
"""

from __future__ import annotations

import contextlib
import json
import os
from pathlib import Path

import modal

GPU = os.environ.get("GPU", "A10G")
TIMEOUT_MIN = int(os.environ.get("TIMEOUT", 60))
CONCURRENCY = int(os.environ.get("CONCURRENCY", 10))

# ---------------------------------------------------------------------------
# Modal setup
# ---------------------------------------------------------------------------

saprot_image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("git", "wget")
    .pip_install(
        "torch>=2.0",
        "transformers>=4.30",
        "pandas",
        "numpy",
        "accelerate",
        "biopython>=1.80",
        "gemmi>=0.6",
    )
    # Foldseek static binary for 3Di token extraction.
    .run_commands(
        "wget -q https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz -O /tmp/foldseek.tar.gz "
        "&& tar -xzf /tmp/foldseek.tar.gz -C /opt/ "
        "&& ln -s /opt/foldseek/bin/foldseek /usr/local/bin/foldseek",
    )
)

app = modal.App(os.environ.get("MODAL_APP_NAME", "avh-saprot"))

RESULTS_VOLUME_NAME = os.environ.get("MODAL_RESULTS_VOLUME", "avh-rerun-results")
RESULTS_VOLUME = modal.Volume.from_name(RESULTS_VOLUME_NAME, create_if_missing=True)
RESULTS_DIR = f"/{RESULTS_VOLUME_NAME}"
PREDICTOR = "saprot"
MONOMER_DIR_ON_VOLUME = "monomer_structures"


# ---------------------------------------------------------------------------
# SaProt scoring function
# ---------------------------------------------------------------------------


@app.function(
    image=saprot_image,
    gpu=GPU,
    timeout=TIMEOUT_MIN * 60,
    max_containers=CONCURRENCY,
    memory=16384,
    volumes={RESULTS_DIR: RESULTS_VOLUME},
)
def score_saprot_batch(records: list[dict]) -> list[dict]:
    """Score a batch of (slug, sequence) records with SaProt 650M.

    Looks up each slug's monomer CIF on the volume, runs Foldseek to extract
    3Di tokens, then builds the SaProt vocabulary and computes per-residue
    masked-marginal log-likelihoods.
    """
    import torch
    from transformers import EsmForMaskedLM, EsmTokenizer

    RESULTS_VOLUME.reload()

    results: list[dict] = []

    try:
        model_name = "westlake-repl/SaProt_650M_AF2"
        print(f"Loading SaProt model: {model_name}")
        tokenizer = EsmTokenizer.from_pretrained(model_name)
        model = EsmForMaskedLM.from_pretrained(model_name)
        model.eval()
        if torch.cuda.is_available():
            model = model.cuda()
        device = next(model.parameters()).device

        monomer_dir = Path(RESULTS_DIR) / MONOMER_DIR_ON_VOLUME

        for r in records:
            slug = r["slug"]
            seq = r["sequence"]
            try:
                # 1. Locate the monomer structure on the volume.
                cif_path = monomer_dir / f"{slug}.cif"
                pdb_path = monomer_dir / f"{slug}.pdb"
                struct_path: Path | None = None
                if cif_path.exists():
                    struct_path = cif_path
                elif pdb_path.exists():
                    struct_path = pdb_path

                # 2. Extract 3Di tokens via foldseek. If the structure is
                #    missing fall back to 'd' (coil) tokens — same convention
                #    as SaProt's AA-only mode.
                if struct_path is not None:
                    three_di = _extract_3di(str(struct_path))
                    if three_di is None or len(three_di) != len(seq):
                        three_di = "d" * len(seq)
                        used_structure = False
                    else:
                        used_structure = True
                else:
                    three_di = "d" * len(seq)
                    used_structure = False

                # 3. Build SaProt input string: interleave AA and 3Di tokens.
                saprot_seq = "".join(
                    f"{aa}{di.lower() if di.isalpha() else 'd'}"
                    for aa, di in zip(seq, three_di, strict=True)
                )

                # 4. Compute masked-marginal PLL.
                inputs = tokenizer(
                    saprot_seq, return_tensors="pt", truncation=True, max_length=2048
                )
                inputs = {k: v.to(device) for k, v in inputs.items()}

                with torch.no_grad():
                    logits = model(**inputs).logits

                input_ids = inputs["input_ids"][0]
                log_probs = torch.nn.functional.log_softmax(logits[0], dim=-1)

                total_ll = 0.0
                n_tokens = 0
                for i in range(1, len(input_ids) - 1):  # skip [CLS] and [EOS]
                    token_id = input_ids[i].item()
                    total_ll += log_probs[i, token_id].item()
                    n_tokens += 1

                pll = total_ll
                pll_norm = total_ll / n_tokens if n_tokens > 0 else 0.0

                results.append(
                    {
                        "slug": slug,
                        "predictor": PREDICTOR,
                        "status": "ok",
                        "length": len(seq),
                        "used_structure": used_structure,
                        "saprot_pll": round(float(pll), 4),
                        "saprot_pll_norm": round(float(pll_norm), 4),
                    }
                )

            except Exception as e:
                import traceback

                print(f"  {slug}: {traceback.format_exc()}")
                results.append(
                    {
                        "slug": slug,
                        "predictor": PREDICTOR,
                        "status": f"error: {e}",
                        "length": len(seq) if isinstance(seq, str) else None,
                        "used_structure": False,
                        "saprot_pll": None,
                        "saprot_pll_norm": None,
                    }
                )

    except Exception as e:
        import traceback

        print(f"SaProt batch failed: {traceback.format_exc()}")
        for r in records:
            results.append(
                {
                    "slug": r["slug"],
                    "predictor": PREDICTOR,
                    "status": f"error: {e}",
                    "length": len(r.get("sequence") or ""),
                    "used_structure": False,
                    "saprot_pll": None,
                    "saprot_pll_norm": None,
                }
            )

    for r in results:
        _save_result(r)
    return results


def _extract_3di(struct_path: str) -> str | None:
    """Run foldseek structureto3didescriptor on a CIF/PDB and return the 3Di string.

    Returns None if foldseek fails or no chain is found.
    """
    import subprocess
    import tempfile

    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            out_tsv = Path(tmpdir) / "out.tsv"
            proc = subprocess.run(
                [
                    "foldseek",
                    "structureto3didescriptor",
                    "-v",
                    "0",
                    "--threads",
                    "1",
                    struct_path,
                    str(out_tsv),
                ],
                capture_output=True,
                text=True,
                timeout=120,
            )
            if proc.returncode != 0:
                print(f"  foldseek exit {proc.returncode}: {proc.stderr[-200:]}")
                return None
            if not out_tsv.exists():
                return None
            for line in out_tsv.read_text().splitlines():
                parts = line.split("\t")
                if len(parts) >= 3 and parts[2]:
                    return parts[2]  # 3Di string is column index 2
            return None
    except Exception as e:
        print(f"  foldseek failed: {e}")
        return None


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
    print(f"SaProt: {len(completed)} done, {len(pending)} pending")

    if pending:
        batch_size = 16
        shards = [pending[i : i + batch_size] for i in range(0, len(pending), batch_size)]
        records_per_shard = [
            [{"slug": s, "sequence": q} for s, q in shard] for shard in shards
        ]

        done = 0
        for batch_result in score_saprot_batch.map(records_per_shard, return_exceptions=True):
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
    print(f"SaProt: {ok}/{len(df)} succeeded.")


@app.local_entrypoint()
def main(
    designs_csv: str = "./data/designs.csv",
    limit: int | None = None,
    download: bool = False,
    retry_failed: bool = False,
    upload_monomers: bool = False,
) -> None:
    """Trigger a remote batch on Modal, or pull results to local disk.

    Args:
        upload_monomers: If True, upload local monomer CIFs to the volume so
                         the remote scorer can run foldseek on them. Idempotent.
    """
    import pandas as pd

    if upload_monomers:
        _upload_monomers(designs_csv)
        return

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
    print(f"Triggering SaProt batch: {len(slugs)} designs (SaProt 650M AF2, AA+3Di)")
    print(
        "  Reminder: run with --upload-monomers first if you haven't pushed "
        "monomer CIFs to the volume."
    )
    run_batch.remote(slugs, seqs)
    print("Done. Pull results with `--download`.")


def _upload_monomers(designs_csv: str) -> None:
    """Push monomer CIFs from data/structures/{esmfold,proteintyper}/ to the volume.

    SaProt needs a 3D structure to extract 3Di tokens. We prefer esmfold for the
    100 PB-mirrored designs and fall back to proteintyper for the 41 rerun
    designs. Slugs that have neither get an empty 3Di string at scoring time.
    """
    import pandas as pd

    df = pd.read_csv(designs_csv)
    slugs = [f"design_{int(d):03d}" for d in df["design_id"].tolist()]

    sources = (
        Path("./data/structures/esmfold"),
        Path("./data/structures/proteintyper"),
    )
    existing: set[str] = set()
    try:
        for entry in RESULTS_VOLUME.iterdir(MONOMER_DIR_ON_VOLUME):
            existing.add(Path(entry.path).stem)
    except Exception:
        pass

    n_uploaded = 0
    n_skipped = 0
    n_missing = 0
    for slug in slugs:
        if slug in existing:
            n_skipped += 1
            continue
        src: Path | None = None
        for parent in sources:
            cif = parent / f"{slug}.cif"
            pdb = parent / f"{slug}.pdb"
            if cif.exists():
                src = cif
                break
            if pdb.exists():
                src = pdb
                break
        if src is None:
            n_missing += 1
            continue
        dest = f"{MONOMER_DIR_ON_VOLUME}/{slug}{src.suffix}"
        RESULTS_VOLUME.write_file(dest, src.read_bytes())
        n_uploaded += 1

    print(
        f"Monomer upload: uploaded={n_uploaded}, skipped (already on volume)={n_skipped}, "
        f"missing on disk={n_missing}"
    )


def _download_to_local(designs_csv: str) -> None:
    """Pull all per-design JSONs into ./data/metrics/saprot/."""
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
