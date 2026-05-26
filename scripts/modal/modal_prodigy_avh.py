# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "modal>=1.0",
#     "pandas",
# ]
# ///
"""PRODIGY binding-affinity scoring on pre-computed TREM2 complex structures.

Reads complex structures persisted on the Modal Volume by the complex
predictors (``structures/{boltz2,protenix,chai,af2m}/<slug>.{cif,pdb}``) and
runs PRODIGY's contact-based ΔG / Kd / pKd predictor on each.

Per-model outputs land on the volume:

* ``prodigy_{model}/{slug}.json`` — ``{prodigy_kd, prodigy_pkd, prodigy_dg, prodigy_temperature}``

Usage::

    cd <repo_root>
    # All four predictors
    modal run --detach scripts/modal/modal_prodigy_avh.py

    # Just one predictor
    modal run --detach scripts/modal/modal_prodigy_avh.py --predictors boltz2

    # Download finished results to disk
    modal run scripts/modal/modal_prodigy_avh.py --download
"""

from __future__ import annotations

import contextlib
import json
import os
import re
from pathlib import Path

import modal

TIMEOUT_MIN = int(os.environ.get("TIMEOUT", 15))
CONCURRENCY = int(os.environ.get("CONCURRENCY", 20))

# ---------------------------------------------------------------------------
# Modal setup
# ---------------------------------------------------------------------------

prodigy_image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install("biopython>=1.80", "prodigy-prot>=1.0", "pandas")
)

app = modal.App(os.environ.get("MODAL_APP_NAME", "avh-prodigy"))

RESULTS_VOLUME_NAME = os.environ.get("MODAL_RESULTS_VOLUME", "avh-rerun-results")
RESULTS_VOLUME = modal.Volume.from_name(RESULTS_VOLUME_NAME, create_if_missing=True)
RESULTS_DIR = f"/{RESULTS_VOLUME_NAME}"

ALL_PREDICTORS = ("boltz2", "protenix", "chai", "af2m")


# ---------------------------------------------------------------------------
# PRODIGY scoring function
# ---------------------------------------------------------------------------


@app.function(
    image=prodigy_image,
    volumes={RESULTS_DIR: RESULTS_VOLUME},
    timeout=TIMEOUT_MIN * 60,
    retries=1,
    max_containers=CONCURRENCY,
)
def score_prodigy(predictor: str, slug: str) -> dict:
    """Run PRODIGY on a persisted complex structure file.

    Selection: chain A (target, TREM2) vs chain B (binder).
    """
    import subprocess

    import numpy as np

    RESULTS_VOLUME.reload()

    struct_dir = Path(RESULTS_DIR) / "structures" / predictor
    pdb_path = struct_dir / f"{slug}.pdb"
    cif_path = struct_dir / f"{slug}.cif"

    if pdb_path.exists():
        struct_file = str(pdb_path)
    elif cif_path.exists():
        struct_file = str(cif_path)
    else:
        result = {
            "slug": slug,
            "predictor": predictor,
            "status": "error: no structure file",
            "prodigy_kd": None,
            "prodigy_pkd": None,
            "prodigy_dg": None,
            "prodigy_temperature": None,
        }
        _save_result(result, predictor)
        return result

    try:
        if struct_file.endswith(".cif"):
            from Bio.PDB import PDBIO, MMCIFParser

            parser = MMCIFParser(QUIET=True)
            structure = parser.get_structure("s", struct_file)
            pdb_tmp = struct_file.replace(".cif", "_prodigy.pdb")
            io = PDBIO()
            io.set_structure(structure)
            io.save(pdb_tmp)
            struct_file = pdb_tmp

        proc = subprocess.run(
            ["prodigy", struct_file, "--selection", "A", "B"],
            capture_output=True,
            text=True,
            timeout=120,
        )
        output = (proc.stdout or "") + (proc.stderr or "")

        # Parse Kd
        m = re.search(
            r"dissociation constant \(M\) at\s*([0-9.]+)[^\n]*C:\s*([0-9.eE+-]+)",
            output,
        )
        temperature: float | None = None
        if m:
            temperature = float(m.group(1))
            kd_str = m.group(2)
        else:
            m2 = re.search(r"Kd\s*=\s*([0-9.eE+-]+)", output)
            kd_str = m2.group(1) if m2 else None

        if kd_str is None:
            result = {
                "slug": slug,
                "predictor": predictor,
                "status": f"error: parse failed: {output[:200]}",
                "prodigy_kd": None,
                "prodigy_pkd": None,
                "prodigy_dg": None,
                "prodigy_temperature": temperature,
            }
            _save_result(result, predictor)
            return result

        kd = float(kd_str)
        pkd = float(-np.log10(kd)) if kd > 0 else None

        dg_match = re.search(r"binding (?:free )?energy.*?(-?[0-9.]+)\s*kcal", output)
        dg = float(dg_match.group(1)) if dg_match else None

        result = {
            "slug": slug,
            "predictor": predictor,
            "status": "ok",
            "prodigy_kd": kd,
            "prodigy_pkd": pkd,
            "prodigy_dg": dg,
            "prodigy_temperature": temperature,
        }
        _save_result(result, predictor)
        return result

    except Exception as e:
        import traceback

        print(f"PRODIGY failed for {predictor}/{slug}: {traceback.format_exc()}")
        result = {
            "slug": slug,
            "predictor": predictor,
            "status": f"error: {e}",
            "prodigy_kd": None,
            "prodigy_pkd": None,
            "prodigy_dg": None,
            "prodigy_temperature": None,
        }
        _save_result(result, predictor)
        return result


def _save_result(result: dict, predictor: str) -> None:
    out_dir = Path(RESULTS_DIR) / f"prodigy_{predictor}"
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / f"{result['slug']}.json").write_text(json.dumps(result))
    RESULTS_VOLUME.commit()


# ---------------------------------------------------------------------------
# Remote orchestrator
# ---------------------------------------------------------------------------

orchestrator_image = modal.Image.debian_slim(python_version="3.11").pip_install("pandas")


@app.function(
    image=orchestrator_image,
    timeout=12 * 3600,
    volumes={RESULTS_DIR: RESULTS_VOLUME},
)
def run_batch(predictors: list[str], slugs: list[str]) -> None:
    """Run PRODIGY on every persisted complex for each requested predictor."""
    import pandas as pd

    RESULTS_VOLUME.reload()

    for predictor in predictors:
        struct_dir = Path(RESULTS_DIR) / "structures" / predictor
        if not struct_dir.exists():
            print(f"No structures dir for {predictor}, skipping.")
            continue

        available_slugs: list[str] = []
        for f in struct_dir.iterdir():
            if f.suffix in (".pdb", ".cif"):
                stem = f.stem
                for suffix in ("_prodigy", "_foldx"):
                    if stem.endswith(suffix):
                        stem = stem[: -len(suffix)]
                if stem in slugs:
                    available_slugs.append(stem)
        available_slugs = sorted(set(available_slugs))
        print(f"PRODIGY/{predictor}: {len(available_slugs)} structures found")

        result_dir = Path(RESULTS_DIR) / f"prodigy_{predictor}"
        completed: set[str] = set()
        if result_dir.exists():
            for f in result_dir.glob("*.json"):
                with contextlib.suppress(Exception):
                    data = json.loads(f.read_text())
                    if data.get("status") == "ok":
                        completed.add(data["slug"])

        pending = [s for s in available_slugs if s not in completed]
        print(f"  {len(completed)} done, {len(pending)} pending")
        if not pending:
            continue

        for done, result in enumerate(
            score_prodigy.map(
                [predictor] * len(pending),
                pending,
                return_exceptions=True,
                wrap_returned_exceptions=False,
            ),
            start=1,
        ):
            if isinstance(result, Exception):
                print(f"  exception: {result}")
            elif done % 50 == 0:
                status = result.get("status", "unknown") if isinstance(result, dict) else "?"
                print(f"  {done}/{len(pending)} ({status})")

        RESULTS_VOLUME.reload()
        rows: list[dict] = []
        if result_dir.exists():
            for f in result_dir.glob("*.json"):
                with contextlib.suppress(Exception):
                    rows.append(json.loads(f.read_text()))

        if rows:
            df = pd.DataFrame(rows)
            csv_path = Path(RESULTS_DIR) / f"prodigy_{predictor}_summary.csv"
            df.to_csv(csv_path, index=False)
            RESULTS_VOLUME.commit()
            ok = int((df["status"] == "ok").sum())
            print(f"  prodigy_{predictor}: {ok}/{len(df)} ok.")


# ---------------------------------------------------------------------------
# Local entrypoint
# ---------------------------------------------------------------------------


@app.local_entrypoint()
def main(
    designs_csv: str = "./data/designs.csv",
    predictors: str = "all",
    limit: int | None = None,
    download: bool = False,
    retry_failed: bool = False,
) -> None:
    """Trigger PRODIGY scoring on persisted complexes, or pull results to disk.

    Args:
        predictors: comma-separated subset of {boltz2, protenix, chai, af2m},
                    or "all" for the full set.
    """
    import pandas as pd

    pred_list = (
        list(ALL_PREDICTORS) if predictors == "all" else [p.strip() for p in predictors.split(",")]
    )

    if download:
        _download_to_local(designs_csv, pred_list)
        return

    if retry_failed:
        for p in pred_list:
            _clear_failed(p)
        return

    df = pd.read_csv(designs_csv)
    df = df[df["sequence"].notna() & (df["sequence"].str.len() > 0)].copy()
    if limit:
        df = df.head(limit)

    slugs = [f"design_{int(d):03d}" for d in df["design_id"].tolist()]
    print(f"Triggering PRODIGY batch: {len(slugs)} slugs across predictors={pred_list}")
    run_batch.remote(pred_list, slugs)
    print("Done. Pull results with `--download`.")


def _download_to_local(designs_csv: str, predictors: list[str]) -> None:
    """Pull all per-design PRODIGY JSONs into ./data/metrics/prodigy_<model>/."""
    import pandas as pd

    df = pd.read_csv(designs_csv)
    expected = {f"design_{int(d):03d}" for d in df["design_id"].tolist()}

    total = 0
    for predictor in predictors:
        out_dir = Path("./data/metrics") / f"prodigy_{predictor}"
        out_dir.mkdir(parents=True, exist_ok=True)
        n = 0
        try:
            for entry in RESULTS_VOLUME.iterdir(f"prodigy_{predictor}"):
                if not entry.path.endswith(".json"):
                    continue
                payload = b"".join(RESULTS_VOLUME.read_file(entry.path))
                data = json.loads(payload)
                slug = data.get("slug") or Path(entry.path).stem
                if slug not in expected:
                    continue
                (out_dir / f"{slug}.json").write_bytes(payload)
                n += 1
        except Exception as e:
            print(f"  prodigy_{predictor}: no results yet ({e})")
        print(f"  prodigy_{predictor}: downloaded {n} JSONs.")
        total += n
    print(f"Downloaded {total} total JSONs.")


def _clear_failed(predictor: str) -> None:
    cleared = 0
    try:
        for entry in RESULTS_VOLUME.iterdir(f"prodigy_{predictor}"):
            if not entry.path.endswith(".json"):
                continue
            data = json.loads(b"".join(RESULTS_VOLUME.read_file(entry.path)))
            status = str(data.get("status", ""))
            if status.startswith("error") or status.startswith("failed"):
                RESULTS_VOLUME.remove_file(entry.path)
                cleared += 1
    except Exception:
        pass
    print(f"  prodigy_{predictor}: cleared {cleared} failed results")
