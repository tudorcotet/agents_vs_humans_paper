"""Re-run ProteinTyper on designs missing ESMFold predictions.

Submits each sequence to a ProteinTyper-compatible HTTP endpoint with
the ``full_monomer`` recipe (same recipe ProteinBase invokes when it
stamps `pb_esmfold_plddt` and friends). Polls the retrieve endpoint
until the job lands, then writes:

    data/structures/proteintyper/design_NNN.cif   ESMFold monomer
    data/metrics/proteintyper/design_NNN.json     full TyperJobOutput
    data/images/design_NNN.png                    stylised render

The targets are the 41 designs that ranked below the top-100 ipSAE
cutoff and so never went to ProteinBase / the wet lab. ProteinBase only
mirrored the 100 screened designs; this script fills in the rest so
analyses can quote ESMFold / typer-derived numbers across all 141.

Sequence-only input → ProteinTyper's default MSA strategy applies (HHblits
on the binder; no target chain in this monomer recipe — the BLI target
TREM2 is folded separately as part of the Boltz-2 complex pipeline,
which we do *not* rerun here).

Configuration (all via environment variables — no hosted defaults so
this script is safe to ship publicly):

    PROTEINTYPER_SUBMIT_URL   — full HTTPS URL of the submit endpoint
    PROTEINTYPER_RETRIEVE_URL — full HTTPS URL of the retrieve endpoint
    PROTEINTYPER_API_TOKEN    — Bearer token sent on every request

Adaptyv employees: the production values are documented in the private
``docs/INTERNAL.md`` (gitignored). External users running this against
their own ProteinTyper deployment should set the same three variables.

Run via::

    mise run rerun:typer
    # or:
    uv run python scripts/folding/run_proteintyper.py

The script is idempotent — designs whose typer JSON already exists on
disk are skipped. Pass ``--force`` to re-submit anyway. Pass
``--design-id 19 --design-id 32`` to target specific rows.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any

import pandas as pd

from scripts.utils import repo_root

POLL_INTERVAL_S = 30
POLL_TIMEOUT_S = 60 * 60   # 60-min hard cap per design (cold start + full_monomer recipe)
INITIAL_DELAY_S = 60       # don't poll for the first minute; level-0 alone is ~60–90 s


def _endpoint(name: str) -> str:
    url = os.environ.get(name)
    if not url:
        sys.stderr.write(
            f"[typer] {name} is unset. Set the three ProteinTyper env vars\n"
            f"(PROTEINTYPER_SUBMIT_URL, PROTEINTYPER_RETRIEVE_URL,\n"
            f"PROTEINTYPER_API_TOKEN). Adaptyv employees: see docs/INTERNAL.md\n"
            f"for the production values.\n"
        )
        sys.exit(1)
    return url


def _submit_url() -> str:
    return _endpoint("PROTEINTYPER_SUBMIT_URL")


def _retrieve_url() -> str:
    return _endpoint("PROTEINTYPER_RETRIEVE_URL")


def _auth_token() -> str:
    token = os.environ.get("PROTEINTYPER_API_TOKEN")
    if not token:
        sys.stderr.write(
            "[typer] PROTEINTYPER_API_TOKEN is unset. Set the three ProteinTyper\n"
            "env vars (PROTEINTYPER_SUBMIT_URL, PROTEINTYPER_RETRIEVE_URL,\n"
            "PROTEINTYPER_API_TOKEN). Adaptyv employees: see docs/INTERNAL.md\n"
            "for the production values.\n"
        )
        sys.exit(1)
    return token


def _post_json(url: str, payload: dict[str, Any], token: str, timeout: int = 60) -> dict[str, Any]:
    data = json.dumps(payload).encode()
    req = urllib.request.Request(
        url,
        data=data,
        headers={
            "Content-Type": "application/json",
            "Authorization": f"Bearer {token}",
        },
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read().decode())


def _submit(sequence: str, slug: str, token: str) -> str:
    """Submit a sequence to the configured ProteinTyper submit endpoint.

    Payload mirrors ProteinBase's own ``workers/proteintyper-receiver``: no
    ``target``, no ``webhook_*``. We pin ``recipe.template = full_monomer``
    so the full set of metrics ProteinBase ships (esmfold + proteinmpnn +
    novelty + classification + domainmatch) comes back instead of the
    minimal 2-metric default from the OpenAPI schema. Returns the typer_key.
    """
    payload = {
        "sequence": sequence,
        "client_id": slug,
        "start_typing_job": True,
        "recipe": {"template": "full_monomer"},
    }
    resp = _post_json(_submit_url(), payload, token)
    sequence_id = resp.get("sequence_id")
    if not sequence_id:
        raise RuntimeError(f"submit returned no sequence_id: {resp}")
    print(
        f"[typer] {slug}: submitted  sequence_id={sequence_id}  "
        f"job_started={resp.get('job_started')}  "
        f"requested_outputs={len(resp.get('requested_outputs') or [])} metrics"
    )
    return sequence_id


def _retrieve(typer_key: str, token: str) -> dict[str, Any] | None:
    """One retrieve call. Returns the TyperJobOutput dict when ready, else None.

    Current ProteinTyper schema::

        {"value": "<JSON string>", "schema": "..."}

    inner JSON shape::

        {"sequence": {"sequence": str, "metrics": [{metric_type, value}, ...]},
         "proteindomains": [...],
         "structures":      [{"file": {"url": "s3://..."}, "img": {...}, ...}]}

    Readiness signal: ``structures`` non-empty + first entry has ``file.url``.
    Falls back to the legacy dict-with-flat-fields shape so older deployments
    still work.

    The endpoint returns 500 during cold start and while the job is still
    running — both are transient. 404 / 425 / 5xx → None (poll again). Any
    other HTTP error raises.
    """
    try:
        resp = _post_json(_retrieve_url(), {"typer_key": typer_key}, token, timeout=30)
    except urllib.error.HTTPError as e:
        if e.code in (404, 425, 500, 502, 503, 504):
            return None
        body = e.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"retrieve HTTP {e.code}: {body}") from e

    raw = resp.get("value")
    if not raw:
        return None
    if isinstance(raw, str):
        try:
            value = json.loads(raw)
        except json.JSONDecodeError:
            return None
    elif isinstance(raw, dict):
        value = raw
    else:
        return None

    # New schema: structures[] gets populated when ESMFold finishes.
    if isinstance(value.get("structures"), list) and value["structures"]:
        first = value["structures"][0]
        if isinstance(first, dict) and isinstance(first.get("file"), dict) and first["file"].get("url"):
            return value

    # Legacy schema fallback.
    if value.get("esmfold_plddt") or value.get("design_class") or value.get("molecular_weight"):
        return value

    return None


def _poll_until_ready(typer_key: str, token: str, slug: str) -> dict[str, Any]:
    deadline = time.time() + POLL_TIMEOUT_S
    if INITIAL_DELAY_S:
        time.sleep(INITIAL_DELAY_S)
    waited = INITIAL_DELAY_S
    while time.time() < deadline:
        out = _retrieve(typer_key, token)
        if out is not None:
            populated = sum(1 for v in out.values() if v is not None)
            print(f"[typer] {slug}: ready after {waited}s ({populated} fields)")
            return out
        time.sleep(POLL_INTERVAL_S)
        waited += POLL_INTERVAL_S
        if waited % 120 == 0:
            print(f"[typer] {slug}: still running after {waited}s")
    raise RuntimeError(f"typer timed out after {POLL_TIMEOUT_S}s for {slug}")


def _download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": "agents_vs_humans_paper/0.1"})
    with urllib.request.urlopen(req, timeout=120) as r, dest.open("wb") as out:
        out.write(r.read())


_S3_PROTEINBASE_MIRROR = "https://proteinbase-pub.t3.storage.dev"


def _resolve_url(raw: Any) -> str | None:
    """Normalize a ProteinTyper file URL. The current API returns S3 URIs
    in the `proteinbase-pub` bucket (`s3://proteinbase-pub/<oid>.cif`);
    the public HTTP mirror behind it is `proteinbase-pub.t3.storage.dev`.
    Older deployments returned a full https:// URL directly."""
    if not isinstance(raw, str) or not raw:
        return None
    if raw.startswith("http"):
        return raw
    if raw.startswith("s3://proteinbase-pub/"):
        return f"{_S3_PROTEINBASE_MIRROR}/{raw.removeprefix('s3://proteinbase-pub/')}"
    return None


def _extract_url(value: Any) -> str | None:
    """Pull a URL from one of the legacy-schema fields (used as fallback)."""
    if isinstance(value, dict):
        return _resolve_url(value.get("url"))
    return _resolve_url(value) if isinstance(value, str) else None


def _structure_urls(typer_output: dict[str, Any]) -> tuple[str | None, str | None]:
    """Return (cif_url, png_url) for the predicted structure + stylised render.

    Tries the current schema first (`structures[0].{file,img}.url`), falls
    back to the legacy `esmfold_structure_prediction` / `esmfold_stylized_image`
    fields if present.
    """
    structures = typer_output.get("structures")
    if isinstance(structures, list) and structures:
        first = structures[0] if isinstance(structures[0], dict) else {}
        cif = _resolve_url((first.get("file") or {}).get("url"))
        png = _resolve_url((first.get("img") or {}).get("url"))
        if cif or png:
            return cif, png
    return (
        _extract_url(typer_output.get("esmfold_structure_prediction")),
        _extract_url(typer_output.get("esmfold_stylized_image")),
    )


def _slug(design_id: int) -> str:
    return f"design_{int(design_id):03d}"


def _save_artifacts(typer_output: dict[str, Any], slug: str, data_dir: Path) -> dict[str, Path]:
    """Save the TyperJobOutput JSON + download CIF / PNG artifacts.

    Layout::

        data/metrics/proteintyper/<slug>.json
        data/structures/proteintyper/<slug>.cif
        data/images/<slug>.png
    """
    out: dict[str, Path] = {}

    out_json = data_dir / "metrics" / "proteintyper" / f"{slug}.json"
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(typer_output, indent=2) + "\n")
    out["typer_output"] = out_json

    cif_url, png_url = _structure_urls(typer_output)
    if cif_url:
        dest = data_dir / "structures" / "proteintyper" / f"{slug}.cif"
        _download(cif_url, dest)
        out["esmfold_cif"] = dest

    if png_url:
        dest = data_dir / "images" / f"{slug}.png"
        _download(png_url, dest)
        out["stylized_png"] = dest

    return out


def _designs_missing_typer(df: pd.DataFrame) -> pd.DataFrame:
    """Designs without ProteinBase ESMFold scalars. 41 rows on the public release."""
    return df[df["pb_esmfold_plddt"].isna()].copy()


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--force", action="store_true",
        help="Resubmit even if data/metrics/proteintyper/<slug>.json exists",
    )
    parser.add_argument(
        "--design-id", action="append", type=int, default=None,
        help="Run on a specific design_id instead of the missing-typer set. Repeatable.",
    )
    args = parser.parse_args(argv)

    token = _auth_token()
    root = repo_root()
    df = pd.read_csv(root / "data" / "designs.csv")

    if args.design_id:
        targets = df[df["design_id"].isin(args.design_id)]
    else:
        targets = _designs_missing_typer(df)

    if targets.empty:
        print("[typer] no targets — every design already has ProteinTyper output.")
        return

    data_dir = root / "data"
    print(f"[typer] running on {len(targets)} designs")
    for _, row in targets.iterrows():
        slug = _slug(row["design_id"])
        sequence = row.get("sequence")
        if not isinstance(sequence, str) or not sequence:
            print(f"[typer] {slug}: no sequence, skipping")
            continue

        existing = data_dir / "metrics" / "proteintyper" / f"{slug}.json"
        if existing.exists() and not args.force:
            print(f"[typer] {slug}: already has output at {existing}, skipping")
            continue

        try:
            typer_key = _submit(sequence, slug, token)
            output = _poll_until_ready(typer_key, token, slug)
            paths = _save_artifacts(output, slug, data_dir)
            wrote = ", ".join(p.name for p in paths.values())
            print(f"[typer] {slug}: saved {wrote}")
        except Exception as e:
            print(f"[typer] {slug}: FAILED — {e}")


if __name__ == "__main__":
    main()
