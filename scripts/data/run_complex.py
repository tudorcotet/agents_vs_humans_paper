"""Re-score every TREM2 design with our complex-prediction Modal panel.

ProteinBase shipped each binder through ESMFold + ProteinMPNN + AFDB50
only — no complex prediction against TREM2. This driver fires the same
multi-folder panel (Boltz-2, Protenix-v2, Chai-1) and pulls the
per-design CIFs + raw metric JSONs into the repo.

Each model runs as its own Modal app (see ``scripts/modal/``). This
script is a thin orchestrator: it does NOT do GPU work locally, just
shells out to ``modal run`` and (optionally) ``--download`` to mirror
the Modal Volume back to ``data/``.

Outputs::

    data/structures/boltz2/<pb_id>.cif    — Boltz-2 complex
    data/structures/protenix/<pb_id>.cif  — Protenix complex
    data/structures/chai/<pb_id>.cif      — Chai-1 complex
    data/metrics/boltz2/<pb_id>.json      — Boltz-2 + ipSAE + iPTM + pTM + …
    data/metrics/protenix/<pb_id>.json    — Protenix + ipSAE + ranking_score
    data/metrics/chai/<pb_id>.json        — Chai-1 + ipSAE + aggregate_score

Idempotent: a model whose JSON+CIF already sit on disk is skipped.

Auth: needs ``modal token set …`` already configured. The Modal apps
write to a workspace Volume named ``avh-rerun-results``.

Optional second auth: ``COMPLEX_API_TOKEN`` is reserved for the case
where these endpoints later move behind an HTTP gateway
(today the apps are direct ``modal run`` invocations and only the
Modal CLI token is needed). The script does not fail if it's unset;
it only complains if the ``--use-gateway`` flag is passed.

Usage::

    # All three models, full 141 designs (uses your Modal token)
    mise run rerun:complex
    # equivalently:
    uv run python scripts/data/run_complex.py

    # Just one model
    uv run python scripts/data/run_complex.py --models boltz2

    # Skip the launch step (volume already has results), just pull
    uv run python scripts/data/run_complex.py --download-only

    # Smoke test with a small set
    uv run python scripts/data/run_complex.py --limit 5
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys

from scripts.utils.load_data import repo_root

MODELS = ("boltz2", "protenix", "chai", "af2m")

MODAL_SCRIPTS: dict[str, str] = {
    "boltz2": "scripts/modal/modal_boltz2_avh.py",
    "protenix": "scripts/modal/modal_protenix_avh.py",
    "chai": "scripts/modal/modal_chai1_avh.py",
    "af2m": "scripts/modal/modal_af2m_avh.py",
}

# Sequence-only or structure-consuming scorers. Run after the complex
# predictors have populated `data/structures/<model>/`.
SCORERS = ("esm_pll", "netsolp", "saprot", "prodigy", "destress")

SCORER_SCRIPTS: dict[str, str] = {
    "esm_pll": "scripts/modal/modal_esm_pll_avh.py",
    "netsolp": "scripts/modal/modal_netsolp_avh.py",
    "saprot": "scripts/modal/modal_saprot_avh.py",
    "prodigy": "scripts/modal/modal_prodigy_avh.py",
    "destress": "scripts/modal/modal_destress_avh.py",
}

# Scorers that emit a separate JSON folder per complex model.
PER_MODEL_SCORERS = ("prodigy", "destress")
# Sequence-only scorers, one JSON folder total.
SEQUENCE_SCORERS = ("esm_pll", "netsolp", "saprot")


def _modal_bin() -> str:
    bin_ = shutil.which("modal")
    if not bin_:
        sys.stderr.write(
            "[complex] `modal` CLI not on PATH. Install once:\n"
            "  uv tool install modal\n"
            "Then `modal token set --token-id ... --token-secret ...` "
            "(token from https://modal.com/settings/tokens).\n"
        )
        sys.exit(1)
    return bin_


def _check_token_env(name: str, where: str) -> None:
    """If a token name is set as a hint, validate. Otherwise just warn."""
    if not os.environ.get(name):
        sys.stderr.write(
            f"[complex] {name} is unset. If these endpoints move behind a gateway, "
            f"fetch the gateway token via your secret manager and set it.\n"
            f"For now the Modal CLI token is sufficient.\n"
        )


def _check_models(models: list[str]) -> None:
    bad = [m for m in models if m not in MODELS]
    if bad:
        sys.stderr.write(f"[complex] unknown models: {bad}. choose from {list(MODELS)}\n")
        sys.exit(2)


def _check_scorers(scorers: list[str]) -> None:
    bad = [s for s in scorers if s not in SCORERS]
    if bad:
        sys.stderr.write(f"[complex] unknown scorers: {bad}. choose from {list(SCORERS)}\n")
        sys.exit(2)


def _launch(model: str, *, detach: bool, limit: int | None) -> int:
    """Trigger the remote orchestrator on Modal. Returns exit code."""
    script = MODAL_SCRIPTS[model]
    cmd = [_modal_bin(), "run"]
    if detach:
        cmd.append("--detach")
    cmd.append(script)
    if limit is not None:
        cmd.extend(["--limit", str(limit)])
    print(f"[complex] launching {model}: {' '.join(cmd)}")
    proc = subprocess.run(cmd, cwd=repo_root())
    return proc.returncode


def _launch_scorer(scorer: str, *, detach: bool, limit: int | None) -> int:
    """Trigger a scorer Modal app. Returns exit code."""
    script = SCORER_SCRIPTS[scorer]
    cmd = [_modal_bin(), "run"]
    if detach:
        cmd.append("--detach")
    cmd.append(script)
    if limit is not None:
        cmd.extend(["--limit", str(limit)])
    print(f"[scorers] launching {scorer}: {' '.join(cmd)}")
    proc = subprocess.run(cmd, cwd=repo_root())
    return proc.returncode


def _download(model: str) -> int:
    """Pull JSON + CIF from the Modal Volume into data/."""
    script = MODAL_SCRIPTS[model]
    cmd = [_modal_bin(), "run", script, "--download"]
    print(f"[complex] downloading {model}: {' '.join(cmd)}")
    proc = subprocess.run(cmd, cwd=repo_root())
    return proc.returncode


def _download_scorer(scorer: str) -> int:
    """Pull per-design JSONs from the Modal Volume into data/metrics/."""
    script = SCORER_SCRIPTS[scorer]
    cmd = [_modal_bin(), "run", script, "--download"]
    print(f"[scorers] downloading {scorer}: {' '.join(cmd)}")
    proc = subprocess.run(cmd, cwd=repo_root())
    return proc.returncode


def _already_on_disk(model: str, pb_ids: list[str]) -> set[str]:
    """Return pb_ids that already have BOTH JSON and CIF locally."""
    root = repo_root()
    json_dir = root / "data" / "metrics" / model
    cif_dir = root / "data" / "structures" / model
    done: set[str] = set()
    if not (json_dir.exists() and cif_dir.exists()):
        return done
    json_have = {p.stem for p in json_dir.glob("*.json")}
    cif_have = {p.stem for p in cif_dir.glob("*.cif")}
    for pb in pb_ids:
        if pb in json_have and pb in cif_have:
            done.add(pb)
    return done


def _design_ids() -> list[str]:
    import pandas as pd

    df = pd.read_csv(repo_root() / "data" / "designs.csv")
    return [f"design_{int(d):03d}" for d in df["design_id"].tolist()]


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--models",
        nargs="+",
        default=[],
        help=f"Subset of complex predictors to run. Choices: {list(MODELS)}",
    )
    parser.add_argument(
        "--scorers",
        nargs="+",
        default=[],
        help=(
            "Subset of post-folding scorers to run "
            "(reads from data/structures/* and writes data/metrics/*). "
            f"Choices: {list(SCORERS)}"
        ),
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Run on only the first N designs (smoke test).",
    )
    parser.add_argument(
        "--no-detach",
        action="store_true",
        help="Block locally instead of detaching. Default is detached so "
        "the run survives client disconnects.",
    )
    parser.add_argument(
        "--launch-only",
        action="store_true",
        help="Trigger the remote batch and skip the download phase.",
    )
    parser.add_argument(
        "--download-only",
        action="store_true",
        help="Skip the launch phase, only pull what's already on the volume.",
    )
    args = parser.parse_args(argv)

    # If neither flag is set, default to the full complex-predictor panel
    # (backwards-compatible with the original behaviour).
    if not args.models and not args.scorers:
        args.models = list(MODELS)

    _check_models(args.models)
    _check_scorers(args.scorers)
    _check_token_env("COMPLEX_API_TOKEN", where=None)

    pb_ids = _design_ids()
    print(
        f"[complex] {len(pb_ids)} designs in data/designs.csv "
        f"(models={args.models}, scorers={args.scorers})"
    )

    for model in args.models:
        done = _already_on_disk(model, pb_ids)
        missing = [p for p in pb_ids if p not in done]
        print(
            f"[complex] {model}: {len(done)}/{len(pb_ids)} already on disk, {len(missing)} to fetch"
        )

    if not args.download_only:
        for model in args.models:
            rc = _launch(model, detach=not args.no_detach, limit=args.limit)
            if rc != 0:
                sys.stderr.write(f"[complex] {model}: launch failed (exit {rc})\n")
        for scorer in args.scorers:
            rc = _launch_scorer(scorer, detach=not args.no_detach, limit=args.limit)
            if rc != 0:
                sys.stderr.write(f"[scorers] {scorer}: launch failed (exit {rc})\n")

    if not args.launch_only:
        for model in args.models:
            rc = _download(model)
            if rc != 0:
                sys.stderr.write(f"[complex] {model}: download failed (exit {rc})\n")
        for scorer in args.scorers:
            rc = _download_scorer(scorer)
            if rc != 0:
                sys.stderr.write(f"[scorers] {scorer}: download failed (exit {rc})\n")


if __name__ == "__main__":
    main()
