"""Emit Protenix-v2 JSON inputs, one per binder × TREM2 target complex.

Each design becomes a 2-chain complex (chain A = binder, chain B = TREM2
construct). The target chain is identical across all 141 inputs.

Run:
    uv run python scripts/folding/build_input.py
    uv run python scripts/folding/build_input.py --n-recycles 4 --n-samples 5

Output:
    data/folding/protenix_v2/inputs/design_NNN.json
    data/folding/protenix_v2/inputs/_manifest.csv   # one row per design

The JSON schema follows
https://github.com/bytedance/Protenix/blob/main/docs/infer_json_format.md
- top level is a list of one dict (Protenix expects a list even for a single job)
- dict has ``name`` and ``sequences``
- each entry of ``sequences`` is ``{"proteinChain": {"sequence": ..., "count": 1}}``

The CLI flags ``--n-recycles`` and ``--n-samples`` are NOT written into the
JSON (Protenix takes them as CLI flags ``-c`` / ``-p``); they are recorded in
the manifest so the Modal runner picks them up consistently.
"""
from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path

import typer
from loguru import logger

from scripts.utils import repo_root

app = typer.Typer(add_completion=False, no_args_is_help=False)


@dataclass(frozen=True)
class Design:
    """One row from designs.fasta with header ``>{id}|{name}|{team}|{method}``."""

    design_id: int
    name: str
    team: str
    method: str
    sequence: str


def _parse_fasta(path: Path) -> list[tuple[str, str]]:
    """Bare-bones FASTA reader (no Biopython dep required)."""
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks)))
            header = line[1:]
            chunks = []
        else:
            chunks.append(line)
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def load_designs_fasta(path: Path) -> list[Design]:
    out: list[Design] = []
    for header, seq in _parse_fasta(path):
        parts = header.split("|")
        if len(parts) < 4:
            raise ValueError(
                f"Expected FASTA header '>{{id}}|{{name}}|{{team}}|{{method}}', got: {header}"
            )
        out.append(
            Design(
                design_id=int(parts[0]),
                name=parts[1],
                team=parts[2],
                method=parts[3],
                sequence=seq.upper(),
            )
        )
    return out


def load_target_sequence(path: Path) -> str:
    records = _parse_fasta(path)
    if len(records) != 1:
        raise ValueError(f"Expected exactly one record in {path}, got {len(records)}")
    return records[0][1].upper()


def build_payload(design: Design, target_sequence: str) -> list[dict]:
    """One Protenix predict job: binder + TREM2 as a 2-chain complex.

    Returns a list (Protenix expects a top-level list, even for a single job).
    """
    return [
        {
            "name": f"design_{design.design_id:03d}",
            "covalent_bonds": [],
            "sequences": [
                {
                    "proteinChain": {
                        "sequence": design.sequence,
                        "count": 1,
                        "modifications": [],
                    }
                },
                {
                    "proteinChain": {
                        "sequence": target_sequence,
                        "count": 1,
                        "modifications": [],
                    }
                },
            ],
        }
    ]


@app.command()
def main(
    n_recycles: int = typer.Option(3, "--n-recycles", help="Protenix -c (Pairformer cycles)."),
    n_samples: int = typer.Option(5, "--n-samples", help="Protenix -p (diffusion samples per seed)."),
    designs_fasta: Path = typer.Option(None, "--designs", help="Override path to designs.fasta."),
    target_fasta: Path = typer.Option(None, "--target", help="Override path to TREM2 construct FASTA."),
    out_dir: Path = typer.Option(None, "--out-dir", help="Override JSON output dir."),
) -> None:
    """Write one Protenix-v2 JSON per design plus a manifest CSV."""
    root = repo_root()
    designs_path = designs_fasta or root / "data" / "designs.fasta"
    target_path = target_fasta or root / "data" / "target" / "trem2_construct.fasta"
    out_path = out_dir or root / "data" / "folding" / "protenix_v2" / "inputs"
    out_path.mkdir(parents=True, exist_ok=True)

    target = load_target_sequence(target_path)
    designs = load_designs_fasta(designs_path)
    logger.info(f"Loaded {len(designs)} designs and {len(target)}-aa target from {target_path.name}")

    manifest_rows: list[dict] = []
    for d in designs:
        payload = build_payload(d, target)
        json_path = out_path / f"design_{d.design_id:03d}.json"
        json_path.write_text(json.dumps(payload, indent=2))
        manifest_rows.append(
            {
                "design_id": d.design_id,
                "name": d.name,
                "team": d.team,
                "method": d.method,
                "json_path": str(json_path.relative_to(root)),
                "binder_len": len(d.sequence),
                "target_len": len(target),
                "n_recycles": n_recycles,
                "n_samples": n_samples,
            }
        )

    manifest = out_path / "_manifest.csv"
    with manifest.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(manifest_rows[0].keys()))
        writer.writeheader()
        writer.writerows(manifest_rows)
    logger.success(f"Wrote {len(manifest_rows)} JSONs to {out_path} (manifest: {manifest.name})")


if __name__ == "__main__":
    app()
