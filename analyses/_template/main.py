"""
TEMPLATE — copy this into analyses/<your_handle>/ and edit.

What this analysis answers (replace this line): one-sentence framing.

Run:
    mise run analysis:<your_handle>
    # or:
    uv run python analyses/<your_handle>/main.py

Inputs:
- data/designs.csv via scripts.utils.load_designs() — ALWAYS go through
  the loader, never a hard-coded path.

Outputs (written inside this folder only):
- analyses/<your_handle>/report.md   human-readable summary, ≤1 page
- analyses/<your_handle>/summary.json machine-readable headline numbers
"""
from __future__ import annotations

import json
from pathlib import Path

from scripts.utils import load_designs

OUT = Path(__file__).resolve().parent


def main() -> None:
    # 1. Load the canonical data. ALWAYS via load_designs() — never a
    #    hard-coded path. Filter at load time when you can.
    df = load_designs(only_screened=True)

    # 2. Compute headline numbers. (Cast through `astype("boolean")` so missing
    #    values don't trigger pandas' silent-downcast FutureWarning.)
    is_hit = df["is_hit"].astype("boolean").fillna(False)
    summary = {
        "n_screened": int(len(df)),
        "n_human": int(df["is_human"].sum()),
        "n_agent": int((~df["is_human"]).sum()),
        "hit_rate_overall": float(is_hit.mean()),
        "hit_rate_human": float(is_hit[df["is_human"]].mean()),
        "hit_rate_agent": float(is_hit[~df["is_human"]].mean()),
    }

    # 3. Write the machine-readable summary.
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    # 4. Write a human-readable report (markdown).
    lines = [
        "# <Replace with your analysis title>",
        "",
        f"**Headline:** {summary['hit_rate_overall']:.0%} overall hit rate "
        f"(human {summary['hit_rate_human']:.0%}, agent {summary['hit_rate_agent']:.0%}).",
        "",
        "## Method",
        "",
        "Briefly: what was computed, which filter was applied, which stat test was used.",
        "Cite the column names from `data/designs.csv` so the reader can re-derive.",
        "",
        "## Results",
        "",
        "```",
        json.dumps(summary, indent=2),
        "```",
        "",
        "## Caveats",
        "",
        "- Small N: cohort hit rates are computed on the 100 designs sent to the BLI assay.",
        "- The top-100 ipSAE cherry-pick is an in-silico selection; the 41 designs not screened",
        "  are absent from this analysis.",
        "- This template's numbers are illustrative — replace with the test that actually answers",
        "  your question.",
        "",
    ]
    (OUT / "report.md").write_text("\n".join(lines))

    print(f"[template] wrote {OUT/'report.md'} and {OUT/'summary.json'}")


if __name__ == "__main__":
    main()
