"""
Method x outcome cross-tab analysis.

Produces:
- data/processed/method_outcome_xtab.csv (per design_method_normalized)
- analyses/methods/report.md             (winners / losers commentary)
"""
from __future__ import annotations

import math
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=FutureWarning)

from scripts.utils import load_designs, repo_root

ROOT = repo_root()
OUT = Path(__file__).resolve().parent


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    full = load_designs()

    full = full.copy()
    full["in_top100"] = full.binding_label.notna()

    rows = []
    for method, g in full.groupby("design_method_normalized", sort=False):
        n_used_total = len(g)
        in_top = g[g.in_top100]
        n_in_top100 = len(in_top)
        n_expressed = int(in_top.expressed.fillna(False).astype(bool).sum())
        n_binders = int(in_top.is_hit.fillna(False).astype(bool).sum())
        kd_med = float(np.median(in_top.kd_arith_mean_nM_all.dropna())) if in_top.kd_arith_mean_nM_all.notna().any() else float("nan")
        n_human_teams = int(g.loc[g.is_human, "team"].nunique())
        n_agent_teams = int(g.loc[~g.is_human, "team"].nunique())
        rows.append({
            "design_method_normalized": method,
            "n_used_total": n_used_total,
            "n_in_top100": n_in_top100,
            "n_expressed": n_expressed,
            "n_binders": n_binders,
            "hit_rate": (n_binders / n_in_top100) if n_in_top100 else float("nan"),
            "hit_rate_among_expressed": (n_binders / n_expressed) if n_expressed else float("nan"),
            "expression_rate": (n_expressed / n_in_top100) if n_in_top100 else float("nan"),
            "median_kd_nM_binders": kd_med,
            "n_human_teams": n_human_teams,
            "n_agent_teams": n_agent_teams,
        })
    xtab = pd.DataFrame(rows).sort_values("hit_rate", ascending=False).reset_index(drop=True)
    xtab.to_csv(OUT / "method_outcome_xtab.csv", index=False)

    # Identify "winning methods": hit rate > overall mean by >1 SE.
    overall_hit = (
        xtab.n_binders.sum() / xtab.n_in_top100.sum()
    ) if xtab.n_in_top100.sum() else float("nan")
    eligible = xtab[xtab.n_in_top100 >= 3].copy()
    se_per = np.sqrt(eligible.hit_rate * (1 - eligible.hit_rate) / eligible.n_in_top100.clip(lower=1))
    winners = eligible[eligible.hit_rate > overall_hit + se_per].sort_values("hit_rate", ascending=False)

    # "Expression-failure-prone" methods: expression rate < 80%.
    losers_expr = xtab[(xtab.n_in_top100 >= 3) & (xtab.expression_rate < 0.80)].sort_values("expression_rate")

    md = []
    md.append("# Method x outcome report\n")
    md.append(f"Overall lab hit rate: {overall_hit*100:.1f}% ({xtab.n_binders.sum()}/{xtab.n_in_top100.sum()}).\n")
    md.append("## Methods cross-tab (hit-rate sorted)")
    md.append("```")
    md.append(xtab.to_string(index=False))
    md.append("```\n")

    md.append("## Winning methods (hit rate > overall + 1 SE; n_in_top100 >= 3)")
    if len(winners):
        md.append("```")
        md.append(winners[["design_method_normalized", "n_in_top100", "n_binders", "hit_rate", "median_kd_nM_binders"]].to_string(index=False))
        md.append("```")
    else:
        md.append("- None passed the threshold.")
    md.append("")

    md.append("## Expression-failure-prone methods (expression rate < 80%; n_in_top100 >= 3)")
    if len(losers_expr):
        md.append("```")
        md.append(losers_expr[["design_method_normalized", "n_in_top100", "n_expressed", "expression_rate", "n_binders"]].to_string(index=False))
        md.append("```")
    else:
        md.append("- None passed the threshold.")
    md.append("")

    (OUT / "report.md").write_text("\n".join(md) + "\n")

    print(f"[methods] {len(xtab)} methods; overall hit rate {overall_hit:.3f}; {len(winners)} winners, {len(losers_expr)} expression-failure-prone.")
    print("[methods] top 3 by hit rate (n_in_top100>=3):")
    print(eligible.head(3)[["design_method_normalized", "n_in_top100", "n_binders", "hit_rate", "expression_rate"]].to_string(index=False))
    print("[methods] bottom 3 by expression rate (n_in_top100>=3):")
    print(eligible.sort_values("expression_rate").head(3)[["design_method_normalized", "n_in_top100", "n_expressed", "expression_rate"]].to_string(index=False))


if __name__ == "__main__":
    main()
