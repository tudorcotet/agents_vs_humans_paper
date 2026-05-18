"""
Build per-design leaderboards (top-N tables for figures).

Inputs:
- data/processed/all_designs_with_metrics_and_results.parquet

Outputs (analyses/leaderboard/):
- top20_overall.csv
- top10_human.csv
- top10_agent.csv
- method_winners.csv          (best design per design_method_normalized)
- team_winners.csv            (best design per team)
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from scripts.utils import load_designs, repo_root

ROOT = repo_root()
OUT = Path(__file__).resolve().parent

COLS = [
    "rank", "design_id", "name", "team", "is_human",
    "design_method", "design_method_normalized", "sequence_length",
    "kd_arith_mean_nM_all", "pkd",
    "n_replicates_pushed", "weird_replicates_flag",
    "submitted_ipsae", "binding_label", "binding_strength",
]


def rank_table(df: pd.DataFrame) -> pd.DataFrame:
    df = df.dropna(subset=["kd_arith_mean_nM_all"]).copy()
    df["pkd"] = -np.log10(df.kd_arith_mean_nM_all / 1e9)
    df = df.sort_values("kd_arith_mean_nM_all").reset_index(drop=True)
    df.insert(0, "rank", np.arange(1, len(df) + 1))
    return df[COLS]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    full = load_designs()

    overall = rank_table(full).head(20)
    overall.to_csv(OUT / "top20_overall.csv", index=False)

    human = rank_table(full[full.is_human]).head(10)
    human.to_csv(OUT / "top10_human.csv", index=False)

    agent = rank_table(full[~full.is_human]).head(10)
    agent.to_csv(OUT / "top10_agent.csv", index=False)

    # Best design per method family.
    binders = full.dropna(subset=["kd_arith_mean_nM_all"]).copy()
    binders["pkd"] = -np.log10(binders.kd_arith_mean_nM_all / 1e9)
    method_winners = (
        binders.sort_values(["design_method_normalized", "kd_arith_mean_nM_all"])
        .groupby("design_method_normalized", as_index=False)
        .first()
        .sort_values("kd_arith_mean_nM_all")
        .reset_index(drop=True)
    )
    method_winners.insert(0, "rank", np.arange(1, len(method_winners) + 1))
    method_winners[COLS].to_csv(OUT / "method_winners.csv", index=False)

    # Best design per team.
    team_winners = (
        binders.sort_values(["team", "kd_arith_mean_nM_all"])
        .groupby("team", as_index=False)
        .first()
        .sort_values("kd_arith_mean_nM_all")
        .reset_index(drop=True)
    )
    team_winners.insert(0, "rank", np.arange(1, len(team_winners) + 1))
    team_winners[COLS].to_csv(OUT / "team_winners.csv", index=False)

    # Console preview
    print(f"[leaderboard] top-20 overall written to {OUT/'top20_overall.csv'}")
    print(overall[["rank", "design_id", "name", "team", "is_human", "kd_arith_mean_nM_all", "pkd"]].head(5).to_string(index=False))
    print(f"[leaderboard] {len(method_winners)} method winners, {len(team_winners)} team winners")


if __name__ == "__main__":
    main()
