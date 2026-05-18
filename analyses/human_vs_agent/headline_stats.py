"""
Headline statistical comparisons: human vs agent cohort on the TREM2 hackathon.

Tests run (all unadjusted; multiple-testing caveat called out in the report):
1. Hit rate x cohort (Fisher's exact, top-100 only).
2. Expression failure x cohort (Fisher's exact, top-100 only).
3. Mann-Whitney U on pkd, binders only.
4. Spearman rho + bootstrap CI between submitted_ipsae and pkd.
5. Mann-Whitney on sequence_length over the full 141.
6. Chi-square on design_method_family x cohort, top-100.
7. Per-team hit rate ranking.
"""
from __future__ import annotations

import json
import math
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore", category=FutureWarning)

from scripts.utils import load_designs, repo_root

ROOT = repo_root()
OUT = Path(__file__).resolve().parent

CAVEATS = """\
**Caveats**
- Selection bias: hit-rate comparisons are conditional on surviving the in-silico ipSAE triage. 41 designs never reached the wet lab. The cohorts had different selection rates: 65 human / 35 agent into the top-100 from 81 / 60 submitted.
- Small N: 65 human, 35 agent in the screened set. Statistical power is limited; report effect sizes alongside p-values.
- Multiple testing: results are exploratory; report unadjusted p-values, note this.
- Replicate noise: for designs with `weird_replicates_flag=True`, the per-design KD has wider effective uncertainty than the SE alone suggests.
"""


def boot_diff_ci(a: np.ndarray, b: np.ndarray, statistic, n: int = 1000, seed: int = 0):
    rng = np.random.default_rng(seed)
    a, b = np.asarray(a), np.asarray(b)
    diffs = np.empty(n)
    for i in range(n):
        sa = rng.choice(a, size=len(a), replace=True)
        sb = rng.choice(b, size=len(b), replace=True)
        diffs[i] = statistic(sa) - statistic(sb)
    lo, hi = np.percentile(diffs, [2.5, 97.5])
    return float(lo), float(hi), float(diffs.mean())


def boot_spearman_ci(x: np.ndarray, y: np.ndarray, n: int = 1000, seed: int = 0):
    rng = np.random.default_rng(seed)
    rhos = np.empty(n)
    idx = np.arange(len(x))
    for i in range(n):
        s = rng.choice(idx, size=len(idx), replace=True)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            rhos[i] = stats.spearmanr(x[s], y[s]).correlation
    rhos = rhos[~np.isnan(rhos)]
    lo, hi = np.percentile(rhos, [2.5, 97.5])
    return float(lo), float(hi)


def fmt_2x2(name: str, a: int, b: int, c: int, d: int, row_labels=("Yes", "No"), col_labels=("Human", "Agent")) -> str:
    table = (
        f"{name}\n"
        f"                   {col_labels[0]:>8s}   {col_labels[1]:>8s}\n"
        f"  {row_labels[0]:>15s} {a:>8d}   {b:>8d}\n"
        f"  {row_labels[1]:>15s} {c:>8d}   {d:>8d}\n"
    )
    return table


def main():
    full = load_designs()
    full = full.copy()
    full["pkd"] = -np.log10(full.kd_arith_mean_nM_all / 1e9)
    full["in_top100"] = full.binding_label.notna()  # has a wet-lab result row
    top100 = full[full.in_top100].copy()

    out = {}
    contingency_lines = []

    # 1. Hit rate x cohort
    h = top100[top100.is_human]
    a = top100[~top100.is_human]
    h_hits = int(h.is_hit.sum())
    h_miss = int(len(h) - h_hits)
    a_hits = int(a.is_hit.sum())
    a_miss = int(len(a) - a_hits)
    fisher = stats.fisher_exact([[h_hits, a_hits], [h_miss, a_miss]], alternative="two-sided")
    rate_h = h_hits / len(h)
    rate_a = a_hits / len(a)
    lo, hi, mean = boot_diff_ci(h.is_hit.astype(int).to_numpy(), a.is_hit.astype(int).to_numpy(), np.mean)
    out["hit_rate_x_cohort"] = {
        "n_human": len(h), "n_agent": len(a),
        "hits_human": h_hits, "hits_agent": a_hits,
        "rate_human": rate_h, "rate_agent": rate_a,
        "diff_human_minus_agent": rate_h - rate_a,
        "fisher_odds_ratio": float(fisher.statistic), "fisher_p_two_sided": float(fisher.pvalue),
        "boot_diff_mean": mean, "boot_diff_ci_2_5": lo, "boot_diff_ci_97_5": hi, "boot_n": 1000,
    }
    contingency_lines.append(
        fmt_2x2("[1] Hit rate x cohort (top-100)", h_hits, a_hits, h_miss, a_miss, ("Hit", "Miss"))
        + f"  Fisher OR={fisher.statistic:.3f}  p={fisher.pvalue:.4g}\n"
    )

    # 2. Expression failure x cohort
    h_nonex = int((~h.expressed.astype(bool)).sum())
    a_nonex = int((~a.expressed.astype(bool)).sum())
    h_ex = len(h) - h_nonex
    a_ex = len(a) - a_nonex
    fisher_ex = stats.fisher_exact([[h_nonex, a_nonex], [h_ex, a_ex]], alternative="two-sided")
    out["nonexpression_x_cohort"] = {
        "nonexpressed_human": h_nonex, "nonexpressed_agent": a_nonex,
        "expressed_human": h_ex, "expressed_agent": a_ex,
        "rate_nonexpression_human": h_nonex / len(h), "rate_nonexpression_agent": a_nonex / len(a),
        "fisher_odds_ratio": float(fisher_ex.statistic), "fisher_p_two_sided": float(fisher_ex.pvalue),
    }
    contingency_lines.append(
        fmt_2x2("[2] Non-expression x cohort (top-100)", h_nonex, a_nonex, h_ex, a_ex, ("NonExpr", "Expressed"))
        + f"  Fisher OR={fisher_ex.statistic:.3f}  p={fisher_ex.pvalue:.4g}\n"
    )

    # 3. Mann-Whitney on pkd, binders only
    binders = top100.dropna(subset=["pkd"]).copy()
    h_pkd = binders.loc[binders.is_human, "pkd"].to_numpy()
    a_pkd = binders.loc[~binders.is_human, "pkd"].to_numpy()
    mw = stats.mannwhitneyu(h_pkd, a_pkd, alternative="two-sided")
    out["mannwhitney_pkd_binders"] = {
        "n_human_binders": int(len(h_pkd)), "n_agent_binders": int(len(a_pkd)),
        "median_pkd_human": float(np.median(h_pkd)) if len(h_pkd) else None,
        "iqr_pkd_human": [float(np.percentile(h_pkd, 25)), float(np.percentile(h_pkd, 75))] if len(h_pkd) else None,
        "median_pkd_agent": float(np.median(a_pkd)) if len(a_pkd) else None,
        "iqr_pkd_agent": [float(np.percentile(a_pkd, 25)), float(np.percentile(a_pkd, 75))] if len(a_pkd) else None,
        "U": float(mw.statistic), "p_two_sided": float(mw.pvalue),
        "median_kd_nM_human": float(10 ** (9 - np.median(h_pkd))) if len(h_pkd) else None,
        "median_kd_nM_agent": float(10 ** (9 - np.median(a_pkd))) if len(a_pkd) else None,
    }

    # 4. Spearman rho + bootstrap CI: ipSAE vs pkd
    # submitted_ipsae is now populated for 136/141 designs (the 5 without
    # are bottom-ranked non-screened). Binders-only block uses the binder
    # subset; the censored block uses the full top-100 with non-binders
    # left-censored at 10 uM.
    pairs_b = binders.dropna(subset=["submitted_ipsae", "pkd"]).copy()
    if len(pairs_b) >= 3:
        rho_b = stats.spearmanr(pairs_b.submitted_ipsae, pairs_b.pkd)
        lo_b, hi_b = boot_spearman_ci(pairs_b.submitted_ipsae.to_numpy(), pairs_b.pkd.to_numpy())
        binders_block = {
            "n": int(len(pairs_b)),
            "rho": float(rho_b.correlation), "p": float(rho_b.pvalue),
            "ci_2_5": lo_b, "ci_97_5": hi_b,
        }
    else:
        rho_b = None
        lo_b = hi_b = None
        binders_block = {"n": int(len(pairs_b)), "note": "submitted_ipsae missing for too many rows; cannot compute rho"}

    full_top = top100.dropna(subset=["submitted_ipsae"]).copy()
    full_top_excl = full_top.dropna(subset=["pkd"]).copy()
    rho_excl = stats.spearmanr(full_top_excl.submitted_ipsae, full_top_excl.pkd) if len(full_top_excl) >= 3 else None

    floor_kd_nM = 10000.0  # 10 uM left-censoring floor
    if len(full_top) >= 3:
        full_top_cens = full_top.copy()
        full_top_cens["pkd_cens"] = full_top_cens.pkd.fillna(-math.log10(floor_kd_nM / 1e9))
        rho_cens = stats.spearmanr(full_top_cens.submitted_ipsae, full_top_cens.pkd_cens)
        lo_c, hi_c = boot_spearman_ci(full_top_cens.submitted_ipsae.to_numpy(), full_top_cens.pkd_cens.to_numpy())
        cens_block = {
            "n": int(len(full_top_cens)), "floor_kd_nM": floor_kd_nM,
            "rho": float(rho_cens.correlation), "p": float(rho_cens.pvalue),
            "ci_2_5": lo_c, "ci_97_5": hi_c,
        }
    else:
        rho_cens = None
        lo_c = hi_c = None
        cens_block = {"n": int(len(full_top)), "note": "submitted_ipsae missing for too many rows; cannot compute rho"}

    out["spearman_ipsae_vs_pkd"] = {
        "binders_only": binders_block,
        "top100_excluding_nonbinders": {
            "n": int(len(full_top_excl)),
            "rho": float(rho_excl.correlation) if rho_excl else None,
            "p": float(rho_excl.pvalue) if rho_excl else None,
        },
        "top100_left_censored_at_10uM": cens_block,
    }

    # 5. Mann-Whitney on sequence_length, full 141
    h_len = full.loc[full.is_human, "sequence_length"].to_numpy()
    a_len = full.loc[~full.is_human, "sequence_length"].to_numpy()
    mw_len = stats.mannwhitneyu(h_len, a_len, alternative="two-sided")
    out["mannwhitney_sequence_length_all141"] = {
        "n_human": int(len(h_len)), "n_agent": int(len(a_len)),
        "median_human": float(np.median(h_len)), "median_agent": float(np.median(a_len)),
        "U": float(mw_len.statistic), "p_two_sided": float(mw_len.pvalue),
    }

    # 6. Chi-square on design_method_family x cohort, top-100
    ct = pd.crosstab(top100.method_family, top100.is_human)
    chi2 = stats.chi2_contingency(ct)
    contrib = (ct - chi2.expected_freq) ** 2 / chi2.expected_freq
    over_under = []
    for fam in ct.index:
        for cohort_bool in ct.columns:
            c = float(contrib.loc[fam, cohort_bool])
            if c > 1.0:
                obs = int(ct.loc[fam, cohort_bool])
                exp = float(chi2.expected_freq[ct.index.get_loc(fam), list(ct.columns).index(cohort_bool)])
                direction = "over" if obs > exp else "under"
                over_under.append({
                    "method_family": fam,
                    "cohort": "human" if bool(cohort_bool) else "agent",
                    "observed": obs, "expected": exp,
                    "contribution": c, "direction": direction,
                })
    out["chi_square_method_family_x_cohort"] = {
        "n": int(ct.values.sum()),
        "chi2": float(chi2.statistic), "p": float(chi2.pvalue), "dof": int(chi2.dof),
        "table": ct.to_dict(),
        "high_contribution_cells": over_under,
    }
    contingency_lines.append("[6] Method family x cohort (top-100)\n" + ct.to_string() +
                             f"\n  chi2={chi2.statistic:.3f}  dof={chi2.dof}  p={chi2.pvalue:.4g}\n")

    # 7. Per-team hit rate
    team_table = (
        top100.groupby(["team", "is_human"]).agg(
            n_lab=("design_id", "size"),
            n_binders=("is_hit", lambda x: int(x.astype(bool).sum())),
            n_expressed=("expressed", lambda x: int(x.astype(bool).sum())),
        ).reset_index()
    )
    team_table["hit_rate"] = team_table.n_binders / team_table.n_lab
    team_table = team_table.sort_values("hit_rate", ascending=False).reset_index(drop=True)
    out["team_hit_rates"] = team_table.to_dict(orient="records")

    # ---- Persist outputs -----------------------------------------------------
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "headline_stats.json").write_text(json.dumps(out, indent=2, default=str))

    # contingency_tables.txt
    contingency_lines.append("[7] Per-team hit rate (top-100)\n" + team_table.to_string(index=False))
    (OUT / "contingency_tables.txt").write_text("\n".join(contingency_lines) + "\n")

    # markdown report
    md = []
    md.append("# Headline statistics - human vs agent (TREM2 hackathon)\n")
    md.append(CAVEATS + "\n")
    md.append("## 1. Hit rate x cohort (top-100)")
    md.append(f"- Human: {h_hits}/{len(h)} = {rate_h*100:.1f}% hit")
    md.append(f"- Agent: {a_hits}/{len(a)} = {rate_a*100:.1f}% hit")
    md.append(f"- Diff (human - agent) = {(rate_h-rate_a)*100:.1f} pp; bootstrap 95% CI [{lo*100:.1f}, {hi*100:.1f}] pp (n=1000).")
    md.append(f"- Fisher's exact two-sided: OR={fisher.statistic:.2f}, p={fisher.pvalue:.4g}\n")

    md.append("## 2. Expression failure x cohort (top-100)")
    md.append(f"- Human non-expression: {h_nonex}/{len(h)} = {h_nonex/len(h)*100:.1f}%")
    md.append(f"- Agent non-expression: {a_nonex}/{len(a)} = {a_nonex/len(a)*100:.1f}%")
    md.append(f"- Fisher's exact: OR={fisher_ex.statistic:.2f}, p={fisher_ex.pvalue:.4g}\n")

    md.append("## 3. pkd by cohort (binders only)")
    md.append(f"- Human binders n={len(h_pkd)}, median pkd={np.median(h_pkd):.3f} (IQR {np.percentile(h_pkd,25):.3f} - {np.percentile(h_pkd,75):.3f})")
    md.append(f"- Agent binders n={len(a_pkd)}, median pkd={np.median(a_pkd):.3f} (IQR {np.percentile(a_pkd,25):.3f} - {np.percentile(a_pkd,75):.3f})")
    md.append(f"- Mann-Whitney U={mw.statistic:.1f}, p={mw.pvalue:.4g}\n")

    md.append("## 4. Spearman rho: ipSAE (`submitted_ipsae`) vs pkd")
    if rho_b is not None:
        md.append(f"- Binders only (n={len(pairs_b)}): rho={rho_b.correlation:.3f}, p={rho_b.pvalue:.4g}, 95% CI [{lo_b:.3f}, {hi_b:.3f}]")
    else:
        md.append(f"- Binders only (n={len(pairs_b)}): submitted_ipsae not populated; rho cannot be computed.")
    if rho_excl is not None:
        md.append(f"- Top-100 excluding non-binders (n={len(full_top_excl)}): rho={rho_excl.correlation:.3f}, p={rho_excl.pvalue:.4g}")
    if rho_cens is not None:
        md.append(f"- Top-100 left-censored at 10 uM (n={len(full_top)}): rho={rho_cens.correlation:.3f}, p={rho_cens.pvalue:.4g}, 95% CI [{lo_c:.3f}, {hi_c:.3f}]\n")
    else:
        md.append(f"- Top-100 left-censored at 10 uM: submitted_ipsae not populated; rho cannot be computed.\n")

    md.append("## 5. Sequence length human vs agent (full 141)")
    md.append(f"- Human n={len(h_len)}, median={np.median(h_len):.0f}")
    md.append(f"- Agent n={len(a_len)}, median={np.median(a_len):.0f}")
    md.append(f"- Mann-Whitney U={mw_len.statistic:.1f}, p={mw_len.pvalue:.4g}\n")

    md.append("## 6. Method family x cohort (top-100)")
    md.append("```")
    md.append(ct.to_string())
    md.append("```")
    md.append(f"chi2={chi2.statistic:.3f}, dof={chi2.dof}, p={chi2.pvalue:.4g}.")
    if over_under:
        md.append("Cells with chi-square contribution > 1.0:")
        for c in sorted(over_under, key=lambda d: -d["contribution"]):
            md.append(f"  - {c['method_family']} x {c['cohort']}: observed {c['observed']} vs expected {c['expected']:.1f} ({c['direction']}-represented, contribution {c['contribution']:.2f})")
    md.append("")

    md.append("## 7. Per-team hit rate (lab-tested designs only)")
    md.append("```")
    md.append(team_table.to_string(index=False))
    md.append("```\n")

    (OUT / "headline_stats.md").write_text("\n".join(md) + "\n")

    # ---- Console summary -----------------------------------------------------
    print("[stats] Hit rate human vs agent:")
    print(f"        human {h_hits}/{len(h)} = {rate_h:.3f}, agent {a_hits}/{len(a)} = {rate_a:.3f}")
    print(f"        Fisher OR={fisher.statistic:.3f}, p={fisher.pvalue:.4g}")
    print(f"[stats] pkd Mann-Whitney p={mw.pvalue:.4g}; med_h={np.median(h_pkd):.3f}, med_a={np.median(a_pkd):.3f}")
    if rho_b is not None:
        print(f"[stats] ipSAE vs pkd Spearman (binders): rho={rho_b.correlation:.3f}, p={rho_b.pvalue:.4g}, CI [{lo_b:.3f}, {hi_b:.3f}]")
    else:
        print(f"[stats] ipSAE vs pkd Spearman: submitted_ipsae is empty in all_designs.parquet; rho not computable.")
    print(f"[stats] outputs written to {OUT}")


if __name__ == "__main__":
    main()
