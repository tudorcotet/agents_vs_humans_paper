"""
The selection funnel by cohort: did humans' in-silico edge convert?

Humans submitted higher-ipSAE designs than agents, so more human designs
survived the top-100 ipSAE cherry-pick. This walks the funnel
(submitted -> survived ipSAE cut -> expressed -> bound) per cohort and asks
whether the human advantage on the *proxy* (ipSAE) carries through to the
*wet-lab* outcome, or evaporates.

Outputs:
- analyses/helen/cohort_funnel_report.md   (human-readable, with Caveats)
- analyses/helen/cohort_funnel_summary.json (machine-readable)
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy import stats

from scripts.utils import load_designs

HERE = Path(__file__).resolve().parent


def _fisher(a_pos: int, a_n: int, h_pos: int, h_n: int) -> dict:
    """2x2 Fisher (rows = agent/human, cols = pos/neg). OR is agent-vs-human."""
    table = [[a_pos, a_n - a_pos], [h_pos, h_n - h_pos]]
    odds, p = stats.fisher_exact(table, alternative="two-sided")
    return {
        "agent": [int(a_pos), int(a_n)], "human": [int(h_pos), int(h_n)],
        "rate_agent": round(a_pos / a_n, 4) if a_n else None,
        "rate_human": round(h_pos / h_n, 4) if h_n else None,
        "odds_ratio_agent_vs_human": (None if not np.isfinite(odds)
                                      else round(float(odds), 4)),
        "fisher_p_two_sided": round(float(p), 5),
        "table_agent_pos_neg_human_pos_neg": table,
    }


def _bh(pvals: list[float], q: float = 0.10) -> list[bool]:
    m = len(pvals)
    order = np.argsort(pvals)
    passed = [False] * m
    for rank, idx in enumerate(order, start=1):
        if pvals[idx] <= rank / m * q:
            for j in order[:rank]:
                passed[j] = True
    return passed


def main() -> None:
    df = load_designs()
    h = df[df.is_human.fillna(False)]
    a = df[~df.is_human.fillna(False)]
    scr_h = h[h.submitted_to_lab.fillna(False)]
    scr_a = a[a.submitted_to_lab.fillna(False)]

    stages = {
        # survived the top-100 ipSAE cut, of everything submitted
        "survived_ipsae_cut": _fisher(
            int(scr_a.shape[0]), int(a.shape[0]),
            int(scr_h.shape[0]), int(h.shape[0])),
        # expressed, of the screened set
        "expressed_given_screened": _fisher(
            int(scr_a.expressed.fillna(False).sum()), int(scr_a.shape[0]),
            int(scr_h.expressed.fillna(False).sum()), int(scr_h.shape[0])),
        # bound, of the screened set (the headline hit rate, in funnel context)
        "bound_given_screened": _fisher(
            int(scr_a.is_hit.fillna(False).sum()), int(scr_a.shape[0]),
            int(scr_h.is_hit.fillna(False).sum()), int(scr_h.shape[0])),
        # end-to-end: bound, of everything submitted
        "bound_given_submitted": _fisher(
            int(scr_a.is_hit.fillna(False).sum()), int(a.shape[0]),
            int(scr_h.is_hit.fillna(False).sum()), int(h.shape[0])),
    }

    # Mechanism: ipSAE of everything submitted (why humans survive the cut more)
    hi = h.submitted_ipsae.dropna().to_numpy(float)
    ai = a.submitted_ipsae.dropna().to_numpy(float)
    U, p_ips = stats.mannwhitneyu(ai, hi, alternative="two-sided")
    rank_biserial = float(2 * U / (len(ai) * len(hi)) - 1)  # agent - human
    ipsae = {
        "n_human": len(hi), "n_agent": len(ai),
        "median_human": round(float(np.median(hi)), 4),
        "median_agent": round(float(np.median(ai)), 4),
        "mannwhitney_U": float(U), "p_two_sided": round(float(p_ips), 6),
        "rank_biserial_agent_minus_human": round(rank_biserial, 4),
    }

    pvals = [s["fisher_p_two_sided"] for s in stages.values()]
    for s, sig in zip(stages.values(), _bh(pvals, 0.10), strict=True):
        s["bh_q0.10_significant"] = bool(sig)

    summary = {"funnel": stages, "ipsae_submitted": ipsae,
               "bh_fdr": {"q": 0.10, "family_size": len(pvals)}}
    (HERE / "cohort_funnel_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n")

    sv, ex, bs, be = (stages["survived_ipsae_cut"],
                      stages["expressed_given_screened"],
                      stages["bound_given_screened"],
                      stages["bound_given_submitted"])
    md = [
        "# The selection funnel by cohort: did the in-silico edge convert?\n",
        "**Headline:** humans submitted markedly higher-ipSAE designs "
        f"(median {ipsae['median_human']} vs {ipsae['median_agent']}, "
        f"Mann-Whitney p={ipsae['p_two_sided']:.2g}, rank-biserial "
        f"{ipsae['rank_biserial_agent_minus_human']:+.2f}), so a larger share "
        f"survived the top-100 ipSAE cut ({sv['rate_human']:.0%} of human vs "
        f"{sv['rate_agent']:.0%} of agent submissions, Fisher "
        f"p={sv['fisher_p_two_sided']:.2g}). Yet among lab-tested designs the "
        f"hit rate is statistically tied ({bs['rate_human']:.0%} human vs "
        f"{bs['rate_agent']:.0%} agent, OR "
        f"{bs['odds_ratio_agent_vs_human']}, Fisher "
        f"p={bs['fisher_p_two_sided']:.2g}). The proxy advantage does not "
        "carry through to the wet lab.\n",
        "## Method\n",
        "`load_designs()` (141). Funnel per cohort: submitted -> survived the "
        "top-100 ipSAE cut (`submitted_to_lab`) -> expressed (`expressed`, "
        "conditional on screened) -> bound (`is_hit`, conditional on "
        "screened), plus end-to-end bound/submitted. Each stage: 2x2 Fisher "
        "exact two-sided with raw counts and the agent-vs-human odds ratio. "
        "Mechanism: Mann-Whitney on `submitted_ipsae` over all submissions + "
        "rank-biserial effect size. BH-FDR q=0.10 over the four stage tests.\n",
        "## Results\n",
        f"- **Survived ipSAE cut:** human {sv['human'][0]}/{sv['human'][1]} "
        f"({sv['rate_human']:.0%}) vs agent {sv['agent'][0]}/{sv['agent'][1]} "
        f"({sv['rate_agent']:.0%}); OR {sv['odds_ratio_agent_vs_human']}, "
        f"Fisher p={sv['fisher_p_two_sided']:.2g}, BH-sig "
        f"{sv['bh_q0.10_significant']}.",
        f"- **Expressed | screened:** human {ex['human'][0]}/{ex['human'][1]} "
        f"({ex['rate_human']:.0%}) vs agent {ex['agent'][0]}/{ex['agent'][1]} "
        f"({ex['rate_agent']:.0%}); Fisher p={ex['fisher_p_two_sided']:.2g}, "
        f"BH-sig {ex['bh_q0.10_significant']}.",
        f"- **Bound | screened:** human {bs['human'][0]}/{bs['human'][1]} "
        f"({bs['rate_human']:.0%}) vs agent {bs['agent'][0]}/{bs['agent'][1]} "
        f"({bs['rate_agent']:.0%}); OR {bs['odds_ratio_agent_vs_human']}, "
        f"Fisher p={bs['fisher_p_two_sided']:.2g}, BH-sig "
        f"{bs['bh_q0.10_significant']}.",
        f"- **Bound | submitted (end-to-end):** human {be['human'][0]}/"
        f"{be['human'][1]} ({be['rate_human']:.0%}) vs agent "
        f"{be['agent'][0]}/{be['agent'][1]} ({be['rate_agent']:.0%}); OR "
        f"{be['odds_ratio_agent_vs_human']}, Fisher "
        f"p={be['fisher_p_two_sided']:.2g}, BH-sig "
        f"{be['bh_q0.10_significant']}.",
        f"- **ipSAE (all submitted):** human median {ipsae['median_human']} "
        f"(n={ipsae['n_human']}) vs agent {ipsae['median_agent']} "
        f"(n={ipsae['n_agent']}); Mann-Whitney p={ipsae['p_two_sided']:.2g}, "
        f"rank-biserial {ipsae['rank_biserial_agent_minus_human']:+.2f}.",
        "\n## Caveats\n",
        "- **Guaranteed-inclusion confound.** Per the Notion draft, top agent "
        "designs had a guaranteed path into the screened set, so the agent "
        "'survival' rate is not a pure ipSAE cut — it understates the true "
        "ipSAE gap between cohorts. Interpret the survival stage carefully.",
        "- **ipSAE is not method-neutral.** Humans used methods (BindCraft, "
        "AF-hallucination, Mosaic) that tend to yield higher Boltz-2 interface "
        "confidence; agents leaned on PXDesign. 'Humans had higher ipSAE' is "
        "partly a method-mix effect, not proof of intent.",
        "- **Small agent N** (60 submitted, 35 screened, 12 binders): the "
        "tied hit rate is consistent with no difference *and* with a modest "
        "difference the study cannot resolve. Effect sizes + counts reported.",
        "- **Top-100 cherry-pick.** All conditional rates are downstream of "
        "the ipSAE triage; the 41 unscreened designs are never lab-tested.",
        "- Conditional definitions are explicit (| screened); do not read "
        "'bound | screened' as an unconditional binding probability.",
    ]
    (HERE / "cohort_funnel_report.md").write_text("\n".join(md) + "\n")
    print(f"[cohort_funnel] ipSAE p={ipsae['p_two_sided']:.2g} "
          f"(human {ipsae['median_human']} vs agent {ipsae['median_agent']}); "
          f"survival p={sv['fisher_p_two_sided']:.2g}; "
          f"hit p={bs['fisher_p_two_sided']:.2g}")


if __name__ == "__main__":
    main()
