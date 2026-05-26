"""
Sequence diversity analysis for the muni x Adaptyv TREM2 hackathon.

Outputs:
- data/processed/within_team_identity.csv         (per-team pairwise identity)
- data/processed/kmer_entropy_per_design.csv      (k=3 Shannon entropy per design)
- analyses/sequence_diversity/cohort_identity.json (human vs agent diversity)
- analyses/sequence_diversity/report.md
- data/positive_controls/known_binders.fasta + README.md (placeholders if seq lookup fails)
- data/processed/identity_to_known_binders.csv (only if real positive controls found)
"""
from __future__ import annotations

import json
import math
import random
from collections import Counter
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

try:
    import Levenshtein
    HAVE_LEV = True
except ImportError:
    HAVE_LEV = False


from scripts.utils import load_designs, repo_root

ROOT = repo_root()
POS = ROOT / "data/controls"
OUT = Path(__file__).resolve().parent

K = 3
RNG = random.Random(0)


def identity_pair(a: str, b: str) -> float:
    """Levenshtein-distance-based identity, normalized by the longer sequence length."""
    if not a or not b:
        return float("nan")
    if HAVE_LEV:
        d = Levenshtein.distance(a, b)
    else:
        # Fallback: trivial dynamic programming - keep it short for the rare path.
        n, m = len(a), len(b)
        dp = list(range(m + 1))
        for i, ca in enumerate(a, 1):
            new = [i] + [0] * m
            for j, cb in enumerate(b, 1):
                new[j] = min(new[j - 1] + 1, dp[j] + 1, dp[j - 1] + (ca != cb))
            dp = new
        d = dp[m]
    L = max(len(a), len(b))
    return 1.0 - d / L


def kmer_entropy(seq: str, k: int = K) -> float:
    if len(seq) < k:
        return float("nan")
    kmers = [seq[i:i + k] for i in range(len(seq) - k + 1)]
    cnt = Counter(kmers)
    total = sum(cnt.values())
    return -sum((c / total) * math.log2(c / total) for c in cnt.values())


def main():
    POS.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)

    designs = load_designs()

    # --- 1. Within-team pairwise identity ------------------------------------
    rows = []
    for team, g in designs.groupby("team"):
        seqs = g.sequence.tolist()
        if len(seqs) < 2:
            rows.append({"team": team, "n_pairs": 0, "mean_identity": float("nan"),
                         "median_identity": float("nan"), "min_identity": float("nan"),
                         "max_identity": float("nan")})
            continue
        ids = [identity_pair(a, b) for a, b in combinations(seqs, 2)]
        rows.append({
            "team": team, "n_pairs": len(ids),
            "mean_identity": float(np.mean(ids)),
            "median_identity": float(np.median(ids)),
            "min_identity": float(np.min(ids)),
            "max_identity": float(np.max(ids)),
        })
    within = pd.DataFrame(rows).sort_values("mean_identity", ascending=False)
    within.to_csv(OUT / "within_team_identity.csv", index=False)

    # --- 2. Cohort pairwise identity (200 random pairs each) ------------------
    def sample_pairs(seqs: list[str], n: int) -> list[float]:
        seen = set()
        ids = []
        attempts = 0
        while len(ids) < n and attempts < n * 20:
            attempts += 1
            i, j = RNG.sample(range(len(seqs)), 2)
            key = (min(i, j), max(i, j))
            if key in seen:
                continue
            seen.add(key)
            ids.append(identity_pair(seqs[i], seqs[j]))
        return ids

    human_seqs = designs.loc[designs.is_human, "sequence"].tolist()
    agent_seqs = designs.loc[~designs.is_human, "sequence"].tolist()
    human_ids = sample_pairs(human_seqs, 200)
    agent_ids = sample_pairs(agent_seqs, 200)
    mw_id = stats.mannwhitneyu(human_ids, agent_ids, alternative="two-sided")

    # --- 3. K-mer entropy per design + cohort comparison ----------------------
    designs = designs.copy()
    designs["kmer_entropy_k3"] = designs.sequence.apply(kmer_entropy)
    designs[["design_id", "name", "team", "is_human", "sequence_length", "kmer_entropy_k3"]].to_csv(
        OUT / "kmer_entropy_per_design.csv", index=False
    )
    h_ent = designs.loc[designs.is_human, "kmer_entropy_k3"].dropna().to_numpy()
    a_ent = designs.loc[~designs.is_human, "kmer_entropy_k3"].dropna().to_numpy()
    mw_ent = stats.mannwhitneyu(h_ent, a_ent, alternative="two-sided")

    cohort_summary = {
        "pairwise_identity": {
            "n_human_pairs": len(human_ids), "n_agent_pairs": len(agent_ids),
            "median_identity_human": float(np.median(human_ids)),
            "median_identity_agent": float(np.median(agent_ids)),
            "mean_identity_human": float(np.mean(human_ids)),
            "mean_identity_agent": float(np.mean(agent_ids)),
            "mannwhitney_U": float(mw_id.statistic), "p": float(mw_id.pvalue),
            "implementation": "python-Levenshtein" if HAVE_LEV else "DP fallback",
        },
        "kmer_entropy_k3": {
            "n_human": int(len(h_ent)), "n_agent": int(len(a_ent)),
            "median_human": float(np.median(h_ent)), "median_agent": float(np.median(a_ent)),
            "mean_human": float(np.mean(h_ent)), "mean_agent": float(np.mean(a_ent)),
            "mannwhitney_U": float(mw_ent.statistic), "p": float(mw_ent.pvalue),
        },
    }
    (OUT / "cohort_identity.json").write_text(json.dumps(cohort_summary, indent=2))

    # --- 4. Positive controls -------------------------------------------------
    # AL002 / VHB937 sequences are not freely re-distributable from the patents at
    # short notice. Write a placeholder FASTA so the pipeline is deterministic and
    # document the lookup status; identity-to-known-binders is left for a later pass.
    pos_fasta = POS / "known_binders.fasta"
    if not pos_fasta.exists():
        pos_fasta.write_text(
            ">AL002_VH_PLACEHOLDER\n"
            "ACTUAL_SEQUENCE_PENDING_PATENT_WO2019028346A1\n"
            ">VHB937_PLACEHOLDER\n"
            "ACTUAL_SEQUENCE_PENDING_PATENT_WO2022122788A2\n"
        )
    pos_readme = POS / "README.md"
    if not pos_readme.exists():
        pos_readme.write_text(
            "# Positive controls (TREM2 hackathon)\n\n"
            "Two clinical-stage anti-TREM2 antibodies are referenced as assay sanity "
            "checks: **AL002** (Alector, patent WO2019028346A1, failed Phase 2 2024) and "
            "**VHB937** (Novartis, patent WO2022122788A2, ongoing Phase 2).\n\n"
            "## Status\n\n"
            "`known_binders.fasta` currently contains placeholder records. The variable-"
            "domain sequences need to be extracted from the patent SEQ IDs by hand or "
            "from the Adaptyv lab provider; the identity-to-known-binders comparison is "
            "therefore **pending**. Once real sequences land, regenerate this FASTA and "
            "re-run `analyses/sequence_diversity/diversity.py`; the script will produce "
            "`data/processed/identity_to_known_binders.csv` automatically.\n"
        )

    pos_records = []
    if pos_fasta.exists():
        cur_name = None
        cur_seq = []
        for line in pos_fasta.read_text().splitlines():
            if line.startswith(">"):
                if cur_name:
                    pos_records.append((cur_name, "".join(cur_seq)))
                cur_name = line[1:].strip()
                cur_seq = []
            else:
                cur_seq.append(line.strip())
        if cur_name:
            pos_records.append((cur_name, "".join(cur_seq)))

    real_records = [
        (n, s) for n, s in pos_records
        if s and "PENDING" not in s and "PLACEHOLDER" not in s and "_" not in s
    ]
    if real_records:
        rows = []
        for design_id, name, seq in zip(designs.design_id, designs.name, designs.sequence):
            for ctrl_name, ctrl_seq in real_records:
                rows.append({
                    "design_id": int(design_id), "design_name": name,
                    "control_name": ctrl_name,
                    "identity": identity_pair(seq, ctrl_seq),
                })
        pd.DataFrame(rows).to_csv(OUT / "identity_to_known_binders.csv", index=False)
        pc_status = f"computed for {len(real_records)} controls"
    else:
        pc_status = "PENDING - placeholder sequences only"

    # --- 5. Markdown report ---------------------------------------------------
    md = []
    md.append("# Sequence diversity report\n")
    md.append("## Within-team pairwise identity (top 5 most similar / 5 most diverse)")
    md.append("```")
    md.append(within.head(5).to_string(index=False))
    md.append("...")
    md.append(within.tail(5).to_string(index=False))
    md.append("```\n")

    md.append("## Cohort-level diversity")
    md.append(f"- Human cohort: {len(human_ids)} sampled pairs, median pairwise identity {np.median(human_ids):.3f}, mean {np.mean(human_ids):.3f}")
    md.append(f"- Agent cohort: {len(agent_ids)} sampled pairs, median pairwise identity {np.median(agent_ids):.3f}, mean {np.mean(agent_ids):.3f}")
    md.append(f"- Mann-Whitney U={mw_id.statistic:.0f}, p={mw_id.pvalue:.4g}\n")

    md.append("## K-mer (k=3) Shannon entropy")
    md.append(f"- Human n={len(h_ent)}, median {np.median(h_ent):.3f}; agent n={len(a_ent)}, median {np.median(a_ent):.3f}")
    md.append(f"- Mann-Whitney U={mw_ent.statistic:.0f}, p={mw_ent.pvalue:.4g}\n")

    md.append("## Positive controls (AL002, VHB937)")
    md.append(f"- Status: {pc_status}.")
    md.append("- See `data/positive_controls/README.md` for the lookup status; once the patent SEQ IDs are imported the identity table will be regenerated automatically.\n")

    (OUT / "report.md").write_text("\n".join(md) + "\n")

    print(f"[diversity] {len(within)} teams; cohort identity p={mw_id.pvalue:.4g}, kmer entropy p={mw_ent.pvalue:.4g}")
    print(f"[diversity] positive controls: {pc_status}")


if __name__ == "__main__":
    main()
