# `analyses/sequence_diversity/` — within-team & cohort diversity

How different are the sequences within a team? Between human and agent
cohorts?

## What it computes

1. **Within-team pairwise identity** — Levenshtein-based normalised
   identity across every pair within a team. Output: mean / median / min
   / max per team.
2. **Cohort pairwise identity** — 200 random pairs sampled from each
   cohort. Mann-Whitney U on the distributions.
3. **K-mer (k=3) Shannon entropy** — per-design k-mer diversity, then
   cohort comparison.
4. **Identity to known binders** — when the AL002 / VHB937 reference
   sequences land, compute identity from every design to each.

## Inputs

- `data/designs.csv` via `scripts.utils.load_designs()`.
- `data/controls/known_binders.fasta` (placeholder until patent SEQ IDs
  are extracted).

## Outputs

- `within_team_identity.csv` — per-team identity stats.
- `kmer_entropy_per_design.csv` — k-mer entropy column.
- `cohort_identity.json` — headline Mann-Whitney stats for cohort
  diversity.
- `report.md` — markdown summary.
- `identity_to_known_binders.csv` — only emitted when real positive-
  control sequences are present (i.e. not placeholders).

## Run

```bash
mise run analysis:diversity
```

## Notes

- Requires `python-Levenshtein` for speed; falls back to a slow DP
  implementation if missing.
- Pair sampling is seeded (`random.Random(0)`) so re-runs are
  deterministic.
