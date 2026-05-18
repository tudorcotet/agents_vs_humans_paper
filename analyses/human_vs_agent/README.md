# `analyses/human_vs_agent/` — cohort comparison

The headline analysis of the paper: do AI agents and human teams produce
different binders against TREM2?

## What it answers

1. Hit rate × cohort (Fisher's exact, top-100 only).
2. Expression failure × cohort (Fisher's exact, top-100 only).
3. Mann-Whitney U on `pkd`, binders only.
4. Spearman ρ + bootstrap CI between `submitted_ipsae` and `pkd`.
5. Mann-Whitney on `sequence_length` over the full 141.
6. Chi-square on `method_family` × cohort, top-100.
7. Per-team hit-rate ranking.

## Inputs

- `data/designs.csv` via `scripts.utils.load_designs()`.

## Outputs

- `report.md` — human-readable markdown summary.
- `summary.json` — machine-readable headline numbers (used by figures).

## Run

```bash
mise run analysis:human-vs-agent
# or
uv run python analyses/human_vs_agent/headline_stats.py
```

## Caveats baked into the report

- **Selection bias**: hit-rate comparisons are conditional on having
  survived the in-silico ipSAE triage. 41 designs never made it to the
  wet lab.
- **Small N**: 65 human / 35 agent in the screened set. Effect sizes
  matter more than p-values.
- **Multiple testing**: unadjusted. We use BH-FDR when running >5 tests.
