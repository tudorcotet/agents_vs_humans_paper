# `analyses/methods/` — design-method × outcome cross-tabs

Which design methods produce the highest hit rates? Which fail to
express? This analysis builds the per-method cross-tab.

## What it computes

For every `design_method_normalized`:

- `n_used_total` — total designs from this method (any cohort, any stage).
- `n_in_top100` — designs that made the BLI assay.
- `n_expressed` — designs that expressed.
- `n_binders` — designs labelled `binder` / `strong` / `medium` / `weak`.
- `hit_rate` — `n_binders / n_in_top100`.
- `hit_rate_among_expressed` — `n_binders / n_expressed`.
- `expression_rate` — `n_expressed / n_in_top100`.
- `median_kd_nM_binders` — median KD over binders.
- `n_human_teams` / `n_agent_teams` — distinct team counts using this
  method.

## Inputs

- `data/designs.csv` via `scripts.utils.load_designs()`.

## Outputs

- `method_outcome_xtab.csv` — the cross-tab.
- `report.md` — markdown summary, naming "winners" (hit-rate >
  overall + 1 SE) and "expression-failure-prone" methods (< 80%
  expression).

## Run

```bash
mise run analysis:methods
```

## Caveats

- Many methods are used by ≤3 designs. We mark a method as a "winner"
  only when `n_in_top100 ≥ 3`.
- `Foundry` is a wild-west bucket meaning "the agent chose its own
  tools." Don't read it as a single algorithm.
