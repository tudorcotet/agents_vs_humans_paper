# Contributing

The point of this repo is one paper with consistent numbers and a coherent
voice, not the maximum number of features in any single notebook. Read
the rules, then send a PR.

## In 60 seconds

```bash
# one-time setup
curl https://mise.run | sh         # or: brew install mise
git clone <repo> && cd agents_vs_humans_paper
git lfs install && git lfs pull    # ~162 MB of data behind LFS
mise install                       # installs python 3.11 + uv
mise run setup                     # installs Python deps via uv

# verify the existing pipeline runs
mise run analysis:all

# add yours
cp -r analyses/_template analyses/<yourname>
$EDITOR analyses/<yourname>/main.py
$EDITOR mise.toml                  # add a [tasks."analysis:<yourname>"] block
mise run analysis:<yourname>
mise run analysis:all              # confirm nothing else broke
```

That's the whole loop. The rest of this doc is the rules each step
encodes.

## Where to put your code

```
analyses/
├── _template/              ← copy me to start
├── human_vs_agent/         (existing canonical analyses)
├── leaderboard/
├── sequence_diversity/
├── methods/
└── <yourname>/             ← your code lives here, NOWHERE ELSE
    ├── README.md           what your analysis answers (≤300 words)
    ├── main.py             entry point, runs end-to-end
    ├── report.md           auto-generated human-readable summary
    ├── summary.json        auto-generated machine-readable headline numbers
    └── *.csv               any tables your analysis writes
```

**Hard rule:** an analysis only writes inside its own `analyses/<name>/`
subdirectory. Don't touch other people's folders, don't modify
`data/designs.csv`, don't drop output into `scripts/` or the repo root.

## Which data to read

**Always go through `scripts.utils.load_designs()` — never read a file
path directly.**

```python
from scripts.utils import load_designs, load_replicates

df = load_designs()                       # 141 rows, every column
df = load_designs(only_screened=True)     # 100 designs sent to BLI
df = load_designs(only_hits=True)         # 37 binders
df = load_designs(only_human=True)        # 81 human-cohort designs
df = load_designs(only_human=False)       # 60 agent-cohort designs

reps = load_replicates()                  # per-replicate BLI rows
reps = reps[reps.selected & ~reps.excluded & reps.fixed]   # canonical filter
```

The columns are documented in [`docs/DATA.md`](docs/DATA.md). The
loader caches and coerces booleans for you — see
[`scripts/utils/load_data.py`](scripts/utils/load_data.py).

Don't read:
- `data/raw_lab/bli_results.csv` directly for KDs — `designs.csv` already
  has the canonical aggregate (`kd_arith_mean_nM_all`, `pkd_arith_mean`).
- `data/proteinbase/sensorgrams/*.json` unless you're plotting raw curves.
  For per-design replicate counts use `pb_n_bli_curves` and `pb_n_spr_curves`
  in `designs.csv`.

## How to run an analysis

Every analysis is a `mise` task. Install `mise` once (`curl https://mise.run | sh`)
and then:

```bash
mise run analysis:<name>            # run one
mise run analysis:all               # run every canonical analysis
mise run setup                      # install deps (one-time after clone)
mise run lint                       # ruff check
mise run figures:blog               # serve the blog HTML figures locally
```

If you don't want `mise`, the `Makefile` mirrors the same targets:
`make analysis-all`, `make setup`, `make lint`, etc.

## Invariants — don't break these

These are the assumptions every analysis depends on.

1. **`data/designs.csv` is the source of truth.** 141 rows. `design_id`
   is the join key everywhere. New annotations ship as a *separate* CSV
   under your analysis dir, keyed on `design_id`. Don't overwrite the
   canonical file.
2. **Cohort comes from `is_human`.** Never parse team names. Marcel is a
   human team; claude-sonnet-4.6 is an agent team.
3. **Replicate filter is `selected & ~excluded & fixed`.** The per-design
   KD columns in `designs.csv` already apply this. Re-aggregate from
   `bli_replicates.csv` only if you want a different filter, and document
   the change.
4. **Pair effect size with every p-value.** Spearman ρ + bootstrap 95% CI.
   Mann-Whitney U + median difference. Fisher's exact + the raw 2×2.
5. **Palette from `theme/palettes.json`.** The cohort split is
   `pal["agent_vs_human"]`. Same two hues on every figure.
6. **Paper voice is clinical.** See
   [`docs/STYLE_GUIDE.md`](docs/STYLE_GUIDE.md). Blog and social copy
   live elsewhere and follow a different voice — don't mix the two.

## Wiring your analysis into `mise`

Open `mise.toml`, scroll to the bottom (the "ADD YOURS BELOW" comment),
and add a block:

```toml
[tasks."analysis:yourname"]
description = "One-sentence answer to: what does my analysis tell us?"
run = "uv run python analyses/yourname/main.py"
```

If your analysis is load-bearing for the paper, also add
`"analysis:yourname"` to the `depends = [...]` list of
`[tasks."analysis:all"]` so it runs whenever someone re-runs everything.

Now `mise run analysis:yourname` works for everyone.

## Code style

- Python 3.11+. `ruff check .` passes.
- `from __future__ import annotations` at the top.
- Module docstring: what the script does, what it writes, ≤6 lines.
- Functions: `snake_case`, type-hinted. Side effects in `main()`.
- Tables: `*.csv` for human eyes, `*.parquet` when types matter.
- No commits of `.DS_Store`, `__pycache__/`, or `figures/exploration/`.

## Data hygiene

- Don't edit `data/designs.csv` by hand. If a row is wrong, open an issue
  with the `design_id` and the corrected value; the regeneration pipeline
  fixes it upstream.
- All binary data lives behind Git LFS (see [`.gitattributes`](.gitattributes)
  and the `## Cloning` section of [`README.md`](README.md)). New binaries
  > 1 MB → add the pattern to LFS in your PR.
- Per-replicate data lives in `data/raw_lab/bli_replicates.csv`. Reach it
  via `from scripts.utils import load_replicates` and re-apply
  `selected & ~excluded & fixed` yourself.

## Reviewing PRs

Before approving, check:

- Reads from `load_designs()`, not a hard-coded path.
- Writes only inside its own `analyses/<name>/`.
- p-values are paired with effect sizes.
- Cohort split uses `is_human`.
- `report.md` is ≤1 page with a Caveats section.
- A `[tasks."analysis:<name>"]` block exists in `mise.toml`.
- `mise run analysis:all` still passes.

## Communication

- Open a GitHub issue for anything that needs discussion before code.
- @-mention a maintainer on a PR you'd like reviewed promptly.
- Anything that changes the canonical CSV or redefines a cohort needs a
  discussion first. Open an issue, not a PR.
