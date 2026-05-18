# `analyses/` — one subdir per analysis

This is where collaborators add their analyses. **Each subdir owns its
outputs.** No cross-writing into other subdirs, no overwriting the
canonical CSV.

## Layout

```
analyses/
├── _template/              ← copy this to start a new analysis
├── human_vs_agent/         canonical cohort stats
├── leaderboard/            top-N tables
├── sequence_diversity/     within-team identity, k-mer entropy
├── methods/                design-method × outcome cross-tabs
└── <your_handle>/          ← your analysis goes here
```

## Adding yours — the full loop

```bash
# one-time
curl https://mise.run | sh         # install mise
mise install                       # python 3.11 + uv
mise run setup                     # install Python deps

# every analysis
cp -r analyses/_template analyses/<your_handle>
$EDITOR analyses/<your_handle>/main.py
$EDITOR mise.toml                  # add [tasks."analysis:<your_handle>"] block
mise run analysis:<your_handle>    # run it
mise run analysis:all              # confirm nothing else broke
```

Full pattern walkthrough in [`../docs/ANALYSES.md`](../docs/ANALYSES.md).
Conventions in [`../CONTRIBUTING.md`](../CONTRIBUTING.md).

## The data contract

**Every analysis reads `data/designs.csv` via `scripts.utils.load_designs()`.**
Never a hard-coded path. Never `data/raw_lab/bli_results.csv` for headline
KDs — `designs.csv` already has the canonical aggregate.

```python
from scripts.utils import load_designs

df = load_designs()                       # 141 × 123, every column
df = load_designs(only_screened=True)     # 100 designs sent to BLI
df = load_designs(only_hits=True)         # 37 binders
df = load_designs(only_human=True)        # 81 human-cohort designs
```

Column reference: [`../docs/DATA.md`](../docs/DATA.md).

## What lives in each analysis dir

| file | purpose | required? |
|---|---|---|
| `README.md`           | What the analysis answers, in <300 words   | yes |
| `main.py`             | Entry point — re-run with `mise run analysis:<name>` | yes |
| `report.md`           | Auto-generated human summary               | yes (auto) |
| `summary.json`        | Auto-generated headline numbers            | yes (auto) |
| `*.csv` / `*.parquet` | Derived tables this analysis produces      | optional |
| `notebooks/`          | Exploratory notebooks (load-bearing only)  | optional |
| `scratch/`            | Gitignored playground                      | optional |

## Rules

- Don't write outside your own subdir. Reviewers will catch this.
- Don't modify `data/designs.csv`. Open an issue with the offending
  `design_id` instead.
- Don't add helper code that other analyses will need. Put it under
  `scripts/utils/` and import from there.
- Every analysis registers a `mise` task in `mise.toml` so anyone can
  reproduce it with `mise run analysis:<name>`.
