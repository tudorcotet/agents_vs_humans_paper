# `figures/blog/` — the seven blog post figures

Hand-authored inline-SVG / HTML/CSS figures used in the blog post. Each
file is **fully self-contained** — open it in a browser and it renders
with no build step.

There is **no Python renderer** for these. The numbers are baked into the
SVG paths and text nodes. If a data value changes, you edit the HTML by
hand. Deliberate tradeoff: one self-contained file per figure, no
toolchain needed to inspect it.

## The seven

| File | What it shows | Effort to update |
|---|---|---|
| `article-graphical-abstract.html`        | Four-card abstract: target, 141→100→89→37 cascade, top-5 leaderboard, agents-on-PXDesign monoculture | low |
| `article-fig1-study-design.html`          | Sankey funnel + per-team submission bars | medium (Sankey paths are bezier) |
| `article-fig2-headline-distributions.html`| Hit-rate stat card + human/agent pKD kernel densities | **high** (KDE paths are hundreds of coords) |
| `article-fig3-tool-monoculture.html`      | Per-agent tool allocation + cohort tool diversity | low |
| `article-fig4-winners.html`               | Scatter (hit rate × best KD) + per-team portfolio bars | medium |
| `article-fig-top5-binders.html`           | Leaderboard of the five tightest binders | low |
| `article-fig-target.html`                 | TREM2 target spec sheet + AL002 / VHB937 controls | low |

## Layout

```
figures/blog/
├── *.html              the 7 figures
├── fonts/              Geist-{Variable,Regular,Medium,Bold}.woff2, GT-Pressura-Extended-Regular.woff2
└── assets/
    ├── logos/          adaptyv-wordmark.svg, muni-wordmark.svg (+ dark variants)
    └── structures/     trem2-structure.png
```

## Source-of-truth numbers (regenerate before publishing)

When `data/designs.csv` changes, these numbers in the HTML need to track:

| number | columns / scripts |
|---|---|
| 141 / 100 / 37 / 11 / 52       | `load_designs().shape`, `submitted_to_lab`, `is_hit`, `binding_label` |
| 65 human / 35 agent screened   | `df[df.submitted_to_lab].groupby('is_human').size()` |
| 38.5% / 34.3% hit rate         | `analyses/human_vs_agent/summary.json` |
| Fisher OR=1.20, p=0.83         | `analyses/human_vs_agent/summary.json` |
| pKD medians 7.04 / 6.98        | `analyses/human_vs_agent/summary.json` |
| Top-5 KDs (1.11, 1.25, 1.91, 3.64, 6.85 nM) | `analyses/leaderboard/top20_overall.csv` |
| Agent tool shares (PXDesign 53%, …) | `df[~df.is_human].method_family.value_counts(normalize=True)` |

After re-running `mise run analysis:all`, scan each `summary.json` and
update the HTML where the numbers diverge. Keep the diff small — one
search-and-replace per number.

## Opening locally

```bash
open figures/blog/article-graphical-abstract.html
# (or any other; they're all standalone)
```

Or serve the directory to test cross-file font caching:

```bash
python -m http.server -d figures/blog 8080
# http://localhost:8080/article-graphical-abstract.html
```

## Don't

- Don't add a Plotly/D3 pipeline for these. They're hand-tuned for the
  blog; replacing with a generator would lose the typographic polish.
- Don't fork colours. The palette uses the brand tokens from
  `../../theme/palettes.json`; if you need a new hue, add it there
  *and* update the relevant CSS variables here.
- Don't commit the rendered PNG/SVG exports — let the HTML be the
  primary artifact.
