# Figures

Two figure surfaces live in this repo. Different rules for each.

## Where figures live

| path | what's there |
|---|---|
| `figures/blog/`        | The 7 hand-authored HTML figures from the blog post. Static SVG/CSS, no Python. Update by hand when numbers change. See `figures/blog/README.md`. |
| `figures/paper/`       | Matplotlib-rendered figures for the manuscript. Reproducible from `data/designs.csv`. |
| `figures/exploration/` | Gitignored scratch. Sketches, drafts, dead ends. |

## Paper figures — the renderer contract

```python
import matplotlib.pyplot as plt
from scripts.plotting import apply_theme, save_figure

apply_theme()              # loads theme/matplotlibrc

fig, ax = plt.subplots(figsize=(3.5, 2.6))   # default column width
# ... plot ...
save_figure(fig, "paper/fig1_hit_rate")      # writes png/pdf/svg under figures/paper/
```

`apply_theme()` is idempotent. Call it once at the top. The matplotlibrc
sets fonts, sizes, ticks, grid, spines. Don't override manually.

## Colours: only from the palette

```python
import json
from scripts.utils.load_data import repo_root

with open(repo_root() / "theme" / "palettes.json") as f:
    pal = json.load(f)

cohort_colours = pal["agent_vs_human"]   # {"human": "#...", "agent": "#..."}
```

For the cohort split, **always** use `pal["agent_vs_human"]`. Same two
hues on every figure. Third or fourth colour? Take it from the named
brand palette in `palettes.json`. Don't pick hex codes by eye.

Sequential data (binding-strength bins, replicate counts) → named
sequential palette. If your data has no order, don't use a sequential
map. `viridis` is acceptable only after you've exhausted the named palette.

**Don't:** rainbow, jet, 3D bars, gridlines on top of data, drop shadows,
gradient fills.

## Fonts

| use | font |
|---|---|
| Display / panel titles | GT Pressura Extended (Regular) |
| Body / axis labels     | Geist |
| Numbers / data readouts| Geist Mono / IBM Plex Mono |

If the fonts aren't installed locally, `apply_theme()` falls back to
DejaVu Sans / DejaVu Sans Mono. **Don't substitute Helvetica or Arial
manually** — that breaks consistency across figures. See
[`../theme/fonts/README.md`](../theme/fonts/README.md) for install
instructions.

## Sizes

| target | width | height |
|---|---|---|
| Single column        | 3.5 in | 2.6 in |
| 1.5-column           | 5.0 in | 3.4 in |
| Full page width      | 7.2 in | 4.8 in |

Stick to these. Inconsistent figure widths read as sloppy at a glance.

## Formats

`save_figure()` writes PNG (300 dpi), PDF (vector), SVG (vector). The
paper LaTeX picks up the PDF; the blog picks up the SVG; the PNG is for
previews and Slack.

Don't commit raster images > 5 MB. Vectorise instead.

## Axes and labels

- Spines: bottom + left only. Top + right hidden. (matplotlibrc does this.)
- Ticks: outward, 4 pt major / 2 pt minor.
- Axis labels: sentence case, units in parens — `Affinity (nM)`,
  `Sequence length (aa)`.
- Legends: no box, no border, inside the data area unless that blocks the data.

## Statistical annotations

When you draw a p-value bracket, put both the p and the effect size:

```
**  Δmedian = +0.42, p = 0.018
```

Never just `**` or just `p<0.05`. The reader needs to know the effect.

## Caveats panels

If a figure has a caveat (small N for a per-team comparison, KD missing
for one hit), put the n directly under each bar/box. Don't hide it in
the caption.

## File naming

```
figures/paper/fig<N>_<short_name>.{png,pdf,svg}
figures/paper/fig<N>S<L>_<short_name>.{png,pdf,svg}   # SI figures
figures/paper/sx_<short_name>.{png,pdf,svg}           # extended panels
figures/exploration/<whatever>.{png,pdf,svg}          # gitignored scratch
```

`<short_name>` is `snake_case`, ≤ 4 words.

## Blog figures — update by hand

The 7 figures in `figures/blog/` (`article-graphical-abstract.html`,
`article-fig1-study-design.html`, …) are hand-authored standalone HTML
with inline SVG. **There is no Python renderer.** Open one in a browser
and you see the figure exactly as it appears in the blog post.

When `data/designs.csv` changes — and therefore the headline numbers
change — you edit the HTML by hand. `figures/blog/README.md` lists every
number baked into each figure with the column / `summary.json` it should
match.

## Examples

The four canonical analyses each produce a publication figure. Read
`analyses/leaderboard/` for the simplest path; `analyses/human_vs_agent/`
for a figure with statistical annotations.

For multipanel figures use `scripts/plotting/_common.py` (extend it if
you need a new layout helper — PR the helper, don't fork it).
