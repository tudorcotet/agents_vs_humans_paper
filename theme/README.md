# `theme/` — shared visual identity

Every figure in this repo loads its look from here. Don't pick HEX codes
ad hoc.

## Layout

```
theme/
├── palettes.json      # named colour palettes (load me, don't eyeball)
├── matplotlibrc       # matplotlib defaults (loaded by apply_theme())
├── mpl_theme.py       # programmatic theme — extra helpers
├── STYLE.md           # the long-form publication style guide
└── fonts/
    └── README.md      # how to install GT Pressura / Geist / IBM Plex Mono
```

## Usage

```python
from scripts.plotting import apply_theme, save_figure

apply_theme()                 # load matplotlibrc + register palettes
fig, ax = plt.subplots()
# ...
save_figure(fig, "paper/fig1_hit_rate")
```

## Palettes

`palettes.json` defines:

- `agent_vs_human` — two distinct hues for the cohort split.
- `binding_strength` — binder / weak / non-binder / no-expression.
- `method_family` — colour per design method family.
- `brand_sequential` — sequential ramp from light to brand cyan.
- Plus DSSP, secondary-structure, and a handful of named single colours.

Load it like this:

```python
import json
from scripts.utils.load_data import repo_root

with open(repo_root() / "theme" / "palettes.json") as f:
    pal = json.load(f)
colour_human = pal["agent_vs_human"]["human"]
```

## Fonts

GT Pressura Extended (display), Geist (body), Geist Mono / IBM Plex Mono
(numbers). Install instructions in `fonts/README.md`. If the fonts are
missing, `apply_theme()` falls back to DejaVu Sans gracefully.

## Adding to the palette

If you genuinely need a new colour family, open a PR that adds a named
entry to `palettes.json`. Don't fork the palette per analysis.
