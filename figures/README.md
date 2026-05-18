# `figures/`

Output directory for every rendered figure. **Code reads from `data/`,
writes here.**

## Layout

```
figures/
├── paper/         # ⭐ figures used in the manuscript (PNG + PDF + SVG)
└── exploration/   # gitignored scratch — sketches, drafts, dead ends
```

## Naming

```
figures/paper/fig<N>_<short_name>.{png,pdf,svg}
figures/paper/fig<N>S<L>_<short_name>.{png,pdf,svg}     # SI figures
figures/paper/sx_<short_name>.{png,pdf,svg}             # extended panels
```

`<short_name>` is `snake_case`, ≤ 4 words.

## How to render

Every figure script should look like:

```python
from scripts.plotting import apply_theme, save_figure
import matplotlib.pyplot as plt

apply_theme()
fig, ax = plt.subplots(figsize=(3.5, 2.6))
# ... plot ...
save_figure(fig, "paper/fig1_hit_rate")   # writes 3 files under figures/paper/
```

See [`../docs/FIGURES.md`](../docs/FIGURES.md) for the full style guide.
