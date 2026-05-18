# Fonts

Three faces, three jobs. Don't substitute. Don't add a fourth face for
"variety".

## 1. GT Pressura Extended Regular — display only

Wide-display sans by Grilli Type. Used for hero titles, panel headings,
graphical-abstract section markers. **Never** for body, chart labels, or
UI copy.

| Property | Value |
|---|---|
| Weight in use | 400 (Regular) only |
| Tracking | Near zero — `letter-spacing: -1.0px` at 60px, `-0.6px` at 42px |
| Sizes | 60 / 50 / 42 px (`title-xl` / `title-lg` / `title-md`) |
| License | Paid (Grilli Type). The bundled `.woff2` is licensed for this project; don't redistribute beyond `figures/blog/`. |

The font file lives at `figures/blog/fonts/GT-Pressura-Extended-Regular.woff2`.
The blog HTML figures reference it via relative `url("fonts/...")`. If
you need it elsewhere in the repo, copy to `theme/fonts/`.

## 2. Geist — body, UI, charts

Geometric sans by Vercel. The default body face. Use weights 400 / 500 / 700.

| Property | Value |
|---|---|
| Weights bundled | Regular (400), Medium (500), Bold (700), Variable |
| Source | https://github.com/vercel/geist-font |
| License | OFL (free) — safe to bundle anywhere |

Files: `figures/blog/fonts/Geist-{Regular,Medium,Bold,Variable}.woff2`.
If you only have room for one weight, ship `Geist-Variable.woff2` and set
`font-weight` per element.

## 3. SF Mono / JetBrains Mono / Menlo — kinetic readouts

KD / kon / koff numbers, sequence IDs, axis labels, small uppercase
technical metadata. The stack falls back gracefully across macOS (SF Mono
/ Menlo) and downloaded JetBrains Mono.

```css
font-family: "SF Mono", "JetBrains Mono", Menlo, Consolas, monospace;
```

No bundled file required — these are system fonts on macOS. For
server-side rendering (headless Chrome on Linux), download JetBrains Mono
into `theme/fonts/` and add an `@font-face` declaration.

## Plain HTML loader pattern

For static HTML figures (the seven in `figures/blog/`), each file declares
its own `@font-face` rules at the top, pointing at `fonts/*.woff2`
alongside. Open one in a browser and it renders. No build step.

## Plotting (matplotlib)

`scripts.plotting.apply_theme()` loads `theme/matplotlibrc`, which sets:

| First choice | Fallback 1 | Fallback 2 |
|---|---|---|
| Geist (if installed system-wide) | Roboto | DejaVu Sans |

GT Pressura is **not** used in matplotlib figures — it can't be embedded
in PNG/PDF output without a license check at every render. Geist is the
canonical figure face.
