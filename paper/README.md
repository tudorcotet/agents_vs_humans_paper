# `paper/`

LaTeX skeleton + markdown sections for the manuscript.

## Layout

```
paper/
├── main.tex                 # the LaTeX entrypoint
├── references.bib           # bibliography (curated from references/*.bib)
├── figures/                 # symlinks or copies into figures/paper/
└── sections/
    ├── abstract.md
    ├── introduction.md
    ├── methods.md
    ├── results.md
    └── discussion.md
```

## Workflow

Sections are written in markdown. `main.tex` glues them in via
`\input{...}` after a `pandoc` pass (or by hand-translation when the
content stabilises).

```bash
# (optional) one-shot markdown → tex
pandoc sections/results.md -o sections/results.tex
```

Figures live under `../figures/paper/`; create symlinks under
`./figures/` if you want LaTeX to find them with relative paths:

```bash
ln -s ../figures/paper figures
```

## Voice

Methods-paper voice (measured, clinical). See
[`../docs/STYLE_GUIDE.md`](../docs/STYLE_GUIDE.md). Keep the paper in
that register — blog / social copy lives elsewhere.

## Citations

Curated entries live in `references.bib`. Topic-organised libraries
(`references/trem2.bib`, `methods.bib`, `agentic_biology.bib`,
`competitions.bib`) live at the repo root and feed `references.bib` —
add new entries to the topic file first, then re-curate.
