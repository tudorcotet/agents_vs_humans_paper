# `analyses/helen/` — manuscript figures

## What this answers

Renders the figures specified in the Notion manuscript draft as single
self-contained SVGs, styled to match the hand-authored blog figures.

## Run

```bash
uv run python analyses/helen/esm2_compute.py     # ESM-2 650M cache (one-time)
uv run python analyses/helen/foldseek_cluster.py # Foldseek TM cache (one-time)
uv run python analyses/helen/epitope_cluster.py  # epitope cache (one-time)
mise run analysis:helen                          # render every SVG
```

The three cache scripts are heavy/one-time; `main.py` is deterministic and
reads their `*_cache.npz`. Missing caches are skipped, not errors.

## Inputs

- `data/designs.csv` via `scripts.utils.load_designs()`.
- `data/proteinbase/{esmfold,boltz2}/*.cif` (LFS) for structure/epitope.

## Outputs (`figures/paper/`)

- `fig1_ipsae_distribution.svg` — ipSAE by cohort, top-100 cutoff line.
- `fig2_expression_binding.svg` — expression + hit rate, Wilson 95% CI.
- `fig4_design_methods.svg` — 5-category design taxonomy by cohort.
- `fig5_esm2_umap.svg` — UMAP of ESM-2 650M embeddings (cohort + team).
- `fig6_sequence_length.svg` — sequence-length distributions.
- `fig7_metric_roc.svg` — metric ROC incl. ESM-2 pLL.
- `fig_identity_heatmap.svg` — pairwise identity, within + between team.
- `fig_foldseek_clustering.svg` — Foldseek TM-score structural clustering.
- `fig_epitope_regions.svg` — TREM2 epitope usage + distinct-patch count.
- `analyses/helen/{report.md, summary.json}`.

Fig 3 (BLI sensorgrams) excluded by request.

## Caveats

- **Fonts embedded:** each SVG inlines the blog woff2 (`Geist` 400/500/700,
  `GT Pressura Extended` 400) as base64 `@font-face`, so figures are
  self-contained and render in the blog typeface anywhere — overriding the
  `theme/matplotlibrc` Helvetica/DejaVu policy on purpose (blog-match).
  GT Pressura Extended is a paid display face; embedding it in a
  distributed artifact is a licensing consideration, accepted here.
- **Palette deviation:** cohort hues are the blog tokens `core.cyan`
  (#30C5F5) / `core.ink_navy` (#142933), per the user's blog-match
  requirement — not `palettes.agent_vs_human.Human` (#1FE48F, green).
  Both are `theme/palettes.json` brand tokens; deliberate, documented.
- **ESM-2 pLL** is a length-normalized single-pass log-likelihood
  (`esm2_t33_650M_UR50D`); true masked-marginal is infeasible on CPU/MPS.
- **Fig 4 taxonomy** crosswalk is an interpretation (Cotet 2025);
  PXDesign and BindCraft bucket as De novo.
- Foldseek/epitope are screened-only (100/141) and inherit the top-100
  selection imbalance (65 human / 35 agent).
