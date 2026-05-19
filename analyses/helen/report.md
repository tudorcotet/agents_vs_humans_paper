# Manuscript figures (analyses/helen)

**Headline:** the Notion-spec figures, rendered as single SVGs in blog styling (Agent `#30C5F5`, Human `#142933`).

## Method

`load_designs()` → matplotlib (`apply_theme()` + blog tokens), one SVG per figure. ESM-2 (`esm2_t33_650M_UR50D`), Foldseek TM (local static binary, ESMFold CIFs), and epitope footprints (gemmi, Boltz-2 CIFs, 5 Å contacts) are precomputed into `*_cache.npz` by the sibling scripts; figure code is deterministic and reads the caches.

## Results

- `figures/paper/fig1_ipsae_distribution.svg`
- `figures/paper/fig2_expression_binding.svg`
- `figures/paper/fig4_design_methods.svg`
- `figures/paper/fig6_sequence_length.svg`
- `figures/paper/fig7_metric_roc.svg`
- `figures/paper/fig_identity_heatmap.svg`
- `figures/paper/fig_foldseek_clustering.svg`
- `figures/paper/fig_epitope_regions.svg`
- `figures/paper/fig5_esm2_umap.svg`

## Caveats

- **Fonts** are embedded per-SVG as base64 @font-face (blog Geist + GT Pressura Extended woff2); self-contained, overrides the matplotlibrc Helvetica/DejaVu policy by design. GT Pressura is a paid face — embedding in distributed artifacts is a licensing call.
- **Fig 3** (BLI sensorgrams): excluded by request.
- **ESM-2 pLL** is a length-normalized single-pass log-likelihood; true masked-marginal is infeasible on CPU/MPS at 650M. Document if revised.
- **Fig 4 taxonomy** (design_method_normalized → 5-category, Cotet 2025) is an interpretation; PXDesign and BindCraft bucket as De novo (BindCraft per the Notion draft's own grouping).
- **Foldseek/epitope** are screened-only (100/141); the 41 non-screened designs have no ProteinBase CIF.
- Per-cohort epitope/structure counts inherit the top-100 selection imbalance (65 human / 35 agent).
- ESM-2 figures present: True (cache found).
