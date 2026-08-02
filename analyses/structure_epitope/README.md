# Structure & Epitope analysis

Code, results, figures, and the draft + SI for the **Structure & Epitope** contribution.
Owns **main-text Figures 8–12**, **Supplementary Figures S4–S6**, and **Tables S7–S9**.

Tier A = permissive open-source stack only (Biopython, pydssp, scikit-learn, gemmi, scipy, matplotlib).
All seeds fixed; re-running reproduces the committed CSVs + figures modulo float formatting. Inputs come
from the repo root: `data/designs.{csv,parquet}`, `data/structures/<predictor>/`, `data/raw_lab/`, and
citations from `references/*.bib`.

## Files
- `drafts/structure_analysis.md` — section draft. `drafts/supplementary_information.md` — SI (figure/table
  cross-reference + alignment status).
- `build_canonical.py` — stage 0 (canonical design table). `reproduce_all.sh` — one-command Tier A pipeline.
- `jobs/` — 24 analysis + figure/table scripts (pipeline + kinetics, literature, predictor-panel, blender-input prep).
- `figures/` — rendered 2D figures (see map). `figures/blender/` — 3D-figure pipeline (scripts + `FIGURE_GUIDE.md`
  + `sample_renders/`).
- `results/` — 28 results CSVs + `predictor_panel.md` (Table S7) + `metric_predictor_table.md` (Table S8) + `kinetics/`.
- `refs/*.cif` — RCSB crystals 6YYE/6Y6C/6B8O/5ELI (**gitignored**; re-download from RCSB).
- `data/` — `00_inventory.md`, `canonical_schema.json` (provenance/schema; the design data itself is repo-root `data/`).

## Figure / table → file → script
| Item | File | Script |
|---|---|---|
| Fig 8a–d | `figures/blender/sample_renders/fig08_partners.png` | blender Fig-8 pipeline |
| Fig 8e | `figures/fig08e_epitope_landscape.png` | `jobs/figures.py` |
| Fig 9 | `figures/fig09_interface_features.png` | `jobs/figures.py` |
| Fig 10 | `figures/blender/sample_renders/fig10_gallery.png` | blender Fig-10 pipeline |
| Fig 11 | `figures/fig11_metrics_vs_affinity.png` | `jobs/figures.py` |
| Fig 12 | `figures/fig12_kinetic_decomposition.png` | `jobs/kinetics_figure.py` |
| Fig S4 | `figures/figS4_ipsae_reproducibility.png` | `jobs/supp_figures2.py` |
| Fig S5 | `figures/figS5_method_robustness.png` | `jobs/supp_figures.py` |
| Fig S6 | `figures/figS6_literature_kd.png` | `jobs/literature_comparison.py` |
| Table S7 | `results/predictor_panel.md` | `jobs/predictor_panel.py` |
| Table S8 | `results/metric_predictor_table.md` | `jobs/metric_predictor_table.py` |
| Table S9 | `results/metric_predictor_signed_diff.csv` | `jobs/metric_predictor_table.py` |

Figure 8 = crystal-partners panels **a–d** (3D) + per-residue landscape **e** (2D); Figure 10 = the 3D complex
gallery. Filenames, `savefig` calls, and in-figure titles all use the final numbering.

## Reproduce
```bash
bash analyses/structure_epitope/reproduce_all.sh          # Tier A: Figs 8e/9/11, Figs S4/S5, reference epitopes, results CSVs
J=analyses/structure_epitope/jobs
python $J/predictor_panel.py && python $J/metric_predictor_table.py           # Tables S7–S9
python $J/kinetics_prep.py && python $J/kinetics_stats.py && python $J/kinetics_figure.py   # Fig 12
python $J/literature_comparison.py                                            # Fig S6
python $J/export_pdbs_for_blender.py                                          # 3D-figure input PDBs; then figures/blender/FIGURE_GUIDE.md
```
Gitignored intermediates (`designs_canonical.parquet`, `footprints_*.parquet`, blender `pdbs/`, `*.blend`,
`refs/*.cif`) are regenerated / re-downloaded by the steps above.

## Notes
- **Paths:** most `jobs/*.py` hard-code an absolute repo `ROOT`; edit it if the tree isn't at that path.
- **Tier A only:** Tier B (PyRosetta chemistry, TM-align clustering, EvoEF2 alanine scan) and GPU structure
  regeneration are not in this pipeline; `results/alanine_scan.csv` (EvoEF2) ships as committed data.
- **Cohorts:** never collapse human vs. agent — every hit-rate / KD / diversity number splits by cohort.
