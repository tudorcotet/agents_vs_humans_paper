# Data package — every artifact, one page

The full data package for the TREM2 hackathon paper. ~162 MB total
across ~620 files. Binary artifacts are tracked via Git LFS — run
`git lfs install && git lfs pull` after clone.

## Tabular data — `data/`

| path | shape | size | load with |
|---|---|---|---|
| `data/designs.csv`              | 141 × 123    | 196 KB | `load_designs()` |
| `data/designs.parquet`          | 141 × 123    | 128 KB | `load_designs()` (preferred, typed) |
| `data/designs.fasta`            | 141 records  | 24 KB  | `Bio.SeqIO.parse(..., 'fasta')` |
| `data/raw_lab/bli_results.csv`  | 100 rows     | 24 KB  | `pd.read_csv` |
| `data/raw_lab/bli_replicates.csv` | 215 rows   | 16 KB  | `load_replicates()` |
| `data/target/trem2_construct.fasta` | 1 record, 175 aa | 1 KB | `Bio.SeqIO.parse(..., 'fasta')` |
| `data/controls/known_binders.fasta` | 2 records (AL002, VHB937 placeholders) | 1 KB | `Bio.SeqIO.parse(..., 'fasta')` |

The 141 × 123 canonical table includes 42 `pb_*` columns from the
ProteinBase enrichment, populated for the 100 screened designs. See
[`docs/DATA.md`](docs/DATA.md) for every column.

## ProteinBase mirror — `data/proteinbase/`

100 screened designs only. Non-screened designs are absent. Filter on
`submitted_to_lab=True` or `pb_id.notna()` before joining.

| path | count | size | what |
|---|---|---|---|
| `data/proteinbase/boltz2/design_NNN.cif`   | 100 | 14 MB | Boltz-2 complex predictions (re-fold). Chain A = TREM2, chain B = binder. |
| `data/proteinbase/esmfold/design_NNN.cif`  | 100 | 5 MB  | ESMFold single-chain (binder only). |
| `data/proteinbase/pae/design_NNN.json`     | 100 | 88 MB | PAE matrices from Boltz-2, residue × residue, Å. |
| `data/proteinbase/images/design_NNN.png`   |  99 | 27 MB | Stylised cartoon renders of the ESMFold structure. |
| `data/proteinbase/sensorgrams/design_NNN_repNN_{spr,bli}.json` | 215 (193 SPR + 22 BLI) | 27 MB | Raw kinetic-curve traces per replicate. |

The local paths are columns in `designs.csv`:

```python
from scripts.utils import load_designs

df = load_designs(only_screened=True)
df[['design_id', 'pb_boltz2_cif', 'pb_esmfold_cif', 'pb_pae_json', 'pb_stylized_png']]
```

Sensorgrams have multiple files per design — glob the folder:

```python
from pathlib import Path
rep_files = sorted(Path('data/proteinbase/sensorgrams').glob('design_017_rep*.json'))
```

Upstream mirror: `https://proteinbase-pub.t3.storage.dev/`. See
[`data/proteinbase/README.md`](data/proteinbase/README.md) for the JSON
schemas.

## Blog figures — `figures/blog/`

| path | count | size | what |
|---|---|---|---|
| `figures/blog/*.html` | 7 | 904 KB (with fonts + assets) | Hand-authored inline-SVG / HTML/CSS figures from the blog post. Open in a browser, no build step. |

The seven files: `article-graphical-abstract.html`,
`article-fig1-study-design.html`, `article-fig2-headline-distributions.html`,
`article-fig3-tool-monoculture.html`, `article-fig4-winners.html`,
`article-fig-top5-binders.html`, `article-fig-target.html`. See
[`figures/blog/README.md`](figures/blog/README.md).

## Git LFS

All binary artifacts ship via Git LFS. A fresh clone pulls ~162 MB:

```bash
git clone <repo>
cd agents_vs_humans_paper
git lfs install
git lfs pull
```

Without LFS, every CIF / JSON / PNG / parquet / font resolves to a
~130-byte pointer stub and analyses fail to read them. The tracked
patterns are in `.gitattributes`. See
[`data/proteinbase/README.md`](data/proteinbase/README.md#git-lfs).

## What's not in the package

- **AF3 / Chai-1 folds** — only Boltz-2 + ESMFold are mirrored.
- **Per-residue scores** (pLDDT / ipSAE per-residue) — recompute from
  shipped CIFs / PAEs.
- **Raw hackathon submission CSVs** — organiser-private; the pooled
  output is `designs.csv`.

## One-line summary

`data/designs.csv` (141 × 123, 42 `pb_*` cols) joins to 100 Boltz-2
CIFs + 100 ESMFold CIFs + 100 PAEs + 99 PNGs + 215 sensorgrams under
`data/proteinbase/`, plus per-design (100) and per-replicate (215) BLI
tables under `data/raw_lab/`. ~162 MB via Git LFS. Self-contained.
