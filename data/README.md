# `data/`

The canonical inputs every analysis reads. **Don't edit these files by
hand.** If something is wrong, raise it in a PR description.

## Layout

```
data/
├── designs.csv          ⭐ ONE row per design, 141 rows × 356 cols, ALL annotations
│                          (curated metadata + pb_* + tp_* + b2_*/px_*/chai_*/af2m_*
│                          + esm_pll_*/netsolp_*/saprot_* + prodigy_*/destress_* per model)
├── designs.parquet      same, typed (preferred for code)
├── designs.fasta        sequences only; headers: >{id}|{name}|{team}|{method}
├── target/
│   └── trem2_construct.fasta     Acro TR2-H52H5, 175 aa, 19.3 kDa
├── controls/
│   ├── known_binders.fasta       AL002 / VHB937 (placeholders until patent extract)
│   └── README.md
├── raw_lab/
│   ├── bli_results.csv           one row per design (after replicate filter), 100 rows
│   └── bli_replicates.csv        one row per replicate (long form), 215 rows
├── structures/                   one subfolder per model
│   ├── boltz2/         100 × Boltz-2 complex CIFs   (14 MB)   ← ProteinBase mirror
│   ├── esmfold/        100 × ESMFold binder CIFs    (5 MB)    ← ProteinBase mirror
│   └── proteintyper/  ≤41 × ESMFold binder CIFs               ← local rerun (the gap)
├── metrics/                      per-model raw scalars + matrices
│   ├── pae/            100 × Boltz-2 PAE matrices (JSON)    (88 MB)
│   └── proteintyper/  ≤41 × full TyperJobOutput (JSON)
├── images/             99 + ≤41 × stylised PNG renders        (27 MB+)
└── sensorgrams/        215 × kinetic curves (JSON, 193 SPR + 22 BLI, 27 MB)
```

The `pb_*` columns in `designs.csv` (42 of them, 100/141 coverage) carry
ProteinBase's enriched metrics and the relative paths into the
`structures/`, `metrics/`, `images/`, `sensorgrams/` trees. See
[`../docs/DATA.md` §11](../docs/DATA.md#11-proteinbase-enriched-data-100141--screened-only).
The 41 non-screened designs are filled in by
[`scripts/folding/run_proteintyper.py`](../scripts/folding/run_proteintyper.py)
— Modal ProteinTyper with the default `full_monomer` recipe (target-only
MSA, sequence-only input), output into `structures/proteintyper/` +
`metrics/proteintyper/` + `images/`.

### Design rationale

- **`structures/` is keyed by model, not by source.** Consistent
  subfolder-per-model layout so analysis code is portable across paper
  repos.
- **Sensorgrams and renders live *outside* `structures/`.** A sensorgram
  is a kinetic trace, a render is a picture; neither is a model output.
- **`metrics/` holds per-model raw payloads** (PAE matrix, full Typer
  JSON) while scalar summaries land in `designs.csv` `pb_*` columns.

## How to load

Always go through the canonical loaders:

```python
from scripts.utils import load_designs, load_replicates

df = load_designs()                       # 141 rows
df = load_designs(only_screened=True)     # 100 rows
df = load_designs(only_hits=True)         # 37 rows (note: design 5 has no KD fit)
reps = load_replicates()                  # per-replicate
```

See [`../docs/DATA.md`](../docs/DATA.md) for every column.

## The target

`data/target/trem2_construct.fasta` is the exact construct used in the
BLI assay: Acro Biosystems catalogue **TR2-H52H5**, 175 aa, 19.3 kDa.
That's the initiator methionine + TREM2 ectodomain (residues 19–174 of
UniProt Q9NZC2) + GGGSGGGS linker + 10×His. Designers fold against this
sequence, not the bare IgSF domain.

## Provenance

| subtree | source | how regenerated |
|---|---|---|
| `structures/{boltz2,esmfold}/` (100 designs) | ProteinBase public mirror, `https://proteinbase-pub.t3.storage.dev/` | `uv run python scripts/data/mirror_proteinbase.py` (to add — currently hand-fetched) |
| `metrics/pae/` (100 designs)                 | ProteinBase public mirror | same |
| `images/` (99 designs from mirror)           | ProteinBase public mirror | same |
| `sensorgrams/` (215 traces)                  | ProteinBase public mirror | same |
| `structures/proteintyper/` (≤41 designs)     | Modal ProteinTyper rerun (target-only MSA) | `mise run rerun:typer` |
| `metrics/proteintyper/` (≤41 designs)        | Modal ProteinTyper rerun | same |
| `images/` (extra ≤41 from rerun)             | Modal ProteinTyper rerun | same |

## What's NOT here

- **Boltz-2 / Protenix complex CIFs for the 41 non-screened designs** —
  they shipped without complex predictions. ProteinTyper only gives us
  binder-monomer ESMFold. Re-fold with Boltz-2 or Protenix v2 to fill
  the gap; outputs would land in `structures/boltz2/` /
  `structures/protenix/`.
- **AF3 / Chai-1 folds** — regenerable from `designs.fasta`.
- **Per-residue pLDDT / ipSAE** — only complex-level scalars are in the
  CSV. Recompute from the shipped CIFs / PAEs if you need per-residue.
- **Raw hackathon submission CSVs** — held back by the organisers; the
  pooled result is `designs.csv`.

## Regeneration

`data/designs.csv` is built upstream by the muni / Adaptyv joint
pipeline. External collaborators work off the shipped CSV; see
[`../docs/REPRODUCIBILITY.md`](../docs/REPRODUCIBILITY.md) for the full
provenance. The complex-prediction + monomer-typer rerun pipeline is
shipped in this repo and runs on Modal:

```bash
# ProteinTyper — fills the 41 non-screened ESMFold gaps. Reads three env
# vars (PROTEINTYPER_SUBMIT_URL, PROTEINTYPER_RETRIEVE_URL,
# PROTEINTYPER_API_TOKEN). The script errors with a clear message if any
# are unset. Point them at your own ProteinTyper deployment.
mise run rerun:typer

# Complex prediction vs TREM2 — three Modal apps, --detach by default.
# Defaults to protenix-v2 (pre-fetched from HuggingFace
# `TMF001/pxdesign-weights`; bypasses the ByteDance gate on volces.com).
# Modal workspace + volume names are ENV-overridable
# (MODAL_APP_NAME, MODAL_RESULTS_VOLUME, MODAL_PROTENIX_VOLUME).
mise run rerun:complex                     # boltz2 + protenix + chai1
mise run rerun:complex -- --download-only  # pull JSON + CIF off the volume
mise run rerun:protenix                    # just the Protenix branch

# Pool every per-model JSON into one wide CSV.
mise run build:grand                       # → data/grand_metrics.csv

# Override the Protenix weight (default is protenix-v2):
PROTENIX_MODEL_NAME=protenix_base_default_v0.5.0 mise run rerun:protenix
```

`data/grand_metrics.csv` is the long-format companion to `designs.csv`:
one row per design × every metric (`b2_*`, `px_*`, `chai_*`, `tp_*`)
plus per-model `cif_exists` flags and two cross-folder consensus
columns (`ipsae_pass_3folders`, `iptm_pass_3folders`). 141 × ~90.

## Git LFS

Every binary artifact in `structures/`, `metrics/`, `images/`,
`sensorgrams/`, the top-level `designs.parquet`, and `raw_lab/*.csv`
lives behind [Git LFS](https://git-lfs.com). The repo's `.gitattributes`
is the source of truth for the tracked patterns.

A fresh clone needs LFS installed and the blobs pulled:

```bash
git clone <repo>
cd agents_vs_humans_paper
git lfs install
git lfs pull        # ~162 MB across ~620 files
```

Without LFS, every CIF / JSON / PNG / parquet / font is a **pointer stub
of <200 bytes**; any code that tries to parse one fails with
"unexpected end of file". Run `git lfs pull` and retry.
