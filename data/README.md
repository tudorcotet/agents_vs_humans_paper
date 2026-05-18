# `data/`

The canonical inputs every analysis reads. **Don't edit these files by
hand.** If something is wrong, raise it in a PR description.

## Layout

```
data/
├── designs.csv          ⭐ ONE row per design, 141 rows × 123 cols, all annotations
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
└── proteinbase/                  ProteinBase public mirror, screened designs only
    ├── boltz2/      100 × Boltz-2 complex CIFs   (14 MB)
    ├── esmfold/     100 × ESMFold binder CIFs    (5 MB)
    ├── pae/         100 × PAE matrices (JSON)    (88 MB)
    ├── images/       99 × stylised PNG renders   (27 MB)
    ├── sensorgrams/ 215 × kinetic curves (JSON, 193 SPR + 22 BLI, 27 MB)
    └── README.md
```

The `pb_*` columns in `designs.csv` (42 of them, 100/141 coverage) carry
ProteinBase's enriched metrics and the relative paths into the
`proteinbase/` tree. See [`../docs/DATA.md` §11](../docs/DATA.md#11-proteinbase-enriched-data-100141--screened-only)
and [`proteinbase/README.md`](proteinbase/README.md).

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

## What's NOT here

- **AF3 / Chai-1 folds** — only Boltz-2 and ESMFold are mirrored in
  `proteinbase/`. AF3 / Chai are regenerable from `designs.fasta`.
- **Per-residue pLDDT / ipSAE** — only complex-level scalars are in the
  CSV. Recompute from the shipped CIFs / PAEs if you need per-residue.
- **Raw hackathon submission CSVs** — held back by the organisers; the
  pooled result is `designs.csv`.

## Regeneration

`data/designs.csv` is built upstream by the bioArena / Adaptyv joint
pipeline. The pieces (sequence pooling, Boltz-2 folding, ipSAE scoring,
BLI fetch, replicate aggregation) are outlined in
[`../docs/REPRODUCIBILITY.md`](../docs/REPRODUCIBILITY.md). External
collaborators don't run the upstream pipeline — they work off the
shipped CSV.
