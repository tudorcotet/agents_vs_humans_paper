# Reproducibility

How `data/designs.csv` was assembled, and what it takes to rebuild it.

## What's reproducible from this repo alone

Everything under `analyses/` re-runs from the checked-in `data/designs.csv`.
You need `uv sync` (or `make setup`). No credentials, no GPU. **This is
the only reproducibility surface most collaborators ever need.**

```bash
make setup
make analysis-all
# every report.md and summary.json under analyses/ is regenerated
```

## What needs upstream tooling

Regenerating `data/designs.csv` itself involves:

1. **Pool raw submissions.** 16 team CSVs → 141 rows, written by the
   muni hackathon organisers. Raw submission CSVs aren't part of
   the public release.
2. **Folding pipeline.** 141 sequences × 4 models (Boltz-2, AF3, Chai-1
   proxy, ESMFold). ~9 GPU-hours per pass. Outputs ~3 GB of PDB files.
3. **Metric extraction.** ipSAE / pTM / ipTM / pLDDT from the folded
   structures. CPU-only, <2 min.
4. **BLI data fetch.** Per-replicate measurements from the lab LIMS.
5. **Replicate aggregation.** Applies the canonical
   `selected=true AND excluded=false AND fixed=true` filter, computes
   `kd_arith_mean_nM_all` and `pkd_arith_mean`.
6. **Homology evidence.** Motif classifier + mmseqs2 vs SwissProt /
   SAbDab / UniRef50.
7. **Modality classification.** Motif + homology + structure call
   collapses raw sequences into `Peptide` / `Miniprotein` / `Large
   miniprotein` / `Nanobody` / `scFv`.
8. **Build canonical.** Joins all of the above and writes
   `data/designs.csv` and `data/designs.parquet`.

The pieces 2–8 are the upstream pipeline. The shipped CSV is the output;
that's what collaborators read.

## Pinning

- Python 3.11 (in `pyproject.toml` and `mise.toml`).
- Deps locked in `uv.lock` (created by `uv sync`; commit it).
- Seeds: bootstrap CIs use seed=0 (see
  `analyses/human_vs_agent/headline_stats.py`). Diversity sampling uses
  `random.Random(0)`.

## Validating a fresh checkout

```bash
# 1. Install
make setup

# 2. Smoke-test the loader
python -c "from scripts.utils import load_designs; df = load_designs(); print(df.shape, int(df.is_human.sum()), 'humans')"
# Expect: (141, 81) 81 humans

# 3. Run all analyses
make analysis-all

# 4. Diff against the committed reports
git diff --stat analyses/
# Expect: nothing — the reports are deterministic.
```

If diffs appear, the likely culprits:

- Different python/numpy/scipy minor version → tiny floating-point shifts.
- The canonical CSV was regenerated upstream → commit the updated reports
  and `data/designs.csv` together.
- A bug you introduced → revert and try again.

## Data lineage

```
upstream pipeline (held back)              ProteinBase public mirror      ProteinTyper rerun (local)
        │                                          │                              │
        ▼                                          ▼                              ▼
data/designs.csv + designs.parquet      data/structures/{boltz2,esmfold}/   data/structures/proteintyper/
   ⭐ source of truth (shipped)          data/metrics/pae/                   data/metrics/proteintyper/
   141 × 123, 42 pb_* cols              data/images/  data/sensorgrams/     data/images/  (≤41 designs)
        │                                100 CIFs + 100 CIFs + 100 PAEs +    ESMFold CIFs + Typer JSON for
        │                                 99 PNGs + 215 sensorgrams           the 41 non-screened designs
        └────────────────────────────────────┬─────────────────────────────────────┘
                                             │ joined via pb_* path columns
                                             ▼
                                scripts.utils.load_designs()
                                             │
                                             ▼
                                analyses/*/main.py  →  report.md + summary.json
                                             │
                                             ▼
                                figures/paper/*.{png,pdf,svg}   +   figures/blog/*.html (hand-authored)
                                             │
                                             ▼
                                paper/main.tex
```

The three top-of-pipe sources are independent: `data/designs.csv` is
built by the muni / Adaptyv joint pipeline; the ProteinBase mirror
under `data/structures/{boltz2,esmfold}/`, `data/metrics/pae/`,
`data/images/`, and `data/sensorgrams/` is pulled from
`https://proteinbase-pub.t3.storage.dev/` (the same blobs that back
[proteinbase.bio](https://proteinbase.bio)); the ProteinTyper rerun fills
in the 41 non-screened designs by calling a ProteinTyper-compatible HTTP
endpoint (URLs + token come from env vars, see
`scripts/folding/run_proteintyper.py`). The CSV references all three
via the `pb_*` path columns — keep them in sync if you regenerate one.

## What's shipped vs held back

| artifact | shipped? |
|---|---|
| `data/designs.csv` / `.parquet`           | ✅ 141 × 123, including 42 `pb_*` columns |
| `data/designs.fasta`                       | ✅ |
| `data/raw_lab/bli_*.csv`                   | ✅ per-design (100) + per-replicate (215) |
| `data/target/*.fasta`                      | ✅ Acro TR2-H52H5, 175 aa |
| `data/controls/*.fasta`                    | ✅ AL002 / VHB937 placeholders (patent extract pending) |
| `data/structures/boltz2/*.cif`             | ✅ 100 Boltz-2 complex CIFs (14 MB) — re-fold |
| `data/structures/esmfold/*.cif`            | ✅ 100 ESMFold binder CIFs (5 MB) |
| `data/structures/proteintyper/*.cif`       | ✅ ≤41 ESMFold binder CIFs from the rerun |
| `data/metrics/pae/*.json`                  | ✅ 100 PAE matrices (88 MB) |
| `data/metrics/proteintyper/*.json`         | ✅ ≤41 `TyperJobOutput` JSONs from the rerun |
| `data/images/*.png`                        | ✅ 99 + ≤41 stylised renders |
| `data/sensorgrams/*.json`                  | ✅ 215 kinetic-curve traces (193 SPR + 22 BLI, 27 MB) |
| `figures/blog/*.html`                      | ✅ 7 hand-authored blog figures |
| AF3 / Chai-1 folds                         | ❌ only Boltz-2 + ESMFold are mirrored |
| Per-residue scores (pLDDT/ipSAE per-res)   | ❌ recompute from shipped CIFs / PAEs |
| Raw hackathon submission CSVs              | ❌ organiser-private; pooled output is `designs.csv` |

All binary artifacts (CIFs, PAEs, PNGs, sensorgrams, parquet, fonts,
raw-lab CSVs) live behind Git LFS — total ~162 MB. `git lfs pull` after
clone.
