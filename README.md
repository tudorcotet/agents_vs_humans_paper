# Agents vs humans — TREM2 paper

On Feb 28 2026, sixteen teams spent a single day in San Francisco designing
binders against **TREM2**. Ten of those teams were humans. Six were
autonomous AI agents running on bioArena with no human in the loop.

141 designs landed. Top 100 by Boltz-2 ipSAE went to the wet lab (BLI).
The result: **37 binders, 11 non-expressing, 52 non-binders**. This repo
holds the data, the code, and the manuscript source for the paper asking
the only question worth asking: **how did the agents do?**

## The cohort

|              | teams | designs | screened | binders |
|---           |---    |---      |---       |---      |
| **Human**    | 10    | 81      | 65       | 25      |
| **Agent**    |  6    | 60      | 35       | 12      |
| **Total**    | 16    | 141     | 100      | 37      |

Human teams: NovoFy, MRAZS, EuroBros, DeNovo, StanFold, 1000Tokens, crow,
BART.bio, BraiNSEY, Marcel.
Agent teams: claude-sonnet-4.6, qwen-3.5-plus, grok-4.1-fast,
Gemini 3.1 Pro, GLM 5, GPT 5.2.

`is_human` in `data/designs.csv` is the truth. Don't infer cohort from a
team name — Marcel is a human, claude-sonnet-4.6 isn't.

## Cloning

The big binary data (structures, sensorgrams, PAEs, images, parquet,
fonts, raw-lab CSVs) lives behind Git LFS. Pull it after cloning:

```bash
git clone <repo>
cd agents_vs_humans_paper
git lfs install
git lfs pull        # ~162 MB
```

Without `git lfs pull`, every binary path resolves to a ~130-byte
pointer stub and the analyses will fail to read them. See
[`data/proteinbase/README.md`](data/proteinbase/README.md#git-lfs)
for the full list of tracked patterns.

## Quick start

```bash
make setup        # uv sync
make analysis-all # re-runs the four canonical analyses against data/designs.csv
```

No GPU, no credentials, no internal tooling needed. `data/designs.csv`
is self-contained.

```python
from scripts.utils import load_designs
df = load_designs()                       # 141 rows
df = load_designs(only_screened=True)     # 100 rows
df = load_designs(only_hits=True)         # 37 rows
```

## Layout

```
agents_vs_humans_paper/
├── data/                    ⭐ one CSV, every annotation
│   ├── designs.csv          one row per design, 141 × 123 (42 pb_* cols, 100/141)
│   ├── designs.parquet      same, typed
│   ├── designs.fasta        sequences only
│   ├── target/              TREM2 ectodomain construct (Acro TR2-H52H5, 175 aa)
│   ├── controls/            AL002 / VHB937 reference binders
│   ├── raw_lab/             per-design (100) and per-replicate (215) BLI tables
│   └── proteinbase/         ProteinBase mirror: 100 CIFs + 100 CIFs + 100 PAEs
│                            + 99 PNGs + 215 sensorgrams (~161 MB, screened only)
├── docs/                    everything you need to collaborate
│   ├── DATA.md              every column of designs.csv
│   ├── ANALYSES.md          how to plug yours in
│   ├── FIGURES.md           figure style + the rendering pipeline
│   ├── STYLE_GUIDE.md       voice + palette
│   ├── GLOSSARY.md          ipSAE, pKD, modality, …
│   └── REPRODUCIBILITY.md   how the CSV was assembled
├── analyses/                one subdir per analysis (yours goes here)
│   ├── _template/           copy this to start
│   ├── human_vs_agent/      cohort stats (Fisher, MW, Spearman, χ²)
│   ├── leaderboard/         top-N tables
│   ├── sequence_diversity/  identity, k-mer entropy
│   └── methods/             method × outcome cross-tabs
├── scripts/
│   ├── utils/load_data.py   the canonical loader — import this
│   ├── utils/stats.py       Fisher, MW, bootstrap CI, BH-FDR helpers
│   └── plotting/_common.py  apply_theme + save_figure
├── theme/                   palettes.json, matplotlibrc, fonts/
├── figures/
│   ├── paper/               rendered figures for the manuscript
│   ├── blog/                the 7 hand-authored HTML figures from the blog post
│   └── exploration/         gitignored scratch
├── paper/                   LaTeX main.tex + markdown sections
└── references/              BibTeX by topic
```

## Adding your analysis

```bash
cp -r analyses/_template analyses/yourname
$EDITOR analyses/yourname/main.py
# add an [tasks."analysis:yourname"] block in mise.toml
mise run analysis:yourname
```

Full conventions in [`CONTRIBUTING.md`](CONTRIBUTING.md). Pattern walkthrough
in [`docs/ANALYSES.md`](docs/ANALYSES.md). Style guide in
[`docs/STYLE_GUIDE.md`](docs/STYLE_GUIDE.md).

## Known caveats

- **Top-100 cherry-pick.** 41 designs never reached the wet lab. Hit
  rates are conditional on surviving the in-silico ipSAE triage.
- **Selection imbalance.** Humans submitted 81 designs and got 65 into
  the top 100. Agents submitted 60 and got 35. Cohort hit-rate
  comparisons inherit that imbalance.
- **Single target.** TREM2 IgSF only. Generalisation to other targets is
  an assumption.
- **Per-team N is thin.** Marcel = 1 design. DeNovo = 2. Per-team
  rankings are anecdotes; per-cohort numbers are the load-bearing claim.
- **Method labels are self-reported.** `Foundry` means "the agent picked
  its own tools." It is not a single algorithm. Bucket as `Other`.
- **One hit has no KD fit.** Design 5 (NovoFy/RFDiffusion) is labelled
  `weak` by the curator but no replicate yielded a fittable curve.
  Filter on `kd_arith_mean_nM_all.notna()` for KD analyses; on
  `is_hit` for binary hit analyses.
- **Folding metrics are 136/141 (legacy) or 100/141 (re-fold).** 5
  bottom-ranked non-screened designs have null `boltz2_*` and
  `submitted_ipsae` — Boltz-2 didn't produce a useful fold for them. The
  ProteinBase re-fold (`pb_boltz2_*`) is 100/141, screened only. The two
  match closely but not exactly (max |Δ ipSAE| = 0.18). Pick one,
  document which.
- **ProteinBase enrichment is screened-only.** The 42 `pb_*` columns and
  every artifact under `data/proteinbase/` cover the 100 designs that
  made the wet lab. The 41 non-screened designs are null across `pb_*`.
  Filter on `pb_id.notna()` or `submitted_to_lab=True`.

## Citation

```bibtex
@dataset{agents_vs_humans_trem2_2026,
  author    = {{Adaptyv Bio} and {bioArena}},
  title     = {Agents vs humans on TREM2: a one-day binder design hackathon},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {TBD},
}
```

Machine-readable: [`CITATION.cff`](CITATION.cff).

## License

Code is Apache-2.0. Data and figures are CC-BY-4.0. Manuscript text is © the
authors until acceptance, then CC-BY-4.0. See [`LICENSE`](LICENSE).
