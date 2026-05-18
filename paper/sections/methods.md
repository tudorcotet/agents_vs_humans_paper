<!--
Methods — ~1500 words. Subsections in this order:

1. Target preparation
2. Hackathon cohort design (human vs agent teams; how teams self-reported method)
3. In-silico selection (top-100 by Boltz-2 ipSAE)
4. BLI assay
5. Replicate aggregation (selected & ~excluded & fixed rule)
6. Statistical analyses (Fisher, Mann-Whitney, Spearman + bootstrap CI, BH-FDR)
7. Data availability

Cite the column name in `data/designs.csv` next to every quantitative claim,
so a reader running `make analysis-all` can re-derive it.
-->

## Target preparation

TBD. (Acro Biosystems `TR2-H52H5`: initiator Met + residues 19–174 of
UniProt Q9NZC2 + GGGSGGGS linker + 10×His. 175 aa, 19.3 kDa. Verbatim
sequence at `data/target/trem2_construct.fasta`.)

## Hackathon cohort design

TBD. (10 human teams = 81 designs; 6 agent teams = 60 designs; 141 total.
The agent teams were end-to-end autonomous on the bioArena platform with
tool access — they are not raw language models emitting sequences.)

## In-silico selection

TBD. (Top 100 by Boltz-2 ipSAE went to the BLI assay; the remaining 41
designs are not wet-lab characterised.)

## BLI assay

TBD. (Adaptyv standard BLI protocol, target-loaded biosensors,
concentration series 0–1000 nM.)

## Replicate aggregation

Per-design KD is the arithmetic mean of `kd_nM` over replicates passing
the filter `selected=true AND excluded=false AND fixed=true` — the
curator-confirmed set. Replicates flagged `staircase=true`,
`unexpected_order=true`, or `confidence=low` stay in the per-replicate
file (`data/raw_lab/bli_replicates.csv`) for transparency and enter the
mean unless the curator excluded them.

## Statistical analyses

- Hit-rate comparisons: Fisher's exact (two-sided), with the raw $2\times2$ counts.
- Distribution comparisons (`pkd`, `sequence_length`): Mann-Whitney U
  (two-sided), with the median difference reported.
- Correlations (`submitted_ipsae` vs `pkd`): Spearman $\rho$ with
  bootstrap 95% CI (1000 iterations, percentile method).
- Multiple comparisons: Benjamini-Hochberg FDR at q=0.10 when more than
  five tests are run within a single analysis.

All statistical code is in `analyses/human_vs_agent/headline_stats.py`.

## Data availability

The canonical 141-row design table is `data/designs.csv`. Per-replicate
BLI rows are in `data/raw_lab/bli_replicates.csv`. The target sequence is
in `data/target/trem2_construct.fasta`. The full repo (code + data +
figures) is at https://github.com/adaptyvbio/muni_trem2_paper (TBD).
