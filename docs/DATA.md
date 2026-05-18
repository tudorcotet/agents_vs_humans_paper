# `data/designs.csv` — the canonical table

One row per submitted design. **141 rows, 123 columns.** Every analysis
in the repo reads from here. The parquet sibling has typed columns and
is the right thing to load in code; the CSV is the readable copy for
spreadsheet eyes.

42 of those 123 columns are `pb_*` — the ProteinBase enriched data,
populated for the 100 screened designs only. See [§11](#11-proteinbase-enriched-data-100141--screened-only).

## Load it

```python
from scripts.utils import load_designs, load_replicates

df = load_designs()                       # 141 rows, every column
df = load_designs(only_screened=True)     # 100 rows — sent to BLI
df = load_designs(only_hits=True)         # 37 rows — binder + weak
df = load_designs(only_human=True)        # 81 rows — human cohort
df = load_designs(only_human=False)       # 60 rows — agent cohort

reps = load_replicates()                  # per-replicate (long) BLI rows
reps = reps[reps.selected & ~reps.excluded & reps.fixed]   # canonical filter
```

`load_designs()` coerces the CSV's `True`/`False` strings to pandas
nullable `boolean` — you can write `df.is_hit & df.expressed` without
fighting object dtype.

## Keys

- **`design_id`** (int, 1..141) — universal join key. Stable across the
  pipeline. **Use this everywhere.**
- **`name`** (string) — display name, unique, leading number matches
  `design_id`.

## Column groups

### 1. Identity (141/141)

| column | type | notes |
|---|---|---|
| `design_id`         | int     | 1..141. Universal join key. |
| `name`              | string  | Display name, e.g. `13_MRAZS_mosaic`. Unique. |
| `team`              | string  | One of 16 teams. |
| `is_human`          | bool    | `True` for the 10 human teams (81 designs), `False` for the 6 agents (60 designs). |
| `cohort`            | string  | `human` / `agent`, derived from `is_human` for plotting convenience. |

### 2. Submission metadata (141/141)

| column | type | notes |
|---|---|---|
| `design_method`              | string | Self-reported method, verbatim. |
| `design_method_normalized`   | string | Canonicalised label. See [`GLOSSARY.md`](GLOSSARY.md). |
| `method_family`              | string | `PXDesign` (44), `RFDiffusion` (30), `Other` (25), `BoltzGen` (17), `Hallucination` (8), `Mosaic` (6), `PPIFLOW` (6), `RFPeptides` (3), `BindCraft` (2). |
| `sequence`                   | string | Amino acid sequence, uppercase, whitespace stripped. |
| `sequence_length`            | int    | `len(sequence)`. Range 12..247. |

### 3. Modality (141/141)

Classifier ported from the EGFR repo, stripped of EGFR-specific tiers.

| column | type | notes |
|---|---|---|
| `modality`               | string | `Miniprotein` (119), `Nanobody` (9), `Large miniprotein` (8), `Peptide` (3), `scFv` (2). |
| `modality_family`        | string | `Peptide` / `Miniprotein` / `Antibody`. |
| `modality_motif_call`    | string | Per-design motif call before homology overrides. |
| `modality_blast_best_ref`| string | Closest reference (used for `is_literature_copy`). |
| `modality_blast_best_frac`| float | Fractional identity to the closest reference. |
| `modality_blast_overridden`| bool | True if the motif call was overridden by homology evidence. |
| `is_literature_copy`     | bool | ≥95% identical to a published binder. |
| `is_framework_borrowed`  | bool | Antibody framework borrowed (≥80% V-gene identity). |
| `design_class`           | string | `miniprotein`, `nanobody`, `scfv`, `peptide`, or `not_assigned`. |

### 4. Sequence-derived features (141/141)

| column | type | notes |
|---|---|---|
| `cys_count`         | int   | Number of cysteines. |
| `alanine_fraction`  | float | Alanine residues / length. |
| `has_g4s_linker`    | bool  | Contains GGGGS or SGGGGS. |
| `methionine_leader` | bool  | Starts with `M`. |
| `molecular_weight`  | float | Daltons. |
| `pI`                | float | Isoelectric point. |
| `pLDDT_esmfold`     | float | Mean pLDDT from a single-chain ESMFold pass. |

### 5. In-silico folding metrics (136/141)

Boltz-2 metrics pulled from the bioArena selection table. 5 designs are
null — all bottom-ranked non-screened (ranks 137–141, `selected=no`) for
which Boltz-2 didn't produce a useful fold. Null-check before use.

| column | type | notes |
|---|---|---|
| `submitted_ipsae`       | float | Boltz-2 ipSAE used to rank the cherry-pick. **Primary in-silico metric.** |
| `boltz2_iptm`           | float | Boltz-2 interface pTM. |
| `boltz2_plddt`          | float | Mean per-chain pLDDT. |
| `boltz2_complex_plddt`  | float | Complex-level pLDDT. |
| `boltz2_pdockq`         | float | Predicted DockQ from Boltz-2. |
| `boltz2_min_ipsae`      | float | Minimum ipSAE across re-scored ranks. |

### 6. Selection (141/141)

| column | type | notes |
|---|---|---|
| `submitted_to_lab` | bool  | `True` for the 100 designs in the top-100 ipSAE cherry-pick. |
| `selection_rank`   | int?  | 1..100 if `submitted_to_lab`, else null. |

### 7. Wet-lab outcome (100/141 — only for screened designs)

| column | type | notes |
|---|---|---|
| `binding_label`       | string  | Hierarchy: `binder` > `strong` > `medium` > `weak` > `potential_binder` > `non_binder` > `no_expression` > `unknown`. In this dataset only `binder` (36), `weak` (1), `non_binder` (52), `no_expression` (11) actually appear. |
| `binding_strength`    | string  | Bucketised strength: `none`, `weak`, `medium`, `strong`. |
| `expressed`           | bool    | `True` if the construct expressed. 89 / 11 / 41 (T/F/null). |
| `is_hit`              | bool    | `binding_label ∈ {binder, strong, medium, weak}`. 37 hits. |

### 8. Replicate-aggregated KD (36/141 — only for fittable binders)

Aggregated from `data/raw_lab/bli_replicates.csv` after applying the
canonical filter (`selected=true AND excluded=false AND fixed=true`).

> **Edge case:** 1 hit (design 5, NovoFy/RFDiffusion, `binding_label=weak`)
> has no fittable replicate, so 36 KDs cover the 37 hits. Filter on
> `kd_arith_mean_nM_all.notna()` for KD plots; on `is_hit` for binary
> hit-rate analyses.

| column | type | notes |
|---|---|---|
| `n_replicates_total`             | int   | Total replicates (any flag state). |
| `n_replicates_pushed`            | int   | Replicates with `fixed=true`. |
| `n_spr_replicates`               | int   | Replicates on Carterra SPR. |
| `n_bli_replicates`               | int   | Replicates on Gator BLI. |
| `n_with_kd`                      | int   | Replicates that produced a fit KD. |
| `kd_arith_mean_nM_all`           | float | Arithmetic mean of `kd_nM` over the canonical filter. **Primary KD.** |
| `kd_arith_mean_nM_highconf`      | float | Same, `confidence=high` replicates only. |
| `kd_arith_mean_nM_spr_only`      | float | SPR replicates only. |
| `kd_arith_mean_nM_bli_only`      | float | BLI replicates only. |
| `kd_min_replicate_nM`            | float | Tightest per-replicate fit. |
| `kd_max_replicate_nM`            | float | Loosest per-replicate fit. |
| `kd_replicate_cv_pct`            | float | CV% across replicate KDs. |
| `pkd_arith_mean`                 | float | `-log10(kd_arith_mean_nM_all / 1e9)`. **Use for stats.** |
| `weird_replicates_flag`          | bool  | True if any replicate had staircase / unexpected order / low confidence. |
| `assay_methods_mixed`            | bool  | True if SPR and BLI disagreed (13/100 designs). |
| `foundry_kd_nM_for_comparison`   | float | The upstream LIMS aggregate per design, kept for cross-checking. |
| `foundry_vs_recomputed_pct_diff` | float | Percent diff between the upstream aggregate and the one we recompute here. |

### 9. Homology / novelty (varies)

| column | type | notes |
|---|---|---|
| `mmseqs2_top1_target`           | string | Best mmseqs2 hit across SwissProt / SAbDab / etc. |
| `mmseqs2_top1_seqid`            | float  | Sequence identity of best hit. |
| `sim_top1_db`                   | string | Which database the best hit came from. |
| `sim_top1_target_id`            | string | Hit identifier in the source DB. |
| `sim_top1_similarity`           | float  | Same as `mmseqs2_top1_seqid` but normalised. |
| `sabdab_top1_target`            | string | Best SAbDab antibody hit (if any). |
| `sabdab_top1_seqid`             | float  | Identity to closest known antibody. |
| `ab_sim_target_id`              | string | Legacy antibody similarity target id, populated for 3 rows. |
| `final_modality_evidence_chain` | string | Audit string for how `modality` was decided. |
| `novelty_score`                 | float  | Composite novelty score (lower = more derivative). 140/141 populated. |

### 10. Forward-compatibility slots (currently empty)

These columns exist for downstream pipelines but are not yet populated.
Treat as TODO.

`foldseek_top1_pdb`, `foldseek_top1_tm`, `foldseek_top1_seqid`,
`foldseek_top1_organism`, `tm_top1_match`, `tm_score`, `annotation_class`,
`cath_class`, `cath_architecture`, `uniref_top1_id`, `uniref_top1_seqid`,
`n_domains`, `domain_lengths`, `cath_topology_per_domain`.

Don't drop them from the parquet — that breaks stable column ordering for
downstream consumers. Null-check on use.

### 11. ProteinBase enriched data (100/141 — screened only)

Mirrored from the ProteinBase public collection. The 100 screened
designs (`submitted_to_lab=True`) get the full bundle; the 41
non-screened designs are null across every `pb_*` column. Filter on
`pb_id.notna()` if you want only enriched rows.

The `pb_boltz2_*` columns are a **re-fold** of the cherry-pick
predictions. On the 100 overlapping designs they track the legacy
`boltz2_*` columns closely but not exactly (max |Δ ipSAE| = 0.18, max
|Δ iptm| = 0.085) — same Boltz-2 model, fresh seed. Use whichever you
prefer; document which in your analysis.

**Boltz-2 confidence (re-fold)**

| column | type | notes |
|---|---|---|
| `pb_boltz2_ipsae`          | float | Re-fold ipSAE. **Primary re-fold metric.** 100/141. |
| `pb_boltz2_iptm`           | float | Re-fold interface pTM. 100/141. |
| `pb_boltz2_ptm`            | float | Re-fold global pTM. 100/141. |
| `pb_boltz2_plddt`          | float | Re-fold mean pLDDT (binder chain). 100/141. |
| `pb_boltz2_complex_plddt`  | float | Re-fold complex-level pLDDT. 100/141. |
| `pb_boltz2_complex_iplddt` | float | Re-fold interface pLDDT. 100/141. |
| `pb_boltz2_complex_pde`    | float | Re-fold predicted distance error at the interface. 100/141. |
| `pb_boltz2_pdockq`         | float | Re-fold predicted DockQ. 100/141. |
| `pb_boltz2_pdockq2`        | float | Re-fold pDockQ2 (v2 formulation). 100/141. |
| `pb_boltz2_min_ipsae`      | float | Min ipSAE across re-scored ranks. 100/141. |
| `pb_boltz2_lis`            | float | Local Interaction Score from Boltz-2. 100/141. |

**ESMFold**

| column | type | notes |
|---|---|---|
| `pb_esmfold_plddt` | float | Mean pLDDT from a fresh single-chain ESMFold pass on the binder. 100/141. |

**ProteinMPNN**

| column | type | notes |
|---|---|---|
| `pb_proteinmpnn_score`            | float | ProteinMPNN log-likelihood of the submitted sequence given the fold. 100/141. |
| `pb_redesigned_proteinmpnn_score` | float | Score of the ProteinMPNN-redesigned sequence on the same backbone. 100/141. |
| `pb_proteinmpnn_seq_recovery`     | float | Fraction of positions where MPNN redesign agrees with the submitted sequence. 100/141. |

**Sequence / structure features**

| column | type | notes |
|---|---|---|
| `pb_molecular_weight`                       | float  | Daltons (ProteinBase recomputation). 100/141. |
| `pb_isoelectric_point`                      | float  | pI (ProteinBase recomputation). 100/141. |
| `pb_novelty`                                | float  | ProteinBase composite novelty score. 99/141. |
| `pb_ted_confidence`                         | float  | TED domain-annotation confidence. 98/141. |
| `pb_classification`                         | string | CATH-style class: `Mainly Alpha` / `Mainly Beta` / `Alpha Beta`. 98/141. |
| `pb_design_class`                           | string | ProteinBase modality call: `Miniprotein` / `Nanobody` / `scFv` / `Peptide` / `Other`. 100/141. |
| `pb_foldstring`                             | string | Compact secondary-structure string per residue (e.g. `HHHH…`). 98/141. |
| `pb_shape_complimentarity_boltz2_binder_ss` | float  | Shape complementarity of binder secondary structure to the target on the Boltz-2 complex. 100/141. |

**AFDB50 homology** (some binders had no AFDB50 hit)

| column | type | notes |
|---|---|---|
| `pb_tm_score_afdb50`    | float | Best TM-score against AFDB50. 89/141. |
| `pb_seqidentity_afdb50` | float | Sequence identity to the AFDB50 top hit. 89/141. |
| `pb_rmsd_afdb50`        | float | RMSD to the AFDB50 top hit. 89/141. |

**Interface**

| column | type | notes |
|---|---|---|
| `pb_interface_residues_target` | float | Interface residues on the TREM2 side. 100/141. |
| `pb_interface_residues_binder` | float | Interface residues on the binder side. 100/141. |
| `pb_interface_residues_total`  | float | Sum of both sides. 100/141. |

**Per-design wet-lab roll-up** (ProteinBase's own aggregation; cross-check against `kd_arith_mean_nM_all`)

| column | type | notes |
|---|---|---|
| `pb_n_binding_records` | float | Number of binding records in ProteinBase for this design. 100/141. |
| `pb_any_binding`       | bool  | True if any record reports binding. 100/141. |
| `pb_n_kd_records`      | float | Number of KD fits. 100/141. |
| `pb_n_bli_curves`      | float | BLI sensorgrams shipped under `data/proteinbase/sensorgrams/`. 100/141. |
| `pb_n_spr_curves`      | float | SPR sensorgrams shipped under `data/proteinbase/sensorgrams/`. 100/141. |
| `pb_kd_M_mean`         | float | Mean KD across replicates, in **molar** (not nM). 36/141. |
| `pb_kd_M_min`          | float | Tightest replicate KD, molar. 36/141. |
| `pb_kd_M_max`          | float | Loosest replicate KD, molar. 36/141. |

**Local artifact paths** (relative to repo root; populated for the 100 screened designs)

| column | type | notes |
|---|---|---|
| `pb_id`            | string | ProteinBase slug, e.g. `gentle-bear-dust`. 100/141. |
| `pb_boltz2_cif`    | string | `data/proteinbase/boltz2/design_NNN.cif`. 100/141. |
| `pb_esmfold_cif`   | string | `data/proteinbase/esmfold/design_NNN.cif`. 100/141. |
| `pb_pae_json`      | string | `data/proteinbase/pae/design_NNN.json`. 100/141. |
| `pb_stylized_png`  | string | `data/proteinbase/images/design_NNN.png`. 99/141. |

For sensorgrams (multiple per design, no single path per row), glob
`data/proteinbase/sensorgrams/design_NNN_rep*.json` instead.

## Joining with the lab tables

| from | to | key |
|---|---|---|
| `data/designs.csv` | `data/raw_lab/bli_replicates.csv` | `design_id` (1-to-many; 215 replicate rows across 100 designs) |
| `data/designs.csv` | `data/raw_lab/bli_results.csv`    | `design_id` (1-to-1; one row per screened design) |

The canonical KD columns in `designs.csv` already apply the
`selected & ~excluded & fixed` replicate filter. Re-aggregate from
`bli_replicates.csv` only if you need a different filter — confidence
strata, SPR-only, BLI-only — than the four we shipped, and document
which.

## What's NOT in this CSV

The CSV holds scalars only. The bulk binary artifacts ship alongside it
under `data/proteinbase/` (see [`data/proteinbase/README.md`](../data/proteinbase/README.md))
and are addressable via the `pb_boltz2_cif` / `pb_esmfold_cif` /
`pb_pae_json` / `pb_stylized_png` path columns. For the 100 screened
designs you get:

- **Boltz-2 complex CIFs** — `data/proteinbase/boltz2/` (100 files, 14 MB).
- **ESMFold binder CIFs** — `data/proteinbase/esmfold/` (100 files, 5 MB).
- **PAE matrices** — `data/proteinbase/pae/` (100 files, 88 MB).
- **Stylised renders** — `data/proteinbase/images/` (99 files, 27 MB).
- **Sensorgrams** — `data/proteinbase/sensorgrams/` (215 files: 193 SPR + 22 BLI, 27 MB).

Still not shipped:

- **AF3 / Chai-1 folds.** Only Boltz-2 + ESMFold are mirrored. AF3 / Chai
  predictions are regenerable from `designs.fasta`.
- **Per-residue pLDDT / ipSAE.** Only complex-level scalars are in the
  CSV; recompute from the shipped CIFs / PAEs if you need per-residue.
- **Raw hackathon submission CSVs.** Organiser-private; the pooled
  `designs.csv` is the public artifact.

## Regeneration

`data/designs.csv` is rebuilt by the upstream pipeline (sequence pooling
→ folding → ipSAE → BLI fetch → replicate aggregation → modality call).
See [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md). If you spot a wrong row,
**don't edit the CSV in place.** Open an issue with the `design_id` and
the corrected value.
