# Glossary

Define every term that an external reader might trip over.

## Hackathon-specific

- **muni** — the platform that hosted the six AI agent teams
  end-to-end. Each agent runs an autonomous pipeline (sequence sampling
  + structure prediction + scoring) with tool access. "Agent" in the
  paper means "an end-to-end pipeline running on muni," not "a raw
  language model emitting sequences."
- **muni** — the partner platform that ran the hackathon alongside
  Adaptyv.
- **Cohort** — `human` (10 teams, 81 designs) or `agent` (6 teams, 60
  designs). Determined by the `is_human` column in `data/designs.csv`,
  never by parsing team names.
- **Top 100 cherry-pick** — the 100 designs ranked highest by Boltz-2
  ipSAE at submission time, sent to the BLI assay. The remaining 41
  designs were not screened wet-lab.
- **Foundry (method)** — the team-naming for the Baker lab's
  RFdiffusion3 pipeline. 17 designs across 5 agent teams
  (claude-sonnet-4.6, qwen-3.5-plus, Gemini 3.1 Pro, GLM 5, GPT 5.2).
  Buckets into the `RFDiffusion` family.

## Target & controls

- **TREM2** — Triggering Receptor Expressed on Myeloid cells 2. An
  Alzheimer's microglial receptor with clinically validated relevance.
  The construct used in the assay is Acro Biosystems catalogue
  `TR2-H52H5`: the ectodomain (residues 19–174 of UniProt Q9NZC2), an
  initiator Met, a GGGSGGGS linker, and a 10×His tag — 175 aa total,
  19.3 kDa. The exact sequence lives at
  `data/target/trem2_construct.fasta`. Designers fold against this
  whole construct, not the bare IgSF subdomain.
- **AL002** — Alector's anti-TREM2 antibody, failed Phase 2 in 2024.
  Used here as a sequence-level positive control (patent WO2019028346A1).
- **VHB937** — Novartis's anti-TREM2 antibody, ongoing Phase 2. Second
  positive control (patent WO2022122788A2).

## Metrics

- **`pKD`** — `-log10(KD in molar)`. Higher = tighter binder. Use for
  stats; report `KD` in nM/μM in prose for readability.
- **ipSAE** — predicted interface confidence score from Dunbrack lab
  (PDF: 2024). The canonical scorer is `scripts/utils/ipsae.py` (TODO:
  port the upstream implementation). We use Boltz-2 ipSAE for the
  cherry-pick.
- **pTM / ipTM** — predicted template-modelling score and interface pTM,
  AlphaFold2 family metrics (Jumper et al. 2021).
- **pLDDT** — per-residue model confidence (0–100) from AlphaFold2/ESMFold.
- **pDockQ** — predicted DockQ score, structure-level prediction quality
  for docked complexes (Bryant et al. 2022).
- **Modality** — protein archetype: `Peptide` (≤25 aa), `Miniprotein`
  (26–110 aa, non-antibody), `Large miniprotein` (>110 aa, non-antibody),
  `Nanobody` (single-VHH motif, 80–160 aa), `scFv` (VH+VL with G4S linker,
  often ≥200 aa).
- **Modality family** — coarser: `Peptide` / `Miniprotein` / `Antibody`.

## Method normalization

Self-reported `design_method` strings collapse to `design_method_normalized`
and `method_family` via the rules in `data/designs.csv`. Examples:

| raw | normalized | family |
|---|---|---|
| `RFDiffusion + LigandMPNN + Boltz` | `RFDiffusion+LigandMPNN+Boltz` | `RFDiffusion` |
| `mosaic` | `Mosaic` | `Mosaic` |
| `BoltzGen + Boltz2 + ipSAE` | `BoltzGen+Boltz2+ipSAE` | `BoltzGen` |
| `AF Hallucination + ProteinMPNN + Pyrosetta + Boltz2` | `AFHall+MPNN+PyR+Boltz2` | `Hallucination` |
| `Foundry` (any) | `Foundry` | `RFDiffusion` (Baker lab RFdiffusion3) |
| anything else | (verbatim) | `Other` |

## Binding labels

Adaptyv's binding-label hierarchy, finest → coarsest:

```
binder > strong > medium > weak > potential_binder > non_binder > no_expression > unknown
```

- `is_hit` = `binding_label ∈ {binder, strong, medium, weak}`.
- `potential_binder` is borderline and reported separately.
- `non_binder` and `no_expression` are distinct failure modes — track both.
- **In this dataset only four labels appear:** `binder` (36), `weak` (1),
  `non_binder` (52), `no_expression` (11). The curator did not use the
  `strong` / `medium` / `potential_binder` / `unknown` tiers for this
  campaign. The hierarchy is documented as a stable contract; analyses
  should not assume the missing tiers will stay missing in future runs.

## Assay terms

- **BLI** — Bio-Layer Interferometry. Gator-family instruments.
- **SPR** — Surface Plasmon Resonance. Carterra-family instruments.
- **Replicate** — one row in `bli_replicates.csv`. Multiple per design.
- **`selected & ~excluded & fixed`** — the canonical replicate filter.
  "Pushed and confirmed" replicates only — everything else stays in
  the per-replicate file for transparency but doesn't enter the
  per-design KD.
- **Sensorgram** — the raw response curve from BLI/SPR. Not shipped in
  this repo; the per-replicate fits are the public artifact.

## Reproducibility terms

- **`design_id`** — the canonical row id, 1..141. Stable across the
  pipeline. **Use this as the join key everywhere.**
