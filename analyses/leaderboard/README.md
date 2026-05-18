# `analyses/leaderboard/` — per-design top-N tables

The leaderboard analysis. Drives Figure 1 in the paper (the headline
ranking) and several supplementary tables.

## What it produces

- `top20_overall.csv` — the 20 tightest-binding designs across all teams.
- `top10_human.csv` / `top10_agent.csv` — per-cohort top-10.
- `method_winners.csv` — best design per `design_method_normalized`.
- `team_winners.csv` — best design per `team`.

All tables sort by `kd_arith_mean_nM_all` ascending (tightest first).

## Inputs

- `data/designs.csv` via `scripts.utils.load_designs()`.

## Run

```bash
mise run analysis:leaderboard
```

## Columns

The output tables carry a stable column set so figure code can read any
of them interchangeably:

```
rank, design_id, name, team, is_human,
design_method, design_method_normalized, sequence_length,
kd_arith_mean_nM_all, pkd,
n_replicates_pushed, weird_replicates_flag,
submitted_ipsae, binding_label, binding_strength
```

Every design in the leaderboard is in the screened set (`submitted_to_lab=True`),
so it carries the full `pb_*` family — re-fold Boltz-2 confidence
(`pb_boltz2_ipsae`, `pb_boltz2_iptm`, `pb_boltz2_pdockq2`,
`pb_boltz2_plddt`, …), ProteinMPNN scores (`pb_proteinmpnn_score`,
`pb_proteinmpnn_seq_recovery`), interface counts
(`pb_interface_residues_{target,binder,total}`), shape complementarity
(`pb_shape_complimentarity_boltz2_binder_ss`), and the local artifact
paths (`pb_boltz2_cif`, `pb_esmfold_cif`, `pb_pae_json`, `pb_stylized_png`).
KDs are also available in molar via `pb_kd_M_{mean,min,max}` (36/141, same
coverage as `kd_arith_mean_nM_all`). See [`../../docs/DATA.md` §11](../../docs/DATA.md#11-proteinbase-enriched-data-100141--screened-only).
