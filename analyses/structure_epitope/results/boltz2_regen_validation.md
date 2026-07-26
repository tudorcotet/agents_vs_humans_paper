# Boltz-2 regeneration — validation (RunPod, with MSA, kernels-on)

Regenerated all 141 Boltz-2 complexes on RunPod (A40): Boltz-2 v2.2.1, **precomputed target MSA**
(target-only, `--use_msa_server` MSA computed once and reused; de novo binder single-seq),
cuEquivariance kernels on, single model, seed 42. Structures: `data/structures_regen/boltz2/`;
scalars `data/structures_regen/boltz2_regen_scalars.csv`.

**Integrity: 141/141 binder chains match the canonical sequence** (vs. 10/141 in the shipped corrupt
files) — the mapping bug is gone.

**Corrected `b2_iptm` agreement (Spearman):** competition ipSAE **0.80**, legacy `boltz2_iptm` **0.83**,
Protenix **0.62**, Chai **0.53** — but only **0.11** vs. the *corrupt* `b2_iptm`. The corrected values
correlate with every clean reference; the corrupt ones correlated with none.

**Binder discrimination (P_expressed, AUROC of iptm):** corrected `b2_iptm` **0.65** (corrupt was
**0.54** = random; Protenix 0.79, Chai 0.77). The regenerated Boltz-2 is a real, non-random predictor
again.

**Bearing on the section:** none of the section's results used the corrupt `b2_*` (the structural
analysis used Protenix+Chai; the affinity correlations used competition `submitted_ipsae` + legacy
`boltz2_*`, both clean). The regenerated Boltz-2 provides an additional clean predictor and can be
folded into the epitope consensus as a third method. This regeneration also **fixes the data package**
for release (corrected Boltz-2 structures to replace the mislabeled ones).

## Epitope robustness on the corrected Boltz-2 (3rd predictor)
Per-residue TREM2 contact frequency on the corrected Boltz-2 agrees with the Protenix+Chai consensus
at **Spearman ρ=0.78** (both expressed and binder sets). Its top-10 binder-contacted residues are all
apical and include **all six hydrophobic-tip residues (44/69/70/71/74/89)**. The epitope-convergence
finding therefore **holds across three independent predictors** (Protenix, Chai, corrected Boltz-2),
in addition to the EvoEF2 energetic hotspots — a strong robustness result.
(`results/epitope_boltz2_regen_freq.csv`; `jobs/regen_boltz2_footprints.py`.)
