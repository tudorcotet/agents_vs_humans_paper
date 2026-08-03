# Partial affinity correlations controlling for BSA — P_kd, n=36

Spearman ρ of each descriptor with pKd, raw and after partialling out **buried area (BSA)** (phase6's `partial_spearman`; matches the `partial_rho_bsa` column of `metric_correlations.csv`). Permutation p on the partial (10,000) uses a **separate RNG stream**, so the committed FDR family in `metric_correlations.csv` is unchanged.

| descriptor | raw ρ | partial ρ (\| BSA) | perm p |
|---|---|---|---|
| PRODIGY ΔG | +0.523 | +0.224 | 0.197 |
| interface size (residues) | +0.473 | +0.155 | 0.374 |
| H-bonds | +0.252 | -0.028 | 0.874 |
| salt bridges | +0.072 | -0.072 | 0.685 |
| hydrophobic fraction | -0.305 | -0.171 | 0.332 |

**Symmetric (reverse) test — BSA \| PRODIGY:** partial ρ = +0.147 (perm p 0.406). Both directions are small and n.s. (raw Spearman(BSA, PRODIGY) = 0.822), so PRODIGY and BSA are **one signal, not two** — an equivalence, not a reduction of PRODIGY to BSA.
**BSA vs length:** Spearman(BSA, length) = 0.190 — interface size and binder length are two largely independent affinity signals.
