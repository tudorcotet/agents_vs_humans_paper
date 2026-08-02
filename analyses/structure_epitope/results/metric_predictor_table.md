# Per-predictor metric table — AUROC, ρ, and FDR

For every metric named in the Summary. **AUROC** = binder-vs-non-binder discrimination on P_expressed (n=89; 37 binders/52 non); **ρ** = Spearman vs pKd on P_kd (n=36); **q** = BH-FDR within tier, per analysis (this is a broader family than the per-section panels, so these q's run slightly more conservative than the draft's headline values — e.g. Protenix BSA q=0.014 here vs 0.004 in §2's 9-descriptor panel; the raw p-values and effect sizes are identical). `*` marks q<0.10. Cells are `—` where the metric does not exist for that predictor (the axis is ragged by construction — see header of `jobs/metric_predictor_table.py`).

Predictors: **Protenix** is the section primary. Interface geometry + secondary structure are on Protenix/Chai/ESMFold2/Boltz-2/AF2M (Boltz-2 + AF2M = regenerated/PR#6 complexes; AF2M ships as PDB); learned confidence + PRODIGY on Boltz-2/Protenix/Chai/AF2M; ESM-2/SaProt are sequence-based (one value); **contact motifs are population enrichment (odds ratios), not a per-design score, so they are reported separately below, not as AUROC/ρ.** pLDDT is **mean-complex pLDDT** on 0–1 for all three (Protenix rescaled from 0–100) — a slightly different quantity from §3's Boltz-2 *binder* pLDDT (ρ=0.14).

## Binder-vs-non-binder discrimination — AUROC (q)

| tier | metric | Boltz-2 | Protenix | Chai | ESMFold2 | AF2M | sequence |
|---|---|---|---|---|---|---|---|
| structural/physical | buried area (BSA) | 0.65* (0.095) | 0.72* (0.025) | 0.68* (0.048) | 0.69* (0.039) | 0.70* (0.037) | — |
| structural/physical | binder iface residues | 0.62 (0.222) | 0.69* (0.039) | 0.68* (0.048) | 0.67* (0.050) | 0.64 (0.128) | — |
| structural/physical | target iface residues | 0.60 (0.405) | 0.64 (0.128) | 0.64 (0.145) | 0.62 (0.222) | 0.66* (0.087) | — |
| structural/physical | hydrophobic frac | 0.53 (0.768) | 0.56 (0.686) | 0.52 (0.874) | 0.53 (0.768) | 0.51 (0.959) | — |
| structural/physical | polar frac | 0.56 (0.679) | 0.56 (0.679) | 0.57 (0.668) | 0.52 (0.874) | 0.51 (0.902) | — |
| structural/physical | charged frac | 0.55 (0.731) | 0.55 (0.731) | 0.60 (0.415) | 0.54 (0.768) | 0.56 (0.668) | — |
| structural/physical | H-bond count | 0.54 (0.731) | 0.51 (0.959) | 0.56 (0.686) | 0.57 (0.668) | 0.54 (0.732) | — |
| structural/physical | salt bridges | 0.57 (0.668) | 0.61 (0.286) | 0.56 (0.668) | 0.57 (0.668) | 0.60 (0.415) | — |
| structural/physical | H-bond density | 0.50 (0.997) | 0.56 (0.679) | 0.51 (0.956) | 0.50 (0.997) | 0.52 (0.902) | — |
| structural/physical | interface helix frac | 0.58 (0.556) | 0.53 (0.787) | 0.57 (0.668) | 0.53 (0.768) | 0.55 (0.731) | — |
| structural/physical | interface strand frac | 0.54 (0.679) | 0.52 (0.768) | 0.52 (0.768) | 0.53 (0.731) | 0.52 (0.768) | — |
| structural/physical | PRODIGY pKd | 0.58 (0.617) | 0.53 (0.768) | 0.53 (0.768) | — | 0.55 (0.731) | — |
| learned/sequence | ipTM | 0.60 (0.111) | 0.79* (0.000) | 0.77* (0.000) | — | 0.65* (0.026) | — |
| learned/sequence | pLDDT (mean) | 0.63* (0.048) | 0.71* (0.002) | 0.74* (0.000) | — | 0.64* (0.033) | — |
| learned/sequence | ipSAE | 0.65* (0.026) | 0.78* (0.000) | 0.78* (0.000) | — | 0.66* (0.023) | — |
| learned/sequence | ESM-2 pll | — | — | — | — | — | 0.52 (0.705) |
| learned/sequence | SaProt pll | — | — | — | — | — | 0.55 (0.416) |

## Affinity (pKd) correlation — ρ (q)

| tier | metric | Boltz-2 | Protenix | Chai | ESMFold2 | AF2M | sequence |
|---|---|---|---|---|---|---|---|
| structural/physical | buried area (BSA) | +0.20 (0.404) | +0.50* (0.027) | +0.43* (0.068) | -0.00 (0.983) | +0.35 (0.180) | — |
| structural/physical | binder iface residues | +0.32 (0.209) | +0.47* (0.042) | +0.52* (0.027) | +0.04 (0.909) | +0.45* (0.062) | — |
| structural/physical | target iface residues | +0.09 (0.705) | +0.27 (0.286) | +0.39 (0.127) | -0.25 (0.290) | +0.51* (0.027) | — |
| structural/physical | hydrophobic frac | -0.19 (0.417) | -0.30 (0.209) | -0.36 (0.179) | -0.32 (0.209) | -0.31 (0.209) | — |
| structural/physical | polar frac | +0.30 (0.213) | +0.26 (0.289) | +0.28 (0.239) | +0.32 (0.209) | +0.34 (0.209) | — |
| structural/physical | charged frac | +0.03 (0.926) | +0.12 (0.638) | +0.21 (0.396) | +0.12 (0.638) | +0.07 (0.784) | — |
| structural/physical | H-bond count | +0.15 (0.531) | +0.25 (0.290) | +0.33 (0.209) | +0.01 (0.983) | +0.20 (0.404) | — |
| structural/physical | salt bridges | -0.09 (0.705) | +0.07 (0.765) | +0.08 (0.757) | -0.02 (0.964) | +0.12 (0.638) | — |
| structural/physical | H-bond density | +0.12 (0.638) | +0.20 (0.404) | +0.25 (0.290) | -0.01 (0.983) | +0.11 (0.649) | — |
| structural/physical | interface helix frac | -0.20 (0.404) | -0.18 (0.417) | -0.25 (0.290) | -0.18 (0.417) | -0.31 (0.209) | — |
| structural/physical | interface strand frac | +0.24 (0.290) | +0.24 (0.290) | +0.17 (0.472) | +0.24 (0.290) | +0.34 (0.209) | — |
| structural/physical | PRODIGY pKd | -0.28 (0.239) | +0.52* (0.027) | +0.42* (0.075) | — | -0.10 (0.705) | — |
| learned/sequence | ipTM | +0.28 (0.177) | +0.20 (0.332) | +0.35 (0.124) | — | +0.15 (0.461) | — |
| learned/sequence | pLDDT (mean) | +0.36 (0.124) | +0.28 (0.177) | +0.42 (0.115) | — | +0.20 (0.332) | — |
| learned/sequence | ipSAE | +0.40 (0.115) | +0.29 (0.177) | +0.34 (0.124) | — | +0.19 (0.332) | — |
| learned/sequence | ESM-2 pll | — | — | — | — | — | +0.09 (0.661) |
| learned/sequence | SaProt pll | — | — | — | — | — | +0.06 (0.733) |

## Signed difference from Protenix (reference) — median(Protenix − other), per design

Positive = Protenix reads *higher* than the other predictor. `rel%` = median diff ÷ Protenix median. `Prot>other%` = fraction of designs where Protenix is higher. Wilcoxon signed-rank p (paired, all designs with both predictors).

| metric | vs predictor | n | median Δ (Prot−other) | rel% | Prot>other% | Wilcoxon p |
|---|---|---|---|---|---|---|
| buried area (BSA) | Chai | 141 | -418 | -27.5% | 18% | 7.1e-18 |
| buried area (BSA) | ESMFold2 | 100 | -244.9 | -15.1% | 29% | 1.3e-08 |
| buried area (BSA) | Boltz-2 | 141 | -152.7 | -10% | 31% | 7.8e-09 |
| buried area (BSA) | AF2M | 141 | -14.5 | -1% | 48% | 9.3e-01 |
| binder iface residues | Chai | 141 | -3 | -15% | 16% | 3.5e-15 |
| binder iface residues | ESMFold2 | 100 | -2 | -9.5% | 23% | 3.3e-08 |
| binder iface residues | Boltz-2 | 141 | -2 | -10% | 30% | 5.4e-07 |
| binder iface residues | AF2M | 141 | -1 | -5% | 38% | 3.3e-02 |
| target iface residues | Chai | 141 | -4 | -25% | 23% | 1.5e-12 |
| target iface residues | ESMFold2 | 100 | -3 | -18.8% | 33% | 3.0e-05 |
| target iface residues | Boltz-2 | 141 | -2 | -12.5% | 33% | 8.7e-06 |
| target iface residues | AF2M | 141 | +0 | +0% | 47% | 9.7e-01 |
| hydrophobic frac | Chai | 141 | +0.012 | +2.1% | 55% | 3.9e-02 |
| hydrophobic frac | ESMFold2 | 100 | +0 | +0% | 48% | 5.6e-01 |
| hydrophobic frac | Boltz-2 | 141 | +0.006 | +1.1% | 50% | 9.1e-02 |
| hydrophobic frac | AF2M | 141 | -0.005 | -0.9% | 44% | 2.1e-01 |
| polar frac | Chai | 141 | -0.005 | -3.2% | 36% | 1.8e-02 |
| polar frac | ESMFold2 | 100 | -0.001 | -0.7% | 45% | 4.0e-01 |
| polar frac | Boltz-2 | 141 | -0.005 | -3.2% | 38% | 1.7e-02 |
| polar frac | AF2M | 141 | -0.007 | -4.4% | 38% | 8.7e-02 |
| charged frac | Chai | 141 | +0 | +0% | 44% | 5.2e-01 |
| charged frac | ESMFold2 | 100 | +0.005 | +2.3% | 52% | 9.4e-01 |
| charged frac | Boltz-2 | 141 | +0 | +0% | 49% | 9.5e-01 |
| charged frac | AF2M | 141 | +0.017 | +7.3% | 58% | 2.9e-03 |
| H-bond count | Chai | 141 | -1 | -14.3% | 40% | 1.1e-03 |
| H-bond count | ESMFold2 | 100 | -1 | -14.3% | 36% | 4.3e-03 |
| H-bond count | Boltz-2 | 141 | +0 | +0% | 45% | 2.9e-01 |
| H-bond count | AF2M | 141 | -1 | -14.3% | 33% | 5.3e-04 |
| salt bridges | Chai | 141 | -1 | -100% | 28% | 1.3e-04 |
| salt bridges | ESMFold2 | 100 | +0 | +0% | 27% | 8.9e-01 |
| salt bridges | Boltz-2 | 141 | +0 | +0% | 23% | 9.4e-03 |
| salt bridges | AF2M | 141 | +0 | +0% | 23% | 3.1e-02 |
| H-bond density | Chai | 141 | +0.051 | +11.6% | 55% | 2.5e-01 |
| H-bond density | ESMFold2 | 100 | -0.014 | -3.4% | 47% | 3.6e-01 |
| H-bond density | Boltz-2 | 141 | +0.074 | +16.8% | 56% | 4.7e-02 |
| H-bond density | AF2M | 141 | -0.088 | -20% | 40% | 1.5e-04 |
| interface helix frac | Chai | 141 | +0.008 | +0.9% | 53% | 1.0e-04 |
| interface helix frac | ESMFold2 | 100 | +0 | +0% | 47% | 2.0e-01 |
| interface helix frac | Boltz-2 | 141 | +0 | +0% | 41% | 3.2e-01 |
| interface helix frac | AF2M | 141 | +0 | +0% | 30% | 3.9e-02 |
| interface strand frac | Chai | 141 | +0 | — | 14% | 9.8e-01 |
| interface strand frac | ESMFold2 | 100 | +0 | — | 8% | 7.8e-01 |
| interface strand frac | Boltz-2 | 141 | +0 | — | 12% | 7.9e-01 |
| interface strand frac | AF2M | 141 | +0 | — | 15% | 3.4e-01 |
| PRODIGY pKd | Chai | 141 | -0.694 | -11.8% | 28% | 7.0e-12 |
| PRODIGY pKd | Boltz-2 | 141 | +0.774 | +13.1% | 72% | 5.0e-10 |
| PRODIGY pKd | AF2M | 141 | +1.056 | +17.9% | 78% | 1.3e-11 |
| ipTM | Chai | 141 | +0.303 | +38.4% | 99% | 3.3e-24 |
| ipTM | Boltz-2 | 136 | -0.077 | -9.6% | 24% | 2.6e-13 |
| ipTM | AF2M | 141 | +0.056 | +7.1% | 72% | 5.0e-07 |
| pLDDT (mean) | Chai | 141 | +0.113 | +13.8% | 100% | 6.9e-25 |
| pLDDT (mean) | Boltz-2 | 141 | +0.043 | +5.2% | 92% | 1.1e-23 |
| pLDDT (mean) | AF2M | 141 | +0.066 | +8% | 99% | 7.1e-25 |
| ipSAE | Chai | 141 | +0.167 | +32% | 91% | 1.5e-19 |
| ipSAE | Boltz-2 | 136 | -0.126 | -23.1% | 16% | 1.7e-15 |
| ipSAE | AF2M | 141 | +0.012 | +2.2% | 58% | 3.4e-01 |

## Contact motifs (not a per-design AUROC/ρ metric)

Level-2/3 typed contact-pair analysis is population enrichment, not a per-design score: aromatic-to-apex OR 3.1 (p=0.025), Arg/Lys-to-tip OR 2.3 (p=0.085), **none surviving FDR**; no convergent residue-pair motif beyond universal apical engagement (`results/motifs_L2_contactpairs.csv`, `motifs_L3_convergent.csv`). Computed on Protenix+Chai.

