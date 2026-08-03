# Cross-predictor interface panel (SI) + Protenix-primary statistics

Predictors: **Protenix** (primary, 141 designs), **Chai** (141), **ESMFold2** (team-contributed complex panel, 100 screened designs). Descriptors identical to `jobs/interface_descriptors.py`; ESMFold2 read the same way (chain A = 175-aa target).

## Coverage & ESMFold2 interface plausibility

- **protenix**: 141 designs (89 expressed); BSA median 1521 Å²; interface-residue median 20; designs with zero interface contacts: 0.
- **chai**: 141 designs (89 expressed); BSA median 2087 Å²; interface-residue median 24; designs with zero interface contacts: 0.
- **esmfold2**: 100 designs (89 expressed); BSA median 1982 Å²; interface-residue median 24; designs with zero interface contacts: 0.

## Binder-vs-non-binder discrimination, per predictor (expressed designs)

| descriptor | predictor | binder med | non med | Cliff δ | MWU p | AUROC | AP | n(b/n) |
|---|---|---|---|---|---|---|---|---|
| bsa | protenix | 1823.0 | 1515.3 | 0.441 | 0.00042 | 0.72 | 0.622 | 37/52 | ⟵ primary
| bsa | chai | 2376.0 | 1816.8 | 0.353 | 0.0047 | 0.677 | 0.566 | 37/52 |
| bsa | esmfold2 | 2199.4 | 1807.0 | 0.376 | 0.0026 | 0.688 | 0.572 | 37/52 |
| n_iface_binder | protenix | 22.0 | 19.0 | 0.384 | 0.0021 | 0.692 | 0.561 | 37/52 |
| n_iface_binder | chai | 27.0 | 22.0 | 0.351 | 0.0049 | 0.675 | 0.544 | 37/52 |
| n_iface_binder | esmfold2 | 27.0 | 24.0 | 0.343 | 0.006 | 0.672 | 0.552 | 37/52 |
| n_iface_target | protenix | 17.0 | 14.5 | 0.285 | 0.022 | 0.642 | 0.553 | 37/52 |
| n_iface_target | chai | 22.0 | 19.0 | 0.272 | 0.03 | 0.636 | 0.526 | 37/52 |
| n_iface_target | esmfold2 | 22.0 | 18.0 | 0.242 | 0.053 | 0.621 | 0.495 | 37/52 |
| hbonds | protenix | 7.0 | 6.0 | 0.012 | 0.93 | 0.506 | 0.456 | 37/52 |
| hbonds | chai | 10.0 | 8.5 | 0.109 | 0.38 | 0.555 | 0.493 | 37/52 |
| hbonds | esmfold2 | 9.0 | 8.0 | 0.136 | 0.28 | 0.568 | 0.48 | 37/52 |
| salt_bridges | protenix | 2.0 | 1.5 | 0.22 | 0.073 | 0.61 | 0.53 | 37/52 |
| salt_bridges | chai | 3.0 | 2.0 | 0.128 | 0.3 | 0.564 | 0.515 | 37/52 |
| salt_bridges | esmfold2 | 2.0 | 2.0 | 0.136 | 0.26 | 0.568 | 0.51 | 37/52 |
| iface_frac_hydrophobic | protenix | 0.6 | 0.6 | 0.11 | 0.38 | 0.555 | 0.437 | 37/52 |
| iface_frac_hydrophobic | chai | 0.5 | 0.6 | -0.037 | 0.77 | 0.518 | 0.44 | 37/52 |
| iface_frac_hydrophobic | esmfold2 | 0.6 | 0.6 | 0.06 | 0.63 | 0.53 | 0.417 | 37/52 |
| iface_frac_polar | protenix | 0.1 | 0.2 | -0.121 | 0.33 | 0.561 | 0.44 | 37/52 |
| iface_frac_polar | chai | 0.1 | 0.2 | -0.148 | 0.24 | 0.574 | 0.461 | 37/52 |
| iface_frac_polar | esmfold2 | 0.1 | 0.1 | -0.038 | 0.76 | 0.519 | 0.436 | 37/52 |
| iface_frac_charged | protenix | 0.2 | 0.2 | 0.097 | 0.44 | 0.549 | 0.429 | 37/52 |
| iface_frac_charged | chai | 0.3 | 0.2 | 0.191 | 0.13 | 0.596 | 0.477 | 37/52 |
| iface_frac_charged | esmfold2 | 0.2 | 0.2 | 0.074 | 0.56 | 0.537 | 0.435 | 37/52 |
| hbond_density_per100 | protenix | 0.4 | 0.4 | -0.117 | 0.35 | 0.558 | 0.473 | 37/52 |
| hbond_density_per100 | chai | 0.4 | 0.4 | 0.018 | 0.89 | 0.509 | 0.442 | 37/52 |
| hbond_density_per100 | esmfold2 | 0.5 | 0.5 | -0.003 | 0.98 | 0.502 | 0.427 | 37/52 |

## Inter-predictor agreement & signed difference (BSA)

| pair | n | Spearman ρ | median signed Δ (A−B) | A>B % | Wilcoxon p |
|---|---|---|---|---|---|
| protenix − chai | 141 | 0.32 | -418 | 18% | 7.1e-18 |
| protenix − esmfold2 | 100 | 0.50 | -245 | 29% | 1.3e-08 |

## Protenix-primary interface size vs the published Protenix+Chai average

| set | binder med | non med | Cliff δ | MWU p | AUROC | AP |
|---|---|---|---|---|---|---|
| Protenix+Chai avg (published) | 2062.1 | 1741.8 | 0.421 | 0.00076 | 0.71 | 0.624 |
| Protenix only (new primary) | 1823.0 | 1515.3 | 0.441 | 0.00042 | 0.72 | 0.622 |

### Length controls — 'size, not length' (reconstructed; validate avg vs published 0.33/0.002, VE 0.0015, LR 0.006)

| set | partial ρ(BSA,hit|len) | p | van Elteren p (cohort) | van Elteren p (len-tertile) | LR test p |
|---|---|---|---|---|---|
| Protenix+Chai avg | 0.291 | 0.006 | 0.0002 | 0.0016 | 0.0060 |
| Protenix only | 0.333 | 0.001 | 0.0001 | 0.0004 | 0.0015 |
