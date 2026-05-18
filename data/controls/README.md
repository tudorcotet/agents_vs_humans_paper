# Positive controls (TREM2 hackathon)

Two clinical-stage anti-TREM2 antibodies are referenced as assay sanity checks: **AL002** (Alector, patent WO2019028346A1, failed Phase 2 2024) and **VHB937** (Novartis, patent WO2022122788A2, ongoing Phase 2).

## Status

`known_binders.fasta` currently contains placeholder records. The variable-domain sequences need to be extracted from the patent SEQ IDs by hand; the identity-to-known-binders comparison is therefore **pending**. Once real sequences land, regenerate this FASTA and re-run `analyses/sequence_diversity/diversity.py`; the script will produce `data/processed/identity_to_known_binders.csv` automatically.
