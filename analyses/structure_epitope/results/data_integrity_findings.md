# Data-integrity findings — corrupt Boltz-2 / AF2M artifacts (for the release decision)

Diagnostic: `jobs/diagnose_scalars.py` (+ `structure_integrity_audit.csv`). Definitive.

## What is corrupt
The Modal-rerun Boltz-2 and AF2-Multimer artifacts are internally mislabeled: the **complex
structure files** (`data/structures/boltz2/`, `data/structures/af2m/`) carry the wrong binder chain
in most designs (correct in only 10/141 and 15/141), **and the derived scalar columns were computed
from those mislabeled files, so they are corrupt too:**
- **`b2_*` (Boltz-2 Modal rerun):** in all 10 groups of designs that share one mislabeled binder,
  `b2_iptm` is **identical** (16 designs all 0.96, 12 all 0.952, …); `b2_iptm` vs the competition
  `submitted_ipsae` ρ=**0.02**, vs Protenix ρ=−0.07, vs Chai ρ=0.01; binder AUROC **0.54** (~random).
- **`af2m_*`:** same pattern — 14 identical-value groups; vs `submitted_ipsae` ρ=**−0.00**; AUROC 0.57.
- **Consequently the consensus columns `ipsae_pass_4folders` / `iptm_pass_4folders` are partly
  corrupt** (they average in b2 and af2m).

## What is CLEAN (verified)
- **Competition / legacy scalars:** `submitted_ipsae`, and legacy `boltz2_iptm` / `boltz2_plddt` /
  `boltz2_complex_plddt` (from the muni selection table) — `boltz2_iptm` vs `submitted_ipsae` ρ=**0.97**,
  varied (not identical) within mislabel-groups → CLEAN. These are what the competition ranked on.
- **`pb_boltz2_*`** (ProteinBase mirror, screened 100): vs `submitted_ipsae` ρ=**0.87** → CLEAN. (A
  correct Boltz-2 scalar set already exists for the screened designs.)
- **`px_*` (Protenix), `chai_*` (Chai-1)** and their structure files: CLEAN (px vs chai iptm ρ=0.67;
  px vs competition ρ=0.72).

## Impact
- **This section: none.** All structural analysis used Protenix + Chai (clean); the Phase-6 metric
  correlations used `submitted_ipsae` and legacy `boltz2_*` (clean), never `b2_*`/`af2m_*`.
- **Broader paper:** any analysis/figure using `b2_*`, `af2m_*`, or `*_pass_4folders` must be
  re-checked or dropped. The competition ranking (submitted_ipsae) is unaffected.

## Options
1. **Scalars (do now, $0):** drop `b2_*`, `af2m_*`, `*_pass_4folders` from analysis; where a Boltz-2
   scalar is wanted, use legacy `boltz2_*` or `pb_boltz2_*`. Add a one-line note to the data docs.
2. **Structure files (before release):**
   - **(B) Quarantine + document** — remove the corrupt `boltz2/`+`af2m/` complex files from the
     public package; state only Protenix/Chai complexes are shipped, Boltz-2/AF2M regenerable from
     `designs.fasta`. Honest minimum, $0.
   - **(C) Regenerate** Boltz-2 (and AF2M) complexes correctly (GPU/RunPod). Complete but costs
     compute; worth it only if shipped Boltz-2 structures are needed (the competition ranked on
     Boltz-2, so a correct set has value — though `pb_boltz2_*` already covers the scalars).
3. **Flag upstream:** the corruption originates in the muni/Adaptyv Modal-rerun/mirror step, upstream
   of this repo — report to the data-package owner so the source is fixed, not just this copy.
