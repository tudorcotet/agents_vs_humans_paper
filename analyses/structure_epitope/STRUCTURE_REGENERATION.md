# Corrupt Boltz-2 / AF2-Multimer structures: what happened, how we fixed it, and manuscript impact

**Author:** Structure & Epitope Analysis (Amy's section) · **Date:** 2026-07-25
**Scope:** the scrambled binder chains in `data/structures/{boltz2,af2m}/`, their derived scalar
columns, the RunPod regeneration, and the effect (if any) on the manuscript.
**TL;DR:** two of the four predictors' shipped complex files had the **wrong binder chain** for most
designs, and their derived metric columns were corrupt with them. We regenerated both correctly on
GPU (Boltz-2 with MSA; AF2-Multimer via ColabFold) and verified the fix. **The manuscript's current
figures and the Structure & Epitope section are unaffected** — they use the clean competition/legacy
metrics and the clean Protenix/Chai structures. The corrupt files **and their derived `b2_*`/`af2m_*`
columns have now been replaced** with the regenerated structures and recomputed values (§3.5).

---

## 1. What was wrong

Each shipped complex file has two chains: **chain A = the 175-aa TREM2 target** and **chain B = the
designed binder**. A per-file audit (matching every chain B against the submitted sequence in
`designs.parquet`; `jobs/diagnose_scalars.py`, `results/structure_integrity_audit.csv`) found:

| Predictor | Correct binder chain | Verdict |
|---|---|---|
| **Protenix** | **141 / 141** | clean |
| **Chai-1** | **141 / 141** | clean |
| **Boltz-2** (`data/structures/boltz2/`) | **10 / 141** | **corrupt** |
| **AF2-Multimer** (`data/structures/af2m/`) | **15 / 141** | **corrupt** |

In the corrupt files the **target chain A is correct everywhere**; only chain B is wrong. The wrong
binders are duplicated many-to-one onto roughly designs 1–10 (e.g. 16 different designs' Boltz-2 files
all contain design 1's binder). This is a **binder-to-design mapping bug in the upstream muni/Adaptyv
Modal re-fold step**, not a modelling error and not a problem with the design data itself.

**The derived scalar columns are corrupt too**, because they were computed *from* the mislabeled
files:
- `b2_*` (Boltz-2) and `af2m_*` (AF2M) columns in `designs.csv`/`.parquet`. Evidence: designs whose
  files share one wrong binder have **identical** `b2_iptm` (16 rows all 0.960, 12 all 0.952, …);
  `b2_iptm` correlates with the competition metric `submitted_ipsae` at only **ρ = 0.02** and is a
  **random** binder classifier (AUROC 0.54). AF2M is the same (ρ ≈ 0.00, AUROC 0.57).
- The two consensus columns **`ipsae_pass_4folders` / `iptm_pass_4folders`** *count*, across
  {b2, px, chai, af2m}, how many have `ipsae_d0chn_max ≥ 0.4` / `iptm ≥ 0.7` — so with `b2`/`af2m`
  corrupt they were **partly corrupt**.

### What was NOT affected (verified clean)
- **`submitted_ipsae`** — the competition ranking metric (upstream muni pipeline).
- **Legacy `boltz2_iptm` / `boltz2_plddt` / `boltz2_complex_plddt`** — from the muni selection table;
  `boltz2_iptm` vs `submitted_ipsae` ρ = **0.97**.
- **`pb_boltz2_*`** — the ProteinBase mirror's Boltz-2 rerun; vs `submitted_ipsae` ρ = **0.87**.
- **`px_*` (Protenix), `chai_*` (Chai-1)**, and the Protenix/Chai/ESMFold structure files.

---

## 2. How we regenerated them

Everything ran on RunPod (single A40 GPU, community/secure), scripted and reproducible
(`infra/runpod_setup.sh`, `infra/prod_boltz.py`, `infra/af2m_input.csv`), with a **binder-chain
integrity check on every output** so a repeat of the mapping bug would be caught immediately.

**Boltz-2 (v2.2.1), all 141:**
- Config **with MSA** (correcting the manuscript's "run without MSA" note): the TREM2 target chain gets
  a real MSA via the ColabFold MMseqs2 server; the de novo binder is single-sequence (`msa: empty`) —
  the common practice for binder co-folding.
- **MSA precomputed once and reused.** Because chain A is identical for all 141, the target MSA is
  computed a single time and passed to every fold (`infra/target.a3m`); the GPU then does only the
  fold, not 141 redundant server queries. (One gotcha fixed: boltz's mmseqs output contained stray
  null bytes that crashed the a3m parser — stripped with `tr -d '\000'`.)
- cuEquivariance kernels on (Boltz-2's default optimized path), single model, seed 42.

**AF2-Multimer, all 141:** ColabFold `colabfold_batch`, AF2-multimer, single model, 3 recycles, on
GPU (sharded across three pods to parallelize).

Corrected structures land in `data/structures_regen/{boltz2,af2m}/` and are intended to **replace** the
corrupt `data/structures/{boltz2,af2m}/` in the released package.

---

## 3. Validation of the fix (both predictors)

**AF2-Multimer, all 141** (`data/structures_regen/af2m/`, `af2m_regen_scalars.csv`): **integrity
141/141** binder chains match (was 15/141); corrected `af2m_iptm` agrees with the clean refs
(competition ipSAE **0.71**, Protenix **0.70**, Chai **0.67**) vs. only **0.11** with the corrupt
`af2m_iptm`; binder discrimination AUROC **0.65** (corrupt 0.57). Same clean result as Boltz-2.

### Boltz-2

- **Integrity: 141 / 141** binder chains match the submitted sequence (was 10/141). Mapping bug gone.
- **Corrected `b2_iptm` now agrees with every clean reference** (Spearman): competition `submitted_ipsae`
  **0.80**, legacy `boltz2_iptm` **0.83**, Protenix **0.62**, Chai **0.53** — but only **0.11** vs. the
  old *corrupt* `b2_iptm`, confirming the old values were noise.
- **Discriminates binders again:** AUROC **0.65** (corrupt was 0.54 = random; Protenix/Chai 0.77–0.79).
- **Strengthens the epitope result:** the corrected Boltz-2's per-residue contact map agrees with the
  Protenix+Chai consensus at **ρ = 0.78**, with all six hydrophobic-tip residues (W44, L69, W70, L71,
  F74, L89) in its top-10 — so the epitope-convergence finding now holds across **three independent
  structure predictors** (Protenix, Chai, corrected Boltz-2) plus the EvoEF2 energetic hotspots.
  (`results/boltz2_regen_validation.md`.)

## 3.5 Recomputing the derived `b2_*` / `af2m_*` columns

The 19 `b2_*` and 19 `af2m_*` columns (plus the two `*_pass_4folders` consensus columns) in
`data/designs.{csv,parquet}` were recomputed from the regenerated structures and their PAE, using the
repo's **own** upstream scorer — `scripts/modal/modal_boltz2_avh.compute_ipsae` (`PAE_CUTOFF=15`,
`DIST_CUTOFF=8`) imported directly, so there is **zero drift** from the code that produced the clean
`px_*`/`chai_*` columns. Script: `jobs/recompute_scalars.py` (idempotent; run with `--apply`). Only the
40 affected columns are touched — every other column round-trips verbatim.

Validation (`jobs/recompute_scalars.py` prints all of this):
- **iptm reproduces exactly**: recomputed `b2_iptm`/`af2m_iptm` match the independently-extracted native
  confidence iptm to **max |Δ| = 1×10⁻¹⁶ / 0**.
- **Agrees with the clean references** (Spearman): `b2_iptm` vs `submitted_ipsae` **0.80**, legacy
  `boltz2_iptm` **0.83**, `px_iptm` **0.62**, `chai_iptm` **0.53**; `af2m_iptm` vs those **0.71 / 0.68 /
  0.70 / 0.67**. Against the **old corrupt** columns ρ ≈ 0.0–0.11 (i.e. the values genuinely changed).
- **Discriminates binders again**: AUROC 0.65 (`iptm`) / 0.66–0.73 (`ipsae_d0chn_max`); corrupt was ~0.54.

The consensus counts shift accordingly (e.g. `ipsae_pass_4folders` now spans 0–4 instead of 1–4). Record
of the recomputed values: `results/recomputed_scalars.csv`.

---

## 4. Impact on the manuscript

**None on the load-bearing content.** Two independent reasons:

1. **The Structure & Epitope Analysis section never used the corrupt data.** Its structural analysis
   (footprints, interface descriptors, fold clustering, motifs) ran on **Protenix + Chai** (clean); its
   metric-vs-affinity analysis used **`submitted_ipsae` and legacy `boltz2_*`** (clean), never `b2_*`
   or `af2m_*`. The corrected Boltz-2 is now added as a confirming third predictor.

2. **The manuscript's existing figures use the clean competition/legacy metrics, not the corrupt
   rerun columns.** Figure 7 (metric ROC) and the in-silico-score table draw on `submitted_ipsae`,
   legacy `boltz2_iptm`/`boltz2_plddt`, and `esm_pll` — all verified clean. Figures 1, 2, 4, 5, 6 do
   not use the complex structures at all.

### Proof by reproduction (not assertion)

We reproduced the manuscript's headline structure-metric table ("In Silico Scores: Humans vs.
Agents") from the data and checked which column matches. **Every published value reproduces to the
digit from the CLEAN column, and the CORRUPT column gives different, non-significant numbers**
(`jobs/prove_manuscript_source.py`):

| Manuscript metric (published H / A / p) | from CLEAN column | reproduced | from CORRUPT column | gives |
|---|---|---|---|---|
| ipSAE (0.717 / 0.605 / **0.00085**) | `submitted_ipsae` | **0.717 / 0.605 / 0.00085** ✅ | `b2_ipsae_d0chn_max` | 0.699 / 0.702 / **0.76** (null) |
| ipTM (0.858 / 0.832 / **0.0014**) | `boltz2_iptm` | **0.858 / 0.832 / 0.0014** ✅ | `b2_iptm` | 0.940 / 0.941 / **0.80** (null) |
| pLDDT binder (0.900 / 0.873 / **0.00016**) | `boltz2_plddt` | **0.900 / 0.873 / 0.00016** ✅ | — | — |
| pLDDT complex (0.910 / 0.884 / **0.00080**) | `boltz2_complex_plddt` | **0.910 / 0.884 / 0.00080** ✅ | — | — |

Had the table been built on the corrupt `b2_*`, its central claim — *humans scored significantly
higher on every in-silico metric* — would have **collapsed to non-significance** (p ≈ 0.76–0.80). It
did not, which proves the table used the clean competition/legacy columns. Same for Figure 7: the
metric ROC-AUCs are real from the clean columns (ipSAE 0.65, ipTM 0.60, pLDDT 0.64) and **random from
the corrupt ones** (0.54) — so the curves cannot have come from `b2_*`. **No published manuscript
result depends on the corrupt structures.**

### Sensitivity: are the published values affected through any *indirect* channel? No.

Beyond direct use, the corrupt structures could in principle touch a published value only through
shared provenance or through the design selection. Both are ruled out (`jobs/prove_manuscript_source.py`
+ the checks below):

- **Provenance independence.** The clean published columns correlate with the corrupt `b2_*` at
  **ρ ≈ 0.00–0.08** — independent computations (competition/muni pipeline vs. the later Modal rerun).
  The corrupt rerun did not feed the published columns.
- **Selection independence.** The top-100 wet-lab set was chosen by ipSAE: `selection_rank` vs.
  `submitted_ipsae` (clean) is **ρ = −1.000** (a perfect rank match), vs. corrupt `b2_ipsae` ρ = −0.04.
  Mean clean ipSAE selected 0.809 vs. not-selected 0.287 (large separation); mean corrupt b2_ipsae
  0.698 vs. 0.706 (none). The corrupt data had **zero** influence on which designs were tested.
- **Experimental results are structure-independent.** Expression, hit rate, and Kd come from BLI, so
  no predicted structure — corrupt or not — can affect them.

The corruption was therefore fully quarantined to the `b2_*`/`af2m_*` rerun columns and their files —
and both the files (§2) and the columns (§3.5) have now been corrected in `data/`.

**What a co-author should still double-check** (small, precautionary):
- **`data/grand_metrics.csv`**, *if it is regenerated*: it is a derived wide table (not committed) built
  by `scripts/data/build_grand_metrics.py` from the per-model JSONs. Its `b2_*`/`af2m_*`/consensus
  columns will pick up the correction only when rebuilt from corrected per-model metric JSONs; the
  committed source of truth `data/designs.{csv,parquet}` is already fixed.
- **Report the mapping bug upstream** to the muni/Adaptyv Modal re-fold step so the binder→design
  mapping is fixed at the source, not just in this copy.

**Net effect on conclusions:** none. The competition ranking, the hit-rate and affinity results, and
the structural/epitope conclusions all stand; the regeneration additionally converts a liability
(shipping wrong structures) into a strengthened, three-predictor epitope result.

---

## 5. Files

- Corruption evidence: `jobs/diagnose_scalars.py`, `results/structure_integrity_audit.csv`,
  `results/data_integrity_findings.md`
- Regeneration: `infra/runpod_setup.sh`, `infra/prod_boltz.py`, `infra/target.a3m`,
  `infra/af2m_input.csv`; outputs in `data/structures_regen/`
- Validation: `results/boltz2_regen_validation.md`, `results/boltz2_regen_scalars.csv`,
  `results/epitope_boltz2_regen_freq.csv`, `jobs/regen_boltz2_footprints.py`
- Column recompute: `jobs/recompute_scalars.py` (patches `data/designs.{csv,parquet}`),
  `results/recomputed_scalars.csv`
