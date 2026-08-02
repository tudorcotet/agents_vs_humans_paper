# Structure and Epitope Analysis (draft v1)

## Summary of what this section establishes

Across 141 designs modeled as complexes with TREM2, this section establishes three core results and two  
extensions:

- **Convergent, non-discriminating epitope.** Binders and non-binders alike engage the same apical ligand-binding face (hydrophobic tip + CDR-like loops); *epitope choice does not separate active from inactive designs* (§1).
- **Interface size is the discriminator.** Binders bury a larger interface (1,823 vs. 1,515 Å², AUROC 0.72, q=0.004; holds across all five predictors); interface chemistry, secondary structure, and contact motifs add nothing beyond size — even interface H-bonds and salt bridges, re-scored with proper geometry (Rosetta HBondSet) on clash-verified relaxed structures, reduce to interface size (§2). 
- **Affinity — learned scores fail; physical/structural quantities and length succeed.** Neither the folding-model confidence scores (ipTM, pLDDT; ipSAE weakly, ρ=0.40) nor the language-model likelihoods (ESM-2, SaProt) track K_D; interface size (ρ=0.47), buried area (0.50), PRODIGY ΔG (0.52) and binder length (0.65) do (pooled n=36). The human-vs-agent split is mostly statistical power, not real: interface size and buried area are underpowered — not absent — in the 12 agents, and length works in both cohorts; only PRODIGY (human-only) and ipSAE (an agent-only inversion) are genuinely cohort-specific (§3).
- **Kinetics (§3).** K_D variation is dominated by k_off (k_off–K_D ρ=0.92); k_on is not predicted by any descriptor. Humans and agents show no detectable difference in rate space (equivalence for k_on; k_off excludes only a ≥~7-fold gap).
- **Affinity in context (§3).** Best measured affinity 1.11 nM (SPR) — ~72× tighter than the best prior de novo TREM2 binder (Perera et al., 2026; 80 nM by MST) on the identical construct.

---

## 1. Epitope landscape of the TREM2 IgSF domain

In the literature, scFv-2 and scFv-4 bind a *distal* TREM2 epitope, β-strands A/F/G plus the C–C′ loop, as supported directly by the respective co-crystal structures (Szykowska et al., 2021; PDB 6YYE, 6Y6C). Across several independent structure-prediction methods, we found that **our designs do not target these scFv antibody epitopes** (binder overlap ≈0, p≈0.8–0.9) but **do overlap the phosphatidylserine ligand-binding site (Sudom et al., 2018; PDB 6B8O) which is on the face opposite the antibody epitopes**. Superposing each co-crystal onto a common TREM2 frame (rigid-body Cα fit, ≤0.75 Å) makes this geometry explicit: the scFv-2/-4 antibodies engage the distal β-A/F/G face while our tightest binder (design 17) and the phosphatidylserine ligand occupy the opposite apical face (**Fig 8a–d**). 

Both binders and non-binders overwhelmingly engage the TREM2 apex: the hydrophobic tip (W44, L69/W70/L71, F74, L89) and the adjacent CDR-like loops (residues ~67–75); contact frequency in binders is ≥0.97 at L69/L71/L72/F74 (Fig 8e). The apical convergence **holds across three independent structure predictors** (Protenix, Chai-1, Boltz-2). Protenix and Chai concur (ρ=0.88), and Boltz-2 independently reproduces it (ρ=0.78; its top-10 binder-contacted residues are all six hydrophobic-tip positions W44, L69/W70/L71, F74, L89). Epitope binning with Jaccard on consensus footprints within silhouette-selected threshold resolves **one dominant apical epitope bin containing 92/97 screened designs**. The complex gallery in **Fig 10** collects these structures on the shared TREM2 frame — the top designs in leaderboard order (a), the best design's interface with its EvoEF2 hotspots (b), and all 37 hit binders superposed, converging on the apical face (c). This **quantitatively confirms and extends Perera et al.'s (2026) qualitative "convergence on the hydrophobic tip"** on ~7× the designs. However, because tip contact is near-universal (36/37 binders and 47/52 non-binders), it does not separate binders from non-binders. 

An EvoEF2 alanine scan of the ten top binders (437 interface residues) independently confirms the energetic hotspots within the epitope. The TREM2 residues whose alanine substitution most destabilizes binding, averaged across designs, are **dominated by the hydrophobic tip**: L71, L72, L69, L75, L89, F74, W44, W70. The
basic-patch residue **R76** also scores highly on the unrelaxed structures, but it is the one hotspot that does not survive backbone relaxation, likely due to the caveats that EvoEF2 ΔΔG in the ML-predicted complexes can be inflated by backbone clashes that relaxation removes. Therefore, we nominate the tip residues **W44, L69, L71, F74, L75, L89 as the highest-value targets for follow-up mutagenesis** with R76 a tentative secondary for its functionally implicated role in ligand-binding and multimerization (Kober et al., 2016), although there is limited evidence that naturally occurring variation at R76 has a recognized cause of human disease in the same sense as R47H or R62H. 

Additionally, the competition explicitly forbade targeting the stalk (residues 133–174), which is present in the folded construct. Only **8/100 screened designs contact any stalk residue and just one is stalk-dominant** (design 85, an agent), most "contacts" being footprints spilling past the IgSF C-terminus. The human-vs-agent difference (10.8% vs. 2.9%, Fisher p=0.26) is insignificant. 

## 2. Interface features of active vs. inactive designs

Within expressed designs (37 binders and 52 non-binders), we compared a panel of interface descriptors and report the findings based on Protenix, the primary structural model used in this section. 

First of all, we found that the **interface size is the discriminator.** Buried surface area is larger in binders (median 1,823 vs. 1,515 Å²; Mann–Whitney p=4.2×10⁻⁴, BH q=0.004; univariate AUROC 0.72, AP 0.62), as are the binder- and target-side interface residue counts (q=0.009, 0.067). The effect is **predictor-independent**: it reproduces on Chai (AUROC 0.68, p=5×10⁻³) and on the architecturally independent, single-sequence ESMFold2 complex panel (AUROC 0.69, p=2.6×10⁻³), which recovers it with no MSA (SI cross-predictor table). The size effect is not a length artifact because on Protenix it survives controlling for binder length (partial ρ=0.33, p=0.001).

However, among the 90 expressed designs, **we could not resolve an interface-chemistry contribution beyond size.** Neither overall interface composition (the hydrophobic or charged fraction) nor the two specific polar-contact metrics (hydrogen bonds and sidechain salt bridges) clearly distinguished binders from non-binders independently of buried surface area (BSA). We defined hydrogen bonds using Rosetta’s orientation- and energy-dependent **HBondSet** and salt bridges as contacts between oppositely charged sidechain groups (Asp/Glu ↔ Lys/Arg ≤4 Å). Both metrics were evaluated on clash-verified, relaxed structures. Although hydrogen-bond count showed modest discriminatory value in the as-predicted models (AUROC = 0.68, p = 0.003), this association did not persist after structural relaxation. In the clash-free structures, hydrogen-bond count was strongly correlated with BSA (ρ = 0.78), and its association with binding after controlling for BSA was small and nonsignificant (partial ρ = −0.12). Sidechain salt-bridge counts were likewise correlated with interface size and showed no significant size-independent association with binding. Thus, interface size is the only interface-level discriminator supported by this dataset. **This result does not imply that interface chemistry is biologically unimportant**; rather, the limited sample size and chemical diversity of the expressed designs may prevent detection of smaller chemistry-dependent effects.

A secondary-structure and contact-motif analysis reinforces this. Binder and non-binder interfaces are equally helical (interface helix fraction 0.90 vs. 0.91) whereas β-strand augmentation of the TREM2 sheet is rare (31 designs) and not a dominant binder feature (only observed in 2 binders). Among the typed contact-pair enrichment, including the aromatic-to-apex contacts (OR 3.1, p=0.025) and a cation–π-like Arg/Lys-to-tip tendency (OR 2.3, p=0.085), none survive FDR. The convergent-motif search across independent teams recovers no specific residue-pair motif beyond the universal apical engagement itself. The same holds at the level of the whole binder fold: **folds are overwhelmingly α-helical** (median helix fraction 0.82; 105/141 predominantly α), splitting by cohort — **agents concentrated on helical bundles** (32/35 screened) while **humans explored the few β and mixed folds** (9 β/mixed), the structural counterpart of the sequence-level tool monoculture (all six agents used PXDesign) noted in the Abstract; but fold, like epitope, is ~95% one category (138/141 helical-dominated, 92/97 apical), so it too does not discriminate binders. **Together, the interpretation is consistent across every axis: active designs differ from inactive ones in how much interface they bury, not in where they bind, in a specific chemical motif, or in fold**, an independent replication of Cotet et al. (2025; interface area their strongest structural feature) and Perera et al. (2026; small interface → failure to bind).

Furthermore, we classiedfied design outcomes using the Failure taxonomy (Bennett et al., 2023), yeilding a a three-way decomposition: **11 Type 0 expression failures, 17 Type I folding failures, 35 Type II interface failures**, and 37 active binders. Folding failure was defined as a binder Cα-RMSD >3 Å between the free-monomer ESMFold prediction and the binder conformation in the predicted complex; interface failures retained the predicted fold but did not bind.Among the 52 expressed non-binders, interface failures outnumbered folding failures by approximately 2:1. Active binders also showed greater monomer–complex agreement than non-binders (median Cα-RMSD 0.63 vs. 1.27 Å). Together with the interface-size analysis, these results suggest that **insufficient interface formation was a more common limitation than gross folding failure**.

## 3. Affinity: computational predictors and literature context

For affinity (P_kd, n=36), we find that **learned scores do not rank affinity**. Folding-model confidence (ipTM ρ=0.28, pLDDT 0.14) and language-model likelihoods (ESM-2 0.09, SaProt 0.06) have CIs spanning zero, with the lone exception of **ipSAE** (ρ=0.40, q=0.086), whose signal comes *entirely from the agent cohort* (0.66 vs. human 0.30, ns). **Physical/structural quantities, on the other hand, do correlate with experimental affinity measurements**: interface size (0.47), buried area (0.50), and PRODIGY ΔG (0.52) clear FDR, and **binder length is the strongest correlate** (0.65), though interface burial keeps its signal after controlling for length (partial ρ≈0.51), so it is not merely a length proxy. Because interface area and PRODIGY are read off the predicted complex, this is a structural property of the model tracking affinity, not sequence-level K_D prediction.

Regarding the human-vs-agent comparison, we concluded that the contrast is a **decomposition, not a within-human story**. For interface size and buried area the 12 agents match the humans' direction but miss significance (ρ 0.36–0.41, n.s.) purely from n=12 underpower, whereas **PRODIGY is genuinely human-only** (agent ρ=0.04). Finally, **ipSAE is a genuine agent-only inversion**, the only metric that tracks affinity in agents but not humans (SI Tables S8–S9, Fig S4). In rate space, **humans and agents show no detectable difference** (log k_on/k_off/K_D MWU p=0.84/0.47/0.51; Fig 12a), with formal equivalence established for k_on and only a ≥7-fold cohort difference excludable for k_off and K_D.

Our binder with the best measured affinity, **1.11 nM (design 17, SPR-kinetic K_D)**, is competitive with the strongest reported de novo binders across all targets and **≈72× tighter than the best prior de novo TREM2 binder** (Perera et al., 2026; Odesign2, 80 nM), although this comparison should be interpreted cautiously because it contrasts an SPR kinetic (K_D) with an MST equilibrium measurement having an approximately 5 nM detection floor; four of our 36 quantified binders fall below that range. The target constructs are nevertheless closely matched across the directly relevant datasets, all spanning TREM2 residues 19–174 (Fig S6). 

## Methods

**Populations.** P_all=141, P_screened=100 (submitted to lab), P_expressed=89 (the binder/non-binder
denominator), P_kd=36 (binders with a fitted, uncensored K_D). No interval censoring exists; the 37th
binder (design 5) had no fittable K_D and is excluded from affinity analyses only.

**Structures.** **Protenix-v2** is the primary structural model (Apache-2.0, open-weights); **Chai-1**
(Apache-2.0) is the confirming second predictor. **Boltz-2** and **AF2-Multimer** are secondary — Boltz-2
confirms the epitope map (§1), and both enter the SI cross-predictor comparison (§3, S8–S9) — but are
not used in the primary analysis. A team-contributed **ESMFold2** complex panel (100 screened designs)
adds a fifth, architecturally independent (single-sequence, no-MSA) predictor to that comparison.
The TREM2 target chain is the 175-aa BLI construct (residues 19–174 + linker + His); numbering is mature
UniProt Q9NZC2 (= construct position + 17).

**Descriptors (Tier A, permissive).** SASA and ΔSASA by Bio.PDB Shrake–Rupley (Biopython, BSD) — *not*
FreeSASA; interface detection by Biopython NeighborSearch (4.5 Å heavy-atom + ΔASA>1 Å²); secondary
structure by pydssp (MIT) — *not* mkdssp; typed contacts (H-bond/salt-bridge) by a direct geometric
implementation *(a coarse geometric proxy — any N/O pair ≤3.5 Å, no angle; superseded for the §2 chemistry
conclusion by the tightened Rosetta-HBondSet / sidechain-salt-bridge analysis on clash-verified relaxed
structures, "Relaxation" below)*; contacts-based ΔG by PRODIGY (Apache-2.0); interface ΔG / alanine scan by EvoEF2 (MIT) on
the **unrelaxed** Protenix complexes. All Tier A structural descriptors are
computed on the unrelaxed predicted complexes. **All interface descriptors in the main text (§2, §3,
Figs 9/11) are computed on Protenix, the section's single primary structural model;** the same descriptors
are additionally computed on Chai, ESMFold2, Boltz-2, and AF2M complexes, and written
to a per-predictor file for the SI cross-predictor comparison (SI Tables S7–S9), which anchors every comparison to Protenix. The
interface-size discriminator holds in all five predictors, while they differ systematically in absolute
interface size (Protenix tightest, AF2M nearly identical to it, Chai largest; Protenix–Chai Spearman
ρ=0.32, Protenix–ESMFold2 0.50). AF2M complexes ship as PDB (the rest as mmCIF); `af2m_mean_plddt` and
`px_mean_plddt` are on a 0–100 scale and are rescaled for the pLDDT comparison.

**Epitope footprints.** A TREM2 residue is scored as contacted when a target heavy atom lies ≤4.5 Å from a  
binder heavy atom and buries >1 Å² on complex formation; the **consensus footprint** requires both Protenix  
and Chai to agree (median Jaccard 0.435). The per-residue epitope *map* is concordant across the two  
predictors (Spearman ρ=0.88; Fig S5a) even where single-design footprints vary, so epitope-landscape  
conclusions are predictor-robust, while per-design binder/non-binder enrichment (which needs exact  
agreement) is not (Fig S5); cutoff sensitivity (4.0/4.5/5.0 Å → 10/12/14 consensus residues) is in Fig S5c.

**Statistics.** Binder-vs-non-binder: Mann–Whitney U with BH-FDR and Cliff's δ; AP (primary) and AUROC.
Affinity: Spearman ρ on pKd with BCa bootstrap 95% CI (5,000), permutation p (10,000), BH-FDR per tier,
partial correlation controlling for binder length, per-cohort, and literature-copy sensitivity. Fisher
exact with BH-FDR for per-residue epitope enrichment. All seeds recorded.

**Kinetics + literature (extended analyses).** Per-design k_on/k_off aggregated from the per-replicate
BLI records under the canonical replicate filter; QC gates (k_off solver floor, k_on mass-transport
ceiling, model adequacy) applied; correlations as above on log k_on/k_off/K_D with BH-FDR per
(family × outcome). Literature K_D values are curated from the
open-items register and several await citation verification;
scFv/PS reference epitopes are computed directly from PDB coordinates (6YYE/6Y6C, Szykowska et al., 2021; 6B8O, Sudom et al., 2018), not from the
published lists. IgSF numbering note: the domain spans **UniProt 19–133** (115 residues; the "19–132"
form is off by one, immaterial to the epitope sets as no reference residue lies at 133).

**Rules-compliance (Methods note).** By Smith–Waterman local alignment and longest-common-substring vs.  
scFv-2/-4 and MOR44698 CDRs, 134/141 designs are comfortably compliant; seven human antibody-scaffold  
designs share framework identity, and three (BraiNSEY 77/79/80) are near-copies of a known anti-TREM2  
antibody (independently flagged as a literature copy in our sequence audit). **Design 79 — a 21.66 nM binder at leaderboard rank**  
**#7 — contains the exact scFv-4 HCDR3 and 100% IGKV1-39 light-chain identity; its status as a de novo**  
**human design should be reconsidered**.

