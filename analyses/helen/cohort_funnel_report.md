# The selection funnel by cohort: did the in-silico edge convert?

**Headline:** humans submitted markedly higher-ipSAE designs (median 0.8319 vs 0.6859, Mann-Whitney p=0.00085, rank-biserial -0.34), so a larger share survived the top-100 ipSAE cut (80% of human vs 58% of agent submissions, Fisher p=0.0053). Yet among lab-tested designs the hit rate is statistically tied (38% human vs 34% agent, OR 0.8348, Fisher p=0.83). The proxy advantage does not carry through to the wet lab.

## Method

`load_designs()` (141). Funnel per cohort: submitted -> survived the top-100 ipSAE cut (`submitted_to_lab`) -> expressed (`expressed`, conditional on screened) -> bound (`is_hit`, conditional on screened), plus end-to-end bound/submitted. Each stage: 2x2 Fisher exact two-sided with raw counts and the agent-vs-human odds ratio. Mechanism: Mann-Whitney on `submitted_ipsae` over all submissions + rank-biserial effect size. BH-FDR q=0.10 over the four stage tests.

## Results

- **Survived ipSAE cut:** human 65/81 (80%) vs agent 35/60 (58%); OR 0.3446, Fisher p=0.0053, BH-sig True.
- **Expressed | screened:** human 60/65 (92%) vs agent 29/35 (83%); Fisher p=0.19, BH-sig False.
- **Bound | screened:** human 25/65 (38%) vs agent 12/35 (34%); OR 0.8348, Fisher p=0.83, BH-sig False.
- **Bound | submitted (end-to-end):** human 25/81 (31%) vs agent 12/60 (20%); OR 0.56, Fisher p=0.18, BH-sig False.
- **ipSAE (all submitted):** human median 0.8319 (n=80) vs agent 0.6859 (n=56); Mann-Whitney p=0.00085, rank-biserial -0.34.

## Caveats

- **Guaranteed-inclusion confound.** Per the Notion draft, top agent designs had a guaranteed path into the screened set, so the agent 'survival' rate is not a pure ipSAE cut — it understates the true ipSAE gap between cohorts. Interpret the survival stage carefully.
- **ipSAE is not method-neutral.** Humans used methods (BindCraft, AF-hallucination, Mosaic) that tend to yield higher Boltz-2 interface confidence; agents leaned on PXDesign. 'Humans had higher ipSAE' is partly a method-mix effect, not proof of intent.
- **Small agent N** (60 submitted, 35 screened, 12 binders): the tied hit rate is consistent with no difference *and* with a modest difference the study cannot resolve. Effect sizes + counts reported.
- **Top-100 cherry-pick.** All conditional rates are downstream of the ipSAE triage; the 41 unscreened designs are never lab-tested.
- Conditional definitions are explicit (| screened); do not read 'bound | screened' as an unconditional binding probability.
