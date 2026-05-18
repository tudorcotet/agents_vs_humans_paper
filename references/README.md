# `references/` — topic-organised BibTeX libraries

These are the source-of-truth bibliographies. The curated subset cited
in the paper lives in `paper/references.bib`.

## Files

- `trem2.bib` — TREM2 biology, clinical references (AL002, VHB937, …).
- `methods.bib` — design methods (RFDiffusion, ProteinMPNN, BoltzGen,
  Boltz-2, AF3, Chai-1, BindCraft, Mosaic, …).
- `agentic_biology.bib` — agent platforms, multi-agent design papers,
  LLM-driven scientific discovery.
- `competitions.bib` — EGFR 2024 (Adaptyv), Nipah 2026, prior wet-lab
  competition writeups.

## Adding a reference

1. Drop the BibTeX entry in the appropriate topic file.
2. If you cite it in the paper, copy the entry into
   `paper/references.bib` and use the same key in the `\cite{...}` call.
3. Prefer DOIs over URLs. Prefer published versions over preprints
   once available.

## Style

- Keys: `lastauthor_keyword_year` (e.g. `dauparas_proteinmpnn_2022`).
- Always include `doi = {...}` when known.
- Cite the version of the tool you used (e.g. Boltz-2 specifically, not
  generic "Boltz").
