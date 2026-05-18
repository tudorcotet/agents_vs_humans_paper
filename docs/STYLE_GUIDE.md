# Writing & visual style guide

The paper has a single voice and a single look. This is how to keep it.

## Voice — the paper

**Clinical, measured, methods-paper voice.** Short sentences. Numbers
where you'd otherwise reach for an adjective. The reader is a competent
scientist who has no patience for flourish.

What "clinical" means in practice:

- **Cite the column name.** "The hit rate (`is_hit` in `data/designs.csv`)
  was 37% (37/100) overall." Not "We found that 37% of designs bound."
- **State the test and the alternative.** "We tested whether agent and
  human hit rates differ using Fisher's exact (two-sided)." Not "We
  found a significant difference."
- **Always pair effect size and p-value.** "Δhit-rate = +7 pp,
  p = 0.42." Not just "p = 0.42."
- **Caveat upfront, not in the discussion.** Selection bias, small N,
  per-team imbalance — say it in the figure caption and the result
  paragraph, not three pages later.

What to avoid (these are AI tells; the `humanizer` skill catches them):

- "Robust", "leverage", "delve", "comprehensive", "novel insights".
- The rule of three: "X, Y, and Z" repeated as a stylistic crutch.
- Inflated symbolism: "shines a light on", "underscores the importance of".
- Em-dash parallelism overload: "It's not just X — it's Y."
- "Stark", "striking", "remarkable" — let the numbers do the work.

## Voice — adjacent writing (blog, social)

The hackathon also has a blog post and a LinkedIn / X thread. They live
outside this repo and follow a different voice — looser, more narrative,
first-person where useful. Don't blend the two voices. The paper is the
paper; the blog is the blog.

## Visual style — at a glance

- **Palette:** `theme/palettes.json`. Two hues for cohort, brand sequential
  for ordered data. See [`FIGURES.md`](FIGURES.md).
- **Fonts:** GT Pressura Extended (display), Geist (body), Geist Mono (data).
- **Frame:** white canvas, 1px borders, no rainbow, no 3D.
- **Sensorgrams:** use real curves where possible; synthesised
  illustrations should match the visual conventions of the published
  blog figures in `figures/blog/`.

## Captions

Caption structure:

> **Panel A.** What you're looking at, one sentence. Per-cohort split
> shown (human, n=N; agent, n=N). Effect size + 95% CI + p-value where
> a statistical test was run. Caveat if the test is underpowered.

Don't put a sentence-long story in a caption. Put it in the body text.

## Numbers

- Affinities: `<KD>` in nM when ≥1 nM, in μM when ≥1 μM, in pM when ≤1 nM.
  Always include `pKD` in any plot used for fitting (Spearman, Mann-Whitney).
- Percentages: rounded to whole numbers if the denominator is < 100, one
  decimal if > 100. Always include `(n/N)` next to a percentage.
- p-values: report exact (e.g. `p = 0.018`). Use `p < 0.001` only when the
  test is genuinely below 10⁻³.

## Plural vs singular

- "Agents" / "humans" when talking about the cohorts at the level of the
  paper. "Agent designs" / "human designs" when talking about the rows.
  Don't say "agent's" — that's possessive, not what you mean.

## Hackathon teams

- 10 **human teams**: NovoFy, MRAZS, EuroBros, DeNovo, StanFold,
  1000Tokens, crow, BART.bio, BraiNSEY, Marcel.
- 6 **agent teams**: claude-sonnet-4.6, qwen-3.5-plus, grok-4.1-fast,
  Gemini 3.1 Pro, GLM 5, GPT 5.2.
- **MRAZS used Mosaic** (the JAX composite-objective wrapper maintained
  by Escalante Bio). MRAZS ≠ Escalante Bio. Don't conflate.
- **`Foundry`** is the method label for "agent picked its own tools." It
  is not a single algorithm. Treat it as the `Other` family.

## Citations

- Use BibTeX. Add entries to the appropriate file in `references/`
  (`trem2.bib`, `methods.bib`, `agentic_biology.bib`, `competitions.bib`).
- Cite tools (Boltz-2, AF3, ProteinMPNN, etc.) at first mention.
- Cite the AL002 and VHB937 patents / trial papers when invoking them
  as positive controls.

## Glossary

A short glossary lives at [`docs/GLOSSARY.md`](GLOSSARY.md). Add any term
that an external reader might trip over.
