# <Replace with your analysis title>

**Headline:** 37% overall hit rate (human 38%, agent 34%).

## Method

Briefly: what was computed, which filter was applied, which stat test was used.
Cite the column names from `data/designs.csv` so the reader can re-derive.

## Results

```
{
  "n_screened": 100,
  "n_human": 65,
  "n_agent": 35,
  "hit_rate_overall": 0.37,
  "hit_rate_human": 0.38461538461538464,
  "hit_rate_agent": 0.34285714285714286
}
```

## Caveats

- Small N: cohort hit rates are computed on the 100 designs sent to the BLI assay.
- The top-100 ipSAE cherry-pick is an in-silico selection; the 41 designs not screened
  are absent from this analysis.
- This template's numbers are illustrative — replace with the test that actually answers
  your question.
