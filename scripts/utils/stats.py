"""Statistics helpers shared by every analysis.

Use the routines in here instead of re-implementing in your own script.
If you genuinely need a new helper, add it here and let other analyses
benefit too.
"""
from __future__ import annotations

import warnings

import numpy as np
from scipy import stats


def boot_diff_ci(
    a: np.ndarray,
    b: np.ndarray,
    statistic,
    n: int = 1000,
    seed: int = 0,
) -> tuple[float, float, float]:
    """Bootstrap 95% CI on `statistic(a) - statistic(b)`.

    Returns (lo, hi, mean) — the percentile-method 2.5/97.5 bounds and the
    bootstrap mean of the difference.
    """
    rng = np.random.default_rng(seed)
    a, b = np.asarray(a), np.asarray(b)
    diffs = np.empty(n)
    for i in range(n):
        sa = rng.choice(a, size=len(a), replace=True)
        sb = rng.choice(b, size=len(b), replace=True)
        diffs[i] = statistic(sa) - statistic(sb)
    lo, hi = np.percentile(diffs, [2.5, 97.5])
    return float(lo), float(hi), float(diffs.mean())


def boot_spearman_ci(
    x: np.ndarray,
    y: np.ndarray,
    n: int = 1000,
    seed: int = 0,
) -> tuple[float, float, float]:
    """Bootstrap 95% CI on Spearman ρ.

    Returns (rho_observed, lo, hi).
    """
    x, y = np.asarray(x), np.asarray(y)
    rng = np.random.default_rng(seed)
    rhos = np.empty(n)
    idx = np.arange(len(x))
    for i in range(n):
        s = rng.choice(idx, size=len(idx), replace=True)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            rhos[i] = stats.spearmanr(x[s], y[s]).correlation
    rhos = rhos[~np.isnan(rhos)]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rho_obs = float(stats.spearmanr(x, y).correlation)
    lo, hi = np.percentile(rhos, [2.5, 97.5])
    return rho_obs, float(lo), float(hi)


def fisher_2x2(table: np.ndarray, *, alternative: str = "two-sided") -> dict:
    """Fisher's exact on a 2×2 table. Returns p, odds_ratio, and the table."""
    odds, p = stats.fisher_exact(table, alternative=alternative)
    return {
        "table": [[int(table[0, 0]), int(table[0, 1])],
                   [int(table[1, 0]), int(table[1, 1])]],
        "odds_ratio": float(odds),
        "p": float(p),
        "alternative": alternative,
    }


def mann_whitney(a: np.ndarray, b: np.ndarray, *, alternative: str = "two-sided") -> dict:
    """Mann-Whitney U with the median difference and both group sizes."""
    a, b = np.asarray(a), np.asarray(b)
    res = stats.mannwhitneyu(a, b, alternative=alternative)
    return {
        "U": float(res.statistic),
        "p": float(res.pvalue),
        "n_a": int(len(a)),
        "n_b": int(len(b)),
        "median_a": float(np.median(a)) if len(a) else float("nan"),
        "median_b": float(np.median(b)) if len(b) else float("nan"),
        "median_diff": float(np.median(a) - np.median(b)) if len(a) and len(b) else float("nan"),
        "alternative": alternative,
    }


def bh_fdr(pvalues: np.ndarray, q: float = 0.10) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values at FDR `q`.

    Returns a boolean mask of which tests survive the FDR threshold.
    """
    p = np.asarray(pvalues, dtype=float)
    order = np.argsort(p)
    ranks = np.empty_like(order)
    ranks[order] = np.arange(1, len(p) + 1)
    crit = ranks / len(p) * q
    sig_ordered = p[order] <= crit[order]
    if not sig_ordered.any():
        return np.zeros_like(p, dtype=bool)
    last = np.where(sig_ordered)[0].max()
    keep = np.zeros_like(p, dtype=bool)
    keep[order[: last + 1]] = True
    return keep
