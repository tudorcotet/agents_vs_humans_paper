"""EXTRA 1 — kinetic decomposition, statistics stage (runs the pre-registered tests).

Spearman rho of each descriptor vs log_kon / log_koff / log_kd, with BCa bootstrap 95% CI,
permutation p, and BH-FDR within each (family x outcome). Tests H-koff / H-kon / H-decomposition.

Run: uv run --with scipy python analyses/structure_epitope/jobs/kinetics_stats.py
Outputs: results/kinetics/kinetics_correlations.csv + results/kinetics/kinetics_findings.md
"""
from __future__ import annotations
import os
import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
OUT = os.path.join(ROOT, "analyses", "structure_epitope", "results", "kinetics")
N_BOOT, N_PERM = 5000, 10000
rng = np.random.default_rng(20260726)

FAMILIES = {
    "stability": ["bsa", "n_iface_binder", "hbonds", "hbond_density_per100",
                  "iface_frac_hydrophobic", "rosetta_dg", "prodigy_protenix_dg"],
    "electrostatic": ["charge_complementarity", "abs_binder_iface_charge",
                      "salt_bridges", "iface_frac_charged"],
    "confidence": ["submitted_ipsae", "boltz2_iptm", "boltz2_plddt", "boltz2_complex_plddt"],
}
OUTCOMES = ["log_kon", "log_koff", "log_kd"]


def bh_fdr(p):
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        k = n - rank + 1
        prev = min(prev, p[i] * n / k)
        q[i] = prev
    return q


def corr_ci_p(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) < 8 or np.std(x) == 0 or np.std(y) == 0:
        return np.nan, np.nan, np.nan, np.nan, len(x)
    rho = stats.spearmanr(x, y).statistic
    def stat(a, b):
        return stats.spearmanr(a, b).statistic
    try:
        bs = stats.bootstrap((x, y), stat, paired=True, vectorized=False,
                             n_resamples=N_BOOT, method="BCa", random_state=rng)
        lo, hi = bs.confidence_interval.low, bs.confidence_interval.high
    except Exception:
        lo = hi = np.nan
    pt = stats.permutation_test((x, y), stat, permutation_type="pairings",
                                n_resamples=N_PERM, random_state=rng)
    return rho, lo, hi, pt.pvalue, len(x)


def run(df, tag):
    rows = []
    for fam, descs in FAMILIES.items():
        for out in OUTCOMES:
            fam_rows = []
            for d in descs:
                if d not in df.columns:
                    continue
                rho, lo, hi, p, n = corr_ci_p(df[d].to_numpy(float), df[out].to_numpy(float))
                fam_rows.append(dict(family=fam, descriptor=d, outcome=out, n=n,
                                     rho=rho, ci_low=lo, ci_high=hi, p_perm=p))
            ps = [r["p_perm"] for r in fam_rows if np.isfinite(r["p_perm"])]
            qs = bh_fdr(ps) if ps else []
            it = iter(qs)
            for r in fam_rows:
                r["q_bh"] = next(it) if np.isfinite(r["p_perm"]) else np.nan
                r["set"] = tag
            rows += fam_rows
    return pd.DataFrame(rows)


def main():
    df = pd.read_csv(os.path.join(OUT, "kinetics_perdesign.csv"))
    df["abs_binder_iface_charge"] = df["binder_iface_charge"].abs()
    primary = df[df.rate_eligible].copy()
    full = df.copy()

    res_primary = run(primary, "primary_n%d" % len(primary))
    res_full = run(full, "sensitivity_n%d" % len(full))
    res = pd.concat([res_primary, res_full], ignore_index=True)
    res.to_csv(os.path.join(OUT, "kinetics_correlations.csv"), index=False)

    p = res_primary.set_index(["descriptor", "outcome"])

    def g(desc, out, col):
        try:
            return p.loc[(desc, out), col]
        except KeyError:
            return np.nan

    L = []
    L.append(f"# EXTRA 1 — kinetic decomposition findings (primary n={len(primary)})\n")
    L.append("Spearman rho (BCa 95% CI), permutation p, BH-FDR per family x outcome. "
             "Signs: lower log_koff = slower dissociation (tighter); higher log_kon = faster association; "
             "lower log_kd = tighter.\n")

    # H-koff: stability tracks koff more than kd
    L.append("## H-koff — stability descriptors track k_off\n")
    L.append("| descriptor | rho(log_koff) [CI] q | rho(log_kd) [CI] q | tracks koff more? |")
    L.append("|---|---|---|---|")
    hkoff_hits = 0
    for d in FAMILIES["stability"]:
        rk, lk, hk, qk = g(d, "log_koff", "rho"), g(d, "log_koff", "ci_low"), g(d, "log_koff", "ci_high"), g(d, "log_koff", "q_bh")
        rd, ld, hd, qd = g(d, "log_kd", "rho"), g(d, "log_kd", "ci_low"), g(d, "log_kd", "ci_high"), g(d, "log_kd", "q_bh")
        if not np.isfinite(rk):
            continue
        more = abs(rk) > abs(rd) and qk < 0.10
        hkoff_hits += int(more)
        L.append(f"| {d} | {rk:+.2f} [{lk:+.2f},{hk:+.2f}] q={qk:.3f} | {rd:+.2f} [{ld:+.2f},{hd:+.2f}] q={qd:.3f} | {'yes' if more else 'no'} |")
    L.append(f"\n**H-koff verdict:** {hkoff_hits} stability descriptors track k_off (FDR<0.10) more strongly than K_D.\n")

    # H-kon: electrostatic tracks kon
    L.append("## H-kon — electrostatic descriptors track k_on\n")
    L.append("| descriptor | rho(log_kon) [CI] q | rho(log_koff) q | rho(log_kd) q |")
    L.append("|---|---|---|---|")
    hkon_hits = 0
    for d in FAMILIES["electrostatic"]:
        rn, ln, hn, qn = g(d, "log_kon", "rho"), g(d, "log_kon", "ci_low"), g(d, "log_kon", "ci_high"), g(d, "log_kon", "q_bh")
        if not np.isfinite(rn):
            continue
        hkon_hits += int(qn < 0.10)
        L.append(f"| {d} | {rn:+.2f} [{ln:+.2f},{hn:+.2f}] q={qn:.3f} | {g(d,'log_koff','rho'):+.2f} q={g(d,'log_koff','q_bh'):.3f} | {g(d,'log_kd','rho'):+.2f} q={g(d,'log_kd','q_bh'):.3f} |")
    # is electrostatic->kon stronger than stability->kon?
    el_kon = np.nanmax([abs(g(d, "log_kon", "rho")) for d in FAMILIES["electrostatic"]])
    st_kon = np.nanmax([abs(g(d, "log_kon", "rho")) for d in FAMILIES["stability"]])
    L.append(f"\n**H-kon verdict:** {hkon_hits} electrostatic descriptors track k_on (FDR<0.10). "
             f"max|rho| electrostatic->kon = {el_kon:.2f} vs stability->kon = {st_kon:.2f}.\n")

    # H-decomposition: confidence metric sig for koff not kd
    L.append("## H-decomposition — confidence metrics vs k_off despite K_D null\n")
    L.append("| metric | rho(log_koff) q | rho(log_kd) q | koff-only? |")
    L.append("|---|---|---|---|")
    hdec = 0
    for d in FAMILIES["confidence"]:
        qk, qd = g(d, "log_koff", "q_bh"), g(d, "log_kd", "q_bh")
        koff_only = np.isfinite(qk) and qk < 0.10 and (not np.isfinite(qd) or qd >= 0.10)
        hdec += int(koff_only)
        L.append(f"| {d} | {g(d,'log_koff','rho'):+.2f} q={qk:.3f} | {g(d,'log_kd','rho'):+.2f} q={qd:.3f} | {'YES' if koff_only else 'no'} |")
    L.append(f"\n**H-decomposition verdict:** {hdec} confidence metric(s) significant for k_off but not K_D.\n")

    # cohort rate-space
    L.append("## Rate-space: do humans and agents occupy the same region?\n")
    for out in ["log_kon", "log_koff", "log_kd"]:
        h = primary[primary.cohort == "human"][out].dropna()
        a = primary[primary.cohort == "agent"][out].dropna()
        pmw = stats.mannwhitneyu(h, a).pvalue if len(h) > 2 and len(a) > 2 else np.nan
        L.append(f"- {out}: human median {h.median():+.2f} vs agent {a.median():+.2f}  (MWU p={pmw:.3f})")
    L.append(f"\n(n_human={int((primary.cohort=='human').sum())}, n_agent={int((primary.cohort=='agent').sum())}; "
             "design 79 is a literature copy and is retained here but flagged.)")

    open(os.path.join(OUT, "kinetics_findings.md"), "w").write("\n".join(L) + "\n")
    print("\n".join(L))
    print(f"\nwrote kinetics_correlations.csv ({len(res)} rows) + kinetics_findings.md")


if __name__ == "__main__":
    main()
