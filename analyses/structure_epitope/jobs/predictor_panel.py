"""SI cross-predictor panel + Protenix-primary interface-size statistics.

Amy's decisions:
  (1) Report interface descriptors on PROTENIX (the section's primary structural model), not the
      Protenix+Chai average, in the main text.
  (2) Provide the full across-predictor comparison in the SI: per-predictor discrimination, pairwise
      agreement (Spearman) and the signed difference.
  (3) Treat the team-contributed ESMFold2 complex panel (100 screened designs) as one more predictor.

Reads results/interface_descriptors_perpredictor.csv (protenix, chai, esmfold2) + the canonical table.
Writes results/predictor_panel.md, predictor_panel_perpred.csv, predictor_panel_pairs.csv.

Run: uv run --with pandas --with scipy --with scikit-learn --with numpy \
       python analyses/structure_epitope/jobs/predictor_panel.py
"""
from __future__ import annotations
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.linear_model import LogisticRegression

ROOT = Path("/Users/amyhe/Desktop/trem2_2026/agents_vs_humans_paper")
RES = ROOT / "analyses/structure_epitope/results"
DESC = ["bsa", "n_iface_binder", "n_iface_target", "hbonds", "salt_bridges",
        "iface_frac_hydrophobic", "iface_frac_polar", "iface_frac_charged", "hbond_density_per100"]
PREDS = ["protenix", "chai", "esmfold2"]


def cliffs(a, b):
    a, b = np.asarray(a), np.asarray(b)
    gt = sum((x > b).sum() for x in a); lt = sum((x < b).sum() for x in a)
    return (gt - lt) / (len(a) * len(b))


def discrim(sub, col):
    """binder vs non discrimination for one descriptor on an expressed subframe."""
    b = sub[sub.is_hit == True][col].dropna(); n = sub[sub.is_hit == False][col].dropna()  # noqa: E712
    if len(b) < 3 or len(n) < 3:
        return None
    U, p = stats.mannwhitneyu(b, n, alternative="two-sided")
    s = sub.dropna(subset=[col]); y = (s.is_hit == True).astype(int)  # noqa: E712
    auc = max(roc_auc_score(y, s[col]), roc_auc_score(y, -s[col]))
    ap = max(average_precision_score(y, s[col]), average_precision_score(y, -s[col]))
    return dict(binder_med=round(b.median(), 1), non_med=round(n.median(), 1),
                mwu_p=p, cliffs=round(cliffs(b.values, n.values), 3),
                auroc=round(auc, 3), ap=round(ap, 3), n_b=len(b), n_n=len(n))


def partial_spearman(x, y, z):
    """Spearman partial correlation of x,y controlling for z (rank-residualize)."""
    rx = pd.Series(x).rank(); ry = pd.Series(y).rank(); rz = pd.Series(z).rank()
    def resid(a, b):
        b = np.c_[np.ones(len(b)), b]
        beta, *_ = np.linalg.lstsq(b, a, rcond=None)
        return a - b @ beta
    ex = resid(rx.values, rz.values); ey = resid(ry.values, rz.values)
    r, p = stats.spearmanr(ex, ey)
    return r, p


def van_elteren(sub, col, stratum):
    """Van Elteren stratified Wilcoxon (design weights 1/(N_h+1)); groups = is_hit."""
    T = 0.0; V = 0.0
    for _, g in sub.dropna(subset=[col]).groupby(stratum):
        y = (g.is_hit == True).astype(int).values  # noqa: E712
        m = int(y.sum()); N = len(y)
        if m == 0 or m == N or N < 2:
            continue
        r = stats.rankdata(g[col].values)
        R = r[y == 1].sum()
        ER = m * (N + 1) / 2.0
        # tie-corrected variance of rank-sum
        _, cnt = np.unique(g[col].values, return_counts=True)
        tie = (cnt**3 - cnt).sum()
        VarR = m * (N - m) / 12.0 * ((N + 1) - tie / (N * (N - 1)))
        w = 1.0 / (N + 1)
        T += w * (R - ER); V += w**2 * VarR
    if V <= 0:
        return np.nan, np.nan
    z = T / np.sqrt(V)
    return z, 2 * stats.norm.sf(abs(z))


def lr_test(sub, col, ctrl="sequence_length"):
    """Likelihood-ratio test that col adds to a length-only logistic model for is_hit."""
    s = sub.dropna(subset=[col, ctrl]).copy()
    y = (s.is_hit == True).astype(int).values  # noqa: E712
    def ll(X):
        Xs = (X - X.mean(0)) / (X.std(0) + 1e-9)
        m = LogisticRegression(max_iter=1000, C=1e6).fit(Xs, y)
        p = np.clip(m.predict_proba(Xs)[:, 1], 1e-9, 1 - 1e-9)
        return (y * np.log(p) + (1 - y) * np.log(1 - p)).sum()
    ll_full = ll(s[[ctrl, col]].values); ll_red = ll(s[[ctrl]].values)
    chi = 2 * (ll_full - ll_red)
    return chi, stats.chi2.sf(chi, 1)


def main():
    lp = pd.read_csv(RES / "interface_descriptors_perpredictor.csv")
    meta = pd.read_parquet(ROOT / "analyses/structure_epitope/data/designs_canonical.parquet")
    keep = ["design_id", "is_hit", "p_expressed", "cohort", "sequence_length"]
    keep = [c for c in keep if c in meta.columns]
    meta = meta[keep].copy()
    if "cohort" in meta:
        meta["is_human"] = (meta["cohort"].astype(str).str.lower().str.contains("human")).astype(int)
    lp = lp.merge(meta, on="design_id", how="left")
    exp = lp[lp.p_expressed == True].copy()

    L = ["# Cross-predictor interface panel (SI) + Protenix-primary statistics\n",
         "Predictors: **Protenix** (primary, 141 designs), **Chai** (141), **ESMFold2** "
         "(team-contributed complex panel, 100 screened designs). Descriptors identical to "
         "`jobs/interface_descriptors.py`; ESMFold2 read the same way (chain A = 175-aa target).\n"]

    # ---- coverage / ESMFold2 quality ----
    L.append("## Coverage & ESMFold2 interface plausibility\n")
    for pred in PREDS:
        sp = lp[lp.predictor == pred]
        spe = exp[exp.predictor == pred]
        zero = (sp["n_iface_binder"] == 0).sum()
        L.append(f"- **{pred}**: {len(sp)} designs ({len(spe)} expressed); "
                 f"BSA median {sp['bsa'].median():.0f} Å²; interface-residue median "
                 f"{sp['n_iface_binder'].median():.0f}; designs with zero interface contacts: {zero}.")
    L.append("")

    # ---- per-predictor discrimination (binder vs non), primary = Protenix ----
    L.append("## Binder-vs-non-binder discrimination, per predictor (expressed designs)\n")
    L.append("| descriptor | predictor | binder med | non med | Cliff δ | MWU p | AUROC | AP | n(b/n) |")
    L.append("|---|---|---|---|---|---|---|---|---|")
    rows = []
    for col in DESC:
        for pred in PREDS:
            spe = exp[exp.predictor == pred]
            d = discrim(spe, col)
            if not d:
                continue
            rows.append(dict(descriptor=col, predictor=pred, **d))
            star = " ⟵ primary" if pred == "protenix" and col == "bsa" else ""
            L.append(f"| {col} | {pred} | {d['binder_med']} | {d['non_med']} | {d['cliffs']} | "
                     f"{d['mwu_p']:.2g} | {d['auroc']} | {d['ap']} | {d['n_b']}/{d['n_n']} |{star}")
    pd.DataFrame(rows).to_csv(RES / "predictor_panel_perpred.csv", index=False)
    L.append("")

    # ---- pairwise agreement + signed difference (all designs with both predictors) ----
    L.append("## Inter-predictor agreement & signed difference (BSA)\n")
    L.append("| pair | n | Spearman ρ | median signed Δ (A−B) | A>B % | Wilcoxon p |")
    L.append("|---|---|---|---|---|---|")
    prows = []
    wide = lp.pivot_table(index="design_id", columns="predictor", values="bsa")
    # All predictors are compared against PROTENIX (the reference); no Chai-vs-ESMFold2 pair (Amy's call).
    for a, b in [("protenix", "chai"), ("protenix", "esmfold2")]:
        w = wide[[a, b]].dropna()
        if len(w) < 5:
            continue
        rho = stats.spearmanr(w[a], w[b])[0]
        d = w[a] - w[b]
        wilc = stats.wilcoxon(d).pvalue
        prows.append(dict(pair=f"{a}_vs_{b}", n=len(w), spearman=round(rho, 3),
                          median_signed_diff=round(d.median(), 1), a_gt_b_frac=round((d > 0).mean(), 3),
                          wilcoxon_p=wilc))
        L.append(f"| {a} − {b} | {len(w)} | {rho:.2f} | {d.median():+.0f} | {100*(d>0).mean():.0f}% | {wilc:.2g} |")
    pd.DataFrame(prows).to_csv(RES / "predictor_panel_pairs.csv", index=False)
    L.append("")

    # ---- Protenix-primary headline + length controls (validate vs averaged) ----
    L.append("## Protenix-primary interface size vs the published Protenix+Chai average\n")
    # Protenix+Chai average is computed here from the per-predictor file — interface_descriptors.csv is
    # now Protenix-only (the section-wide primary), so it can no longer stand in for the average.
    avg = (lp[lp.predictor.isin(["protenix", "chai"])].groupby("design_id")["bsa"].mean()
           .reset_index().rename(columns={"bsa": "bsa_avg"}))
    base = exp[exp.predictor == "protenix"].merge(avg, on="design_id", how="left")
    L.append("| set | binder med | non med | Cliff δ | MWU p | AUROC | AP |")
    L.append("|---|---|---|---|---|---|---|")
    for label, col in [("Protenix+Chai avg (published)", "bsa_avg"), ("Protenix only (new primary)", "bsa")]:
        d = discrim(base, col)
        L.append(f"| {label} | {d['binder_med']} | {d['non_med']} | {d['cliffs']} | "
                 f"{d['mwu_p']:.2g} | {d['auroc']} | {d['ap']} |")
    L.append("")
    L.append("### Length controls — 'size, not length' (reconstructed; validate avg vs published 0.33/0.002, VE 0.0015, LR 0.006)\n")
    L.append("| set | partial ρ(BSA,hit|len) | p | van Elteren p (cohort) | van Elteren p (len-tertile) | LR test p |")
    L.append("|---|---|---|---|---|---|")
    base["len_tertile"] = pd.qcut(base["sequence_length"], 3, labels=False, duplicates="drop")
    for label, col in [("Protenix+Chai avg", "bsa_avg"), ("Protenix only", "bsa")]:
        s = base.dropna(subset=[col, "sequence_length"])
        pr, pp = partial_spearman(s[col].values, (s.is_hit == True).astype(int).values,  # noqa: E712
                                  s["sequence_length"].values)
        _, ve_c = van_elteren(s, col, "is_human") if "is_human" in s else (np.nan, np.nan)
        _, ve_l = van_elteren(s, col, "len_tertile")
        _, lr = lr_test(s, col)
        L.append(f"| {label} | {pr:.3f} | {pp:.3f} | {ve_c:.4f} | {ve_l:.4f} | {lr:.4f} |")

    (RES / "predictor_panel.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))
    print("\nwrote predictor_panel.md + predictor_panel_perpred.csv + predictor_panel_pairs.csv")


if __name__ == "__main__":
    main()
