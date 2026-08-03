"""EXTRA 1 — kinetic decomposition figures.

Fig 14 (new main-text): (a) rate-rate map — log k_on x log k_off with iso-affinity diagonals,
coloured by cohort; (b) forest plot of Spearman rho for each descriptor against log k_on / log
k_off / log K_D (the visual form of the decomposition argument).
Fig S9: k_off vs the top stability descriptor; k_on vs the top electrostatic descriptor.

Run: uv run --with scipy python analyses/structure_epitope/jobs/kinetics_figure.py
"""
from __future__ import annotations
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
OUT = os.path.join(ROOT, "analyses", "structure_epitope", "results", "kinetics")
FIG = os.path.join(ROOT, "analyses", "structure_epitope", "figures")
AGENT, HUMAN, INK, ACC = "#30C5F5", "#1FE48F", "#0F1419", "#FF4628"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK, "axes.linewidth": 0.6,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 300})

df = pd.read_csv(os.path.join(OUT, "kinetics_perdesign.csv"))
df["cohort"] = df.get("cohort", pd.Series(np.where(df.is_human, "human", "agent")))
prim = df[df.rate_eligible].copy()

# ---------------------------------------------------------------- Fig 14
fig = plt.figure(figsize=(12.5, 5.2))
gs = fig.add_gridspec(1, 2, width_ratios=[1.05, 1.0], wspace=0.28)

# (a) rate-rate map
ax = fig.add_subplot(gs[0, 0])
xr = (prim.log_kon.min() - 0.3, prim.log_kon.max() + 0.3)
for kd_nM, lab in [(1, "1 nM"), (10, "10 nM"), (100, "100 nM"), (1000, "1 µM")]:
    xs = np.array(xr)
    ys = xs + (np.log10(kd_nM) - 9)      # log_koff = log_kon + log10(KD_M)
    ax.plot(xs, ys, ls="--", lw=0.7, color="#C3C9D2", zorder=1)
    ax.text(xr[1], xr[1] + (np.log10(kd_nM) - 9), f" {lab}", fontsize=7, color="#8A93A0",
            va="center", ha="left")
for coh, col in [("human", HUMAN), ("agent", AGENT)]:
    s = prim[prim.cohort == coh]
    ax.scatter(s.log_kon, s.log_koff, s=46, color=col, edgecolor=INK, linewidth=0.4,
               label=f"{coh} (n={len(s)})", zorder=3)
# flag design 79 (literature copy)
d79 = prim[prim.design_id == 79]
if len(d79):
    ax.scatter(d79.log_kon, d79.log_koff, s=150, facecolor="none", edgecolor=ACC,
               linewidth=1.4, zorder=4)
    ax.annotate("79 (lit. copy)", (d79.log_kon.iloc[0], d79.log_koff.iloc[0]),
                fontsize=7, color=ACC, xytext=(5, 6), textcoords="offset points")
ax.set_xlabel("log₁₀ k_on  (M⁻¹s⁻¹)")
ax.set_ylabel("log₁₀ k_off  (s⁻¹)")
ax.set_title("(a) rate–rate map: same affinity, different kinetic routes", fontsize=9, loc="left")
ax.legend(frameon=False, fontsize=8, loc="lower right")

# (b) forest plot of rho across outcomes
ax2 = fig.add_subplot(gs[0, 1])
try:
    corr = pd.read_csv(os.path.join(OUT, "kinetics_correlations.csv"))
    corr = corr[corr.set.str.startswith("primary")]
    show = ["bsa", "hbonds", "iface_frac_hydrophobic", "rosetta_dg",
            "charge_complementarity", "salt_bridges",
            "submitted_ipsae", "boltz2_iptm"]
    show = [d for d in show if d in set(corr.descriptor)]
    ocol = {"log_kon": AGENT, "log_koff": ACC, "log_kd": "#7A5CFF"}
    offs = {"log_kon": +0.22, "log_koff": 0.0, "log_kd": -0.22}
    yl = np.arange(len(show))[::-1]
    for out, col in ocol.items():
        sub = corr[corr.outcome == out].set_index("descriptor")
        ys = [yl[i] + offs[out] for i, d in enumerate(show)]
        xs = [sub.loc[d, "rho"] if d in sub.index else np.nan for d in show]
        lo = [sub.loc[d, "ci_low"] if d in sub.index else np.nan for d in show]
        hi = [sub.loc[d, "ci_high"] if d in sub.index else np.nan for d in show]
        ax2.errorbar(xs, ys, xerr=[np.array(xs) - np.array(lo), np.array(hi) - np.array(xs)],
                     fmt="o", ms=4.5, color=col, ecolor=col, elinewidth=1.1, capsize=2, label=out)
    ax2.axvline(0, color=INK, lw=0.7)
    ax2.set_yticks(yl)
    ax2.set_yticklabels(show, fontsize=8)
    ax2.set_xlabel("Spearman ρ (BCa 95% CI)")
    ax2.set_title("(b) descriptor ↔ rate decomposition", fontsize=9, loc="left")
    ax2.legend(frameon=False, fontsize=7.5, loc="lower left", ncol=1)
except FileNotFoundError:
    ax2.text(0.5, 0.5, "run kinetics_stats.py first", ha="center", transform=ax2.transAxes)

nh = int((prim.cohort == "human").sum()); na = int((prim.cohort == "agent").sum())
fig.suptitle(f"Fig 12. Kinetic decomposition of binding: K_D = k_off / k_on across {len(prim)} rate-eligible designs "
             f"({nh} human / {na} agent, one target, one assay).", fontsize=9.2, x=0.01, ha="left")
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(os.path.join(FIG, "fig12_kinetic_decomposition.png"), bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- Fig S9 secondary scatters
fig2, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 3.9))
def scatter(ax, xcol, ycol, xlabel, ylabel, title):
    for coh, col in [("human", HUMAN), ("agent", AGENT)]:
        s = prim[prim.cohort == coh].dropna(subset=[xcol, ycol])
        ax.scatter(s[xcol], s[ycol], s=40, color=col, edgecolor=INK, linewidth=0.4, label=coh)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel); ax.set_title(title, fontsize=8.5, loc="left")
scatter(a1, "bsa", "log_koff", "interface ΔSASA (Å²)", "log₁₀ k_off",
        "(a) stability → dissociation")
scatter(a2, "charge_complementarity", "log_kon", "interface charge complementarity", "log₁₀ k_on",
        "(b) electrostatics → association")
a1.legend(frameon=False, fontsize=7)
fig2.suptitle("Fig S9. Rate-resolved structure–kinetics: stability vs k_off, electrostatics vs k_on.",
              fontsize=8.6, x=0.01, ha="left")
fig2.tight_layout(rect=[0, 0, 1, 0.94])
None  # dropped former-FigS9 (rate scatters), omitted
plt.close(fig2)
print("wrote fig12_kinetic_decomposition.png")
