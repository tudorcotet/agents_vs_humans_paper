"""Fig S8 — binder/non-binder discrimination (AUROC) by metric x predictor.

The section's central classification result as a figure: learned confidence scores (ipSAE, ipTM, mean pLDDT)
are strongly predictor-dependent, while the structural anchor (buried area) is nearly flat across predictors.
Grouped bars derived entirely from results/metric_predictor_table.csv (P_expressed, n=89). The source carries
point estimates + p-values but NO confidence intervals, so none are drawn (stated in the title). Boltz-2 — the
metric the competition actually selected on, and the weakest column — is hatched to stand out without a legend
lookup. ESMFold2 has no learned confidence scores, shown as an explicit 'n/a', not a missing/zero bar.

Matches jobs/figures.py conventions (palette family, font sizes, loc="left" titles). Nothing is hardcoded:
every bar height and the n come from the CSV. Run: uv run python jobs/fig_auroc_by_predictor.py
"""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from pathlib import Path

ROOT = Path("/Users/amyhe/Desktop/trem2_2026/agents_vs_humans_paper")
RES = ROOT / "analyses/structure_epitope/results"
FIG = ROOT / "analyses/structure_epitope/figures"; FIG.mkdir(exist_ok=True)
INK = "#0F1419"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK, "axes.linewidth": 0.6,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 300})

t = pd.read_csv(RES / "metric_predictor_table.csv")
# learned confidence scores first, structural anchor last: 'predictor-dependent -> robust' reads left->right
METRICS = [("ipSAE", "ipSAE"), ("ipTM", "ipTM"), ("pLDDT (mean)", "mean pLDDT"),
           ("buried area (BSA)", "buried area\n(interface size)")]
PREDS = ["Protenix", "Chai", "Boltz-2", "AF2M", "ESMFold2"]
COL = {"Protenix": "#1B6FB3", "Chai": "#2E9E8F", "Boltz-2": "#E8862C", "AF2M": "#7E86A6", "ESMFold2": "#B9C0CC"}
n_disc = int(t.n_disc.dropna().iloc[0])

def auroc(metric_csv, pred):
    r = t[(t.metric == metric_csv) & (t.predictor == pred)]
    return float(r.auroc.iloc[0]) if len(r) and pd.notna(r.auroc.iloc[0]) else np.nan

nb = len(PREDS); bw = 0.16
xg = np.arange(len(METRICS)) * (nb * bw + 0.22)   # group centres
fig, ax = plt.subplots(figsize=(9.2, 3.7))
ax.axhline(0.5, color=INK, lw=0.9, ls="--", zorder=1)
ax.text(xg[-1] + nb * bw / 2, 0.505, "0.5 = no discrimination", ha="right", va="bottom",
        fontsize=7, color="#5C6773", style="italic")

for gi, (mcsv, _) in enumerate(METRICS):
    for pi, pred in enumerate(PREDS):
        x = xg[gi] + (pi - (nb - 1) / 2) * bw
        v = auroc(mcsv, pred)
        if np.isnan(v):
            ax.text(x, 0.508, "n/a", rotation=90, ha="center", va="bottom", fontsize=6.5, color="#9096A0")
            continue
        ax.bar(x, v, bw * 0.9, color=COL[pred], edgecolor=INK, linewidth=0.5,
               hatch="////" if pred == "Boltz-2" else None, zorder=3)
        ax.text(x, v + 0.008, f"{v:.2f}", ha="center", va="bottom", fontsize=6.1, color=INK)

# region separator + neutral headers (the two takeaways read off the bar heights, no legend needed)
sep = (xg[2] + xg[3]) / 2
ax.axvline(sep, color="#D0D4DB", lw=1.0, ls=":", zorder=1)
ax.text((xg[0] + xg[2]) / 2, 0.88, "learned confidence scores", ha="center", fontsize=8, color="#5C6773")
ax.text(xg[3], 0.88, "structural", ha="center", fontsize=8, color="#5C6773")

ax.set_xticks(xg); ax.set_xticklabels([m[1] for m in METRICS])
ax.set_ylim(0, 0.92); ax.set_ylabel("AUROC (binder vs non-binder)")
handles = [Patch(facecolor=COL[p], edgecolor=INK, linewidth=0.5,
                 hatch="////" if p == "Boltz-2" else None,
                 label=p + ("  (competition selection)" if p == "Boltz-2" else "")) for p in PREDS]
ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.32), ncol=5,
          frameon=False, fontsize=7.5, handlelength=1.2, columnspacing=1.1)
fig.suptitle(f"Fig S8. Binder/non-binder discrimination (AUROC) by metric × predictor — learned confidence "
             f"scores are predictor-dependent, interface size is not (P_expressed, n={n_disc}; point estimates, "
             f"no CIs in source)", fontsize=8.5, x=0.01, ha="left")
fig.tight_layout(rect=[0, 0.02, 1, 0.96])
fig.savefig(FIG / "figS8_auroc_by_predictor.png", bbox_inches="tight")
plt.close(fig)
print("wrote figS8_auroc_by_predictor.png (n_disc =", n_disc, ")")
