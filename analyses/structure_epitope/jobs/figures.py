"""Analytical figures for the Structure & Epitope section (2D panels; 3D renders pending PyMOL)."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
ROOT = Path("/Users/amyhe/Desktop/trem2_2026/agents_vs_humans_paper")
RES = ROOT / "analyses/structure_epitope/results"
FIG = ROOT / "analyses/structure_epitope/figures"; FIG.mkdir(exist_ok=True)
AGENT, HUMAN = "#30C5F5", "#1FE48F"; BIND, NON = "#00D9FF", "#5C6773"; INK = "#0F1419"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK, "axes.linewidth": 0.6,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 300})
meta = pd.read_parquet(ROOT / "analyses/structure_epitope/data/designs_canonical.parquet")

# ---------- Fig 8 (epitope landscape, panel e): per-residue frequency ----------
freq = pd.read_csv(RES / "epitope_per_residue_freq.csv")
per = pd.read_csv(RES / "epitope_HPerera_perresidue.csv")
sig = set(per[per.q_bh < 0.05].uniprot)
f = freq[(freq.uniprot >= 19) & (freq.uniprot <= 132)].sort_values("uniprot")
fig, ax = plt.subplots(figsize=(9, 3.2))
ax.bar(f.uniprot, f.freq_binder, width=0.9, color=BIND, label="binders", alpha=0.9)
ax.plot(f.uniprot, f.freq_nonbinder, color=NON, lw=1.2, label="non-binders", drawstyle="steps-mid")
for u in sig:
    ax.text(u, 1.02, "*", ha="center", va="bottom", color="#FF4628", fontsize=11, fontweight="bold")
regions = {"hydrophobic tip": [44,69,70,71,74,89], "basic patch": [46,47,62,76,77,87,122,123]}
for name, res, y, col in [("tip", regions["hydrophobic tip"], -0.08, "#FFB547"),
                          ("basic patch", regions["basic patch"], -0.14, "#9B5DE5")]:
    for r in res:
        if 19 <= r <= 132: ax.plot([r], [y], marker="s", ms=3, color=col)
ax.set_xlabel("TREM2 residue (UniProt Q9NZC2)"); ax.set_ylabel("contact frequency")
ax.set_title("Fig 8e. Per-residue TREM2 engagement — binders vs non-binders (consensus footprint; * = FDR q<0.05)",
             fontsize=8.5, loc="left")
ax.set_ylim(-0.18, 1.1); ax.legend(loc="upper right", frameon=False, fontsize=8)
ax.text(0.01, -0.16, "▪ tip  ▪ basic patch", transform=ax.transAxes, fontsize=7, color="#666")
fig.tight_layout(); fig.savefig(FIG / "fig08e_epitope_landscape.png", bbox_inches="tight")
plt.close(fig)

# ---------- Fig 9: interface features separate binders ----------
avi = pd.read_csv(RES / "phase2_active_vs_inactive.csv")
desc = pd.read_csv(RES / "interface_descriptors.csv").merge(
    meta[["design_id","is_hit","p_expressed","cohort"]], on="design_id")
exp = desc[desc.p_expressed == True]
fig, axs = plt.subplots(1, 3, figsize=(11, 3.4))
# (a) effect sizes
a = avi.sort_values("cliffs_delta")
cols = ["#FF4628" if q < 0.05 else "#9EA2AF" for q in a.q_bh]
axs[0].barh(a.descriptor, a.cliffs_delta, color=cols)
axs[0].axvline(0, color=INK, lw=0.5); axs[0].set_xlabel("Cliff's δ (binder − non-binder)")
axs[0].set_title("(a) interface descriptor effect sizes\n(red = FDR q<0.05)", fontsize=8, loc="left")
# (b) BSA box
d0 = exp[exp.is_hit==True].bsa.dropna(); d1 = exp[exp.is_hit==False].bsa.dropna()
bp = axs[1].boxplot([d1, d0], labels=["non-binder","binder"], patch_artist=True, widths=0.6)
for patch, c in zip(bp["boxes"], [NON, BIND]): patch.set_facecolor(c); patch.set_alpha(0.85)
for m in bp["medians"]: m.set_color(INK)
axs[1].set_ylabel("buried surface area (Å²)")
axs[1].set_title("(b) interface area\nAUROC 0.71, q=0.007", fontsize=8, loc="left")
# (c) AUROC bars
a2 = avi.sort_values("AUROC")
axs[2].barh(a2.descriptor, a2.AUROC, color="#33C4FF"); axs[2].axvline(0.5, color=INK, lw=0.5, ls="--")
axs[2].set_xlim(0.4, 0.8); axs[2].set_xlabel("AUROC (binder classifier)")
axs[2].set_title("(c) univariate AUROC", fontsize=8, loc="left")
fig.suptitle("Fig 9. Interface size — not epitope or chemistry — separates binders from non-binders (P_expressed, n=89)",
             fontsize=9, x=0.01, ha="left")
fig.tight_layout(rect=[0,0,1,0.96]); fig.savefig(FIG / "fig09_interface_features.png", bbox_inches="tight")
plt.close(fig)

# ---------- Fig 12: metrics vs affinity ----------
d = desc.merge(meta[["design_id","p_kd","pkd_arith_mean"]], on="design_id")
df6 = pd.read_parquet(ROOT/"data/designs.parquet")[["design_id","px_ipsae_d0chn_max","sequence_length"]]
d = d.merge(df6, on="design_id"); kd = d[d.p_kd==True]
from scipy import stats
panels = [("n_iface_binder","interface size (residues)"),("bsa","buried area (Å²)"),
          ("sequence_length","binder length (aa)"),("px_ipsae_d0chn_max","Protenix ipSAE (confidence)")]
fig, axs = plt.subplots(1, 4, figsize=(13, 3.2))
for ax, (col, lab) in zip(axs, panels):
    s = kd.dropna(subset=[col,"pkd_arith_mean"])
    cc = [AGENT if c=="agent" else HUMAN for c in s.cohort]
    ax.scatter(s[col], s.pkd_arith_mean, c=cc, s=28, edgecolor=INK, linewidth=0.3)
    rho, p = stats.spearmanr(s[col], s.pkd_arith_mean)
    ax.set_xlabel(lab); ax.set_title(f"ρ={rho:.2f} (p={p:.3f})", fontsize=8.5)
axs[0].set_ylabel("pKd = −log10(Kd)")
fig.suptitle("Fig 11. Structural interface size & length track affinity; confidence scores do not (P_kd, n=36; ● human ● agent)",
             fontsize=9, x=0.01, ha="left")
fig.tight_layout(rect=[0,0,1,0.95]); fig.savefig(FIG / "fig11_metrics_vs_affinity.png", bbox_inches="tight")
plt.close(fig)

# ---------- Fig 11: fold x cohort (agent monoculture) ----------
fold = pd.read_csv(RES / "fold_topology.csv")
scr = fold[fold.p_screened==True]
order = [">=4-helix","3-helix bundle","2-helix/hairpin","alpha+beta","beta/mixed","loop/other"]
ct = pd.crosstab(scr.cohort, scr.topology).reindex(columns=order, fill_value=0)
frac = ct.div(ct.sum(axis=1), axis=0)
fig, ax = plt.subplots(figsize=(8, 2.8))
left = np.zeros(len(frac)); cmap = ["#142933","#36B7F6","#33C4FF","#9EDFFF","#FFB547","#5C6773"]
for i, topo in enumerate(order):
    ax.barh(frac.index, frac[topo], left=left, color=cmap[i], label=topo)
    left += frac[topo].values
ax.set_xlabel("fraction of screened designs"); ax.legend(ncol=3, fontsize=7, frameon=False, loc="upper center", bbox_to_anchor=(0.5,-0.25))
ax.set_title("Fig 11. Binder fold by cohort — agents concentrate on helical bundles (structural tool-monoculture)",
             fontsize=8.5, loc="left")
fig.tight_layout(); None  # dropped former-Fig11 (fold x cohort), omitted from final set
plt.close(fig)

print("wrote figures:", sorted(p.name for p in FIG.glob("*.png")))
