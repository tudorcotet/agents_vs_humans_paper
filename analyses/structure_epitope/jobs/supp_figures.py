"""Supplementary figures S6 (method robustness) and S8 (rules compliance)."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats
ROOT = Path("/Users/amyhe/Desktop/trem2_2026/agents_vs_humans_paper")
RES = ROOT / "analyses/structure_epitope/results"
FIG = ROOT / "analyses/structure_epitope/figures"
AGENT, HUMAN, BIND, NON, INK = "#30C5F5", "#1FE48F", "#00D9FF", "#5C6773", "#0F1419"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK, "axes.linewidth": 0.6,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 300})
meta = pd.read_parquet(ROOT / "analyses/structure_epitope/data/designs_canonical.parquet")

# ================= S6: method robustness =================
fp = pd.read_parquet(RES / "footprints_long.parquet")
scr = set(meta[meta.p_screened == True].design_id)
fps = fp[fp.design_id.isin(scr) & (fp.region == "IgSF")]

fig, axs = plt.subplots(1, 3, figsize=(12, 3.6))
# (a) per-residue contact frequency: Protenix vs Chai
freq = fps.groupby(["predictor", "uniprot"])["primary"].mean().unstack("predictor")
r = stats.spearmanr(freq["protenix"], freq["chai"])[0]
axs[0].scatter(freq["protenix"], freq["chai"], s=14, color="#33C4FF", edgecolor=INK, linewidth=0.3)
axs[0].plot([0, 1], [0, 1], "--", color=INK, lw=0.6)
axs[0].set_xlabel("Protenix contact freq"); axs[0].set_ylabel("Chai-1 contact freq")
axs[0].set_title(f"(a) per-residue agreement\nSpearman ρ={r:.2f} (each pt = 1 TREM2 residue)", fontsize=8, loc="left")
# (b) per-design Jaccard between predictors
ag = pd.read_csv(RES / "epitope_predictor_agreement.csv")
axs[1].hist([ag[ag.is_hit==True].jaccard.dropna(), ag[ag.is_hit==False].jaccard.dropna()],
            bins=np.linspace(0, 1, 16), color=[BIND, NON], label=["binder", "non-binder"], stacked=True)
axs[1].axvline(ag.jaccard.median(), color=INK, ls="--", lw=0.8)
axs[1].set_xlabel("Protenix–Chai footprint Jaccard"); axs[1].set_ylabel("designs")
axs[1].set_title(f"(b) predictor footprint overlap\nmedian J={ag.jaccard.median():.2f}", fontsize=8, loc="left")
axs[1].legend(frameon=False, fontsize=7)
# (c) cutoff sensitivity: footprint size per design at 4.0/4.5/5.0 (consensus = both predictors)
sizes = {}
for cut in ["contact_40", "contact_45", "contact_50"]:
    piv = fps.pivot_table(index=["design_id","uniprot"], columns="predictor", values=cut, aggfunc="first").fillna(False)
    cons = (piv.get("protenix", False) & piv.get("chai", False))
    sizes[cut.replace("contact_","")+" Å"] = cons.groupby(level=0).sum().reindex(sorted(scr)).fillna(0)
sz = pd.DataFrame({("%s"%k[:3]): v for k, v in [("4.0 Å", sizes["40 Å"]),("4.5 Å", sizes["45 Å"]),("5.0 Å", sizes["50 Å"])]})
bp = axs[2].boxplot([sizes["40 Å"], sizes["45 Å"], sizes["50 Å"]], tick_labels=["4.0","4.5","5.0"], patch_artist=True, widths=0.6)
for patch in bp["boxes"]: patch.set_facecolor("#9EDFFF")
for m in bp["medians"]: m.set_color(INK)
axs[2].set_xlabel("contact cutoff (Å)"); axs[2].set_ylabel("consensus footprint size (residues)")
axs[2].set_title("(c) cutoff sensitivity\n(consensus footprint)", fontsize=8, loc="left")
fig.suptitle("Fig S5. Method robustness — Protenix vs Chai-1 epitope agreement and cutoff sensitivity (screened, n=100). "
             "Boltz-2/AF2M complex CIFs excluded (corrupt).", fontsize=8.5, x=0.01, ha="left")
fig.tight_layout(rect=[0,0,1,0.95]); fig.savefig(FIG / "figS5_method_robustness.png", bbox_inches="tight")
plt.close(fig)

# ================= S8: rules compliance =================
rc = pd.read_csv(RES / "rules_compliance.csv")
rc["max_full_locid"] = rc[["scFv2_full_locid","scFv4_full_locid"]].max(axis=1)
fig, axs = plt.subplots(1, 3, figsize=(12, 3.6))
# (a) max CDR-LCS distribution, flag threshold 5
vals, counts = np.unique(rc.max_cdr_lcs, return_counts=True)
axs[0].bar(vals, counts, color=["#FF4628" if v>=5 else "#5C6773" for v in vals])
axs[0].axvline(4.5, color=INK, ls="--", lw=0.8); axs[0].set_xlabel("longest common substring vs any CDR (aa)")
axs[0].set_ylabel("designs"); axs[0].set_title("(a) CDR LCS distribution\n(red = flag ≥5 aa)", fontsize=8, loc="left")
# (b) max local identity to scFv-2/-4, flag region
axs[1].hist(rc.max_full_locid, bins=20, color="#9EDFFF", edgecolor=INK, linewidth=0.3)
axs[1].axvline(70, color="#FF4628", ls="--", lw=0.8, label="flag ≥70% over ≥10 aa")
axs[1].set_xlabel("max local identity to scFv-2/-4 (%)"); axs[1].set_ylabel("designs")
axs[1].set_title("(b) local identity to reference scFvs", fontsize=8, loc="left"); axs[1].legend(frameon=False, fontsize=7)
# (c) flagged designs table-as-scatter: LCS vs local id, label literature copies
flag = rc[rc.FLAG_any]
axs[2].scatter(rc.max_cdr_lcs, rc.max_full_locid, s=12, color="#CCC")
axs[2].scatter(flag.max_cdr_lcs, flag.max_full_locid, s=40, color="#FF4628", edgecolor=INK, zorder=3)
for _, row in flag.iterrows():
    axs[2].annotate(str(row.design_id), (row.max_cdr_lcs, row.max_full_locid), fontsize=7,
                    xytext=(3,3), textcoords="offset points")
axs[2].set_xlabel("max CDR LCS (aa)"); axs[2].set_ylabel("max local identity (%)")
axs[2].set_title("(c) 7 flagged designs (antibody scaffolds)\n77/79/80 = literature copies", fontsize=8, loc="left")
fig.suptitle("Fig S8. Rules compliance — 134/141 designs comfortably compliant; 7 human antibody-scaffold designs flagged "
             "(design 79 = 21.7 nM binder with exact scFv-4 HCDR3).", fontsize=8.5, x=0.01, ha="left")
fig.tight_layout(rect=[0,0,1,0.95]); None  # dropped former-FigS8 (rules compliance), omitted
plt.close(fig)
print("wrote figS5_method_robustness.png")
print(f"S6a per-residue Protenix–Chai Spearman ρ = {r:.3f}; S6c median footprint 4.0/4.5/5.0 = "
      f"{np.median(sizes['40 Å']):.0f}/{np.median(sizes['45 Å']):.0f}/{np.median(sizes['50 Å']):.0f}")
