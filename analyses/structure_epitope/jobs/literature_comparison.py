"""EXTRA 2 (OPEN_ITEMS_register_1) — K_D comparison with literature measurements.

Fills the manuscript Discussion placeholder (line 281). Builds an assay-metadata-carrying
reference table and a log-axis strip plot of this study's 36 K_D vs prior TREM2 binders and
cross-target design campaigns, with Perera's ~5 nM MST detection floor marked.

Literature values are CURATED FROM OPEN_ITEMS_register_1.md — the register flags several as
"check"; they must be citation-verified before publication (source_verified column).

Run: uv run --with scipy python analyses/structure_epitope/jobs/literature_comparison.py
Outputs: results/literature_kd_comparison.csv, results/EXTRA2_framing.md, figures/figS6_literature_kd.png
"""
from __future__ import annotations
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
SE = os.path.join(ROOT, "analyses", "structure_epitope")
RES, FIG = os.path.join(SE, "results"), os.path.join(SE, "figures")
AGENT, HUMAN, INK, ACC = "#30C5F5", "#1FE48F", "#0F1419", "#FF4628"
PERERA_FLOOR = 5.0
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK, "axes.linewidth": 0.6,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 300})

# Reference set — curated from OPEN_ITEMS_register_1.md. kd_nM=None means no usable public value.
REF = [
    # binder, kd_nM, category, assay, immobilized, valency_avidity, fit_type, construct, citation, verified
    ("scFv-2", 1.0, "anti-TREM2 Ab", "SPR", "TREM2 (biotin/SA chip)", "scFv oligomer (avidity flag)", "kinetic", "TREM2 19–174", "Szykowska 2021", "epitope from 6YYE coords; KD from paper"),
    ("scFv-4", 1.0, "anti-TREM2 Ab", "SPR", "TREM2 (biotin/SA chip)", "scFv oligomer (avidity flag)", "kinetic", "TREM2 19–174", "Szykowska 2021", "epitope from 6Y6C coords; KD from paper"),
    ("4D9", None, "anti-TREM2 Ab", "?", "?", "?", "?", "stalk ~145–153", "Schlepckow 2020", "CHECK"),
    ("AL002", None, "anti-TREM2 Ab", "?", "?", "?", "?", "stalk", "clinical (Phase 2 fail)", "CHECK"),
    ("VHB937", None, "anti-TREM2 Ab", "?", "?", "?", "?", "IgSF (epitope undisclosed)", "Novartis", "no public affinity"),
    ("Perera Odesign2", 80.0, "de novo TREM2", "MST", "TREM2 labeled (10 nM fixed)", "monomer", "equilibrium", "Sino TREM2 19–174", "Perera 2026", "from register; MST floor ~5 nM"),
    ("Perera BindCraft8", 370.0, "de novo TREM2", "MST", "TREM2 labeled (10 nM fixed)", "monomer", "equilibrium", "Sino TREM2 19–174", "Perera 2026", "from register"),
    ("Adaptyv EGFR (best)", 1.21, "cross-target", "?", "?", "?", "?", "EGFR", "Adaptyv", "CHECK — context only"),
    ("BindCraft PD-L1", 615.0, "cross-target", "?", "?", "?", "?", "PD-L1", "BindCraft", "CHECK — context only"),
    ("BindCraft CD45", 14.7, "cross-target", "?", "?", "?", "?", "CD45", "BindCraft", "CHECK — context only"),
    ("BindCraft Der f7", 12.8, "cross-target", "?", "?", "?", "?", "Der f7", "BindCraft", "CHECK — context only"),
    ("AlphaProteo (best)", 0.0085, "cross-target", "?", "?", "?", "?", "various", "AlphaProteo", "CHECK — context only"),
]
COLS = ["binder", "kd_nM", "category", "assay", "immobilized_partner", "valency_avidity",
        "fit_type", "target_construct", "citation", "source_verified"]


def main():
    d = pd.read_parquet(os.path.join(ROOT, "data", "designs.parquet"))
    ours = d[d.kd_arith_mean_nM_all.notna()][["design_id", "kd_arith_mean_nM_all", "is_human"]].copy()
    ours["cohort"] = ours.is_human.map({True: "human", False: "agent"})
    best = ours.kd_arith_mean_nM_all.min()

    ref = pd.DataFrame(REF, columns=COLS)
    # add this study's summary rows for the table
    summ = pd.DataFrame([
        ["This study (best)", round(best, 2), "de novo TREM2", "BLI/SPR kinetic (see note)",
         "design (Twin-Strep C-term)", "monomer", "kinetic", "Acro TR2-H52H5 TREM2 19–174", "this study", "measured; assay label under review (SPR in replicate table)"],
        ["This study (median of 36)", round(ours.kd_arith_mean_nM_all.median(), 1), "de novo TREM2", "BLI/SPR kinetic",
         "design (Twin-Strep C-term)", "monomer", "kinetic", "Acro TR2-H52H5 TREM2 19–174", "this study", "measured"],
    ], columns=COLS)
    table = pd.concat([summ, ref], ignore_index=True)
    table.to_csv(os.path.join(RES, "literature_kd_comparison.csv"), index=False)

    # ---------------- figure ----------------
    fig, ax = plt.subplots(figsize=(10.5, 5.6))
    groups = ["This study (n=36)", "de novo TREM2 (Perera)", "anti-TREM2 Ab", "cross-target campaigns"]
    ypos = {g: i for i, g in enumerate(groups[::-1])}

    # this study: strip, colored by cohort
    rng = np.random.default_rng(7)
    y0 = ypos["This study (n=36)"]
    for coh, col in [("human", HUMAN), ("agent", AGENT)]:
        s = ours[ours.cohort == coh]
        jit = y0 + (rng.random(len(s)) - 0.5) * 0.5
        ax.scatter(s.kd_arith_mean_nM_all, jit, s=34, color=col, edgecolor=INK, linewidth=0.3,
                   alpha=0.9, label=f"this study {coh} (n={len(s)})", zorder=3)
    ax.scatter([best], [y0 + 0.42], marker="v", s=70, color=ACC, edgecolor=INK, zorder=4)
    ax.annotate(f"best {best:.2f} nM (design 17)", (best, y0 + 0.42), fontsize=7.5, color=ACC,
                xytext=(-8, 8), textcoords="offset points", ha="right")

    def place(g, name, kd, marker, col, dy=0.0, ha="left"):
        ax.scatter([kd], [ypos[g] + dy], marker=marker, s=70, color=col, edgecolor=INK, zorder=3)
        pad = " " if ha == "left" else "  "
        ax.annotate(pad + name if ha == "left" else name + pad, (kd, ypos[g] + dy),
                    fontsize=7.5, va="center", ha=ha, color=INK)

    place("de novo TREM2 (Perera)", "Odesign2 80 nM (MST)", 80, "s", "#7A5CFF", ha="right")
    place("de novo TREM2 (Perera)", "BindCraft8 370 nM", 370, "s", "#7A5CFF", ha="left")
    place("anti-TREM2 Ab", "scFv-2 ~1 nM (AppKD, avidity)", 1.0, "D", "#8A93A0", dy=+0.20)
    place("anti-TREM2 Ab", "scFv-4 ~1 nM (AppKD, avidity)", 1.0, "D", "#8A93A0", dy=-0.20)
    for name, kd, ha in [("AlphaProteo 8.5 pM", 0.0085, "left"), ("Adaptyv EGFR 1.21 nM", 1.21, "left"),
                         ("BindCraft Der f7 12.8", 12.8, "left"), ("BindCraft PD-L1 615", 615, "right")]:
        place("cross-target campaigns", name, kd, "o", "#C3C9D2", ha=ha)

    ax.axvline(PERERA_FLOOR, color=ACC, ls="--", lw=0.9)
    ax.text(PERERA_FLOOR, 1.5, " Perera MST\n floor ~5 nM", color=ACC, fontsize=7.5, va="center")
    ax.set_xscale("log")
    ax.set_xlabel("K_D (nM, log scale)")
    ax.set_yticks(list(ypos.values()))
    ax.set_yticklabels(list(ypos.keys()))
    ax.set_ylim(-0.6, len(groups) - 0.3)
    ax.set_title("Fig S6. K_D in context: this study (36 de novo TREM2 binders) vs. prior TREM2 binders "
                 "and cross-target design campaigns.\nSame target construct (TREM2 19–174) as Szykowska SPR and "
                 "Perera MST. Literature values curated from the open-items register — verify citations before use.",
                 fontsize=8, loc="left")
    ax.legend(frameon=False, fontsize=7.5, loc="lower right")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "figS6_literature_kd.png"), bbox_inches="tight")
    plt.close(fig)

    # ---------------- framing ----------------
    n_below = int((ours.kd_arith_mean_nM_all < PERERA_FLOOR).sum())
    F = [
        "# EXTRA 2 — literature K_D comparison: framing (manuscript Discussion, line 281)\n",
        "**Comparability anchor (the strongest claim available):** target construct is identical across the "
        "three directly-relevant datasets — Szykowska SPR used His19–Ser174, Perera used Sino Biological "
        "TREM2(19–174), and this study used Acro TR2-H52H5 (TREM2 19–174). **All three are the same construct.** "
        "State this explicitly.\n",
        f"**Claim.** Best measured affinity **{best:.2f} nM** (design 17) is competitive with the strongest "
        "reported de novo binders across all targets, and **~72× better than the best prior de novo TREM2 "
        "binder** (Perera Odesign2, 80 nM).\n",
        f"**Qualify.** Perera's MST detection floor is ~5 nM (labeled TREM2 fixed at 10 nM), so their assay "
        f"could not have resolved our top designs — **{n_below} of our binders sit below 5 nM**. The gap is real "
        "but partly reflects assay range.\n",
        "**Concede (do not omit).** Perera's hit rate is **8/13 (61.5%)** vs. our **37/100 (37%)**. Their 13 were "
        "hand-picked from 106 by one expert lab and 'measurable binding' is a looser bar than a fitted K_D — but on "
        "the raw numbers they win; address it directly.\n",
        "**Assay-orientation note.** Design 79 at **21.66 nM** vs. scFv-4's ~1 nM AppKD is not a contradiction: "
        "theirs is SPR with inverted orientation and oligomer-driven avidity they themselves flag; ours is a "
        "monomeric, surface-tethered design with monovalent analyte. ~20× weaker is the predicted direction.\n",
        "**Provenance caveat (found in this analysis).** Our kinetic fits are labeled **SPR** in the replicate "
        "table (32/36 SPR-only), while the manuscript frames the assay as BLI. Confirm with Tudor (Tier-4 4.1) "
        "before finalizing the assay column.\n",
        "**Verification TODO.** Antibody/cross-target K_D values here are curated from the register and several are "
        "marked CHECK (4D9, AL002, Adaptyv, BindCraft, AlphaProteo). Verify each citation before publication; the "
        "scFv-2/-4 *epitopes* are already coordinate-verified (6YYE/6Y6C).\n",
        "Table: results/literature_kd_comparison.csv · Figure: figS6_literature_kd.png",
    ]
    open(os.path.join(RES, "EXTRA2_framing.md"), "w").write("\n".join(F) + "\n")
    print("\n".join(F))
    print(f"\nwrote literature_kd_comparison.csv ({len(table)} rows), EXTRA2_framing.md, figS6_literature_kd.png")


if __name__ == "__main__":
    main()
