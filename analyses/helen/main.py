"""
Manuscript figures for the muni x Adaptyv TREM2 hackathon paper.

Each figure is one self-contained SVG in figures/paper/, styled to match the
hand-authored blog figures (figures/blog/*.html): Agent #30C5F5 (cyan),
Human #142933 (navy ink), Geist / GT Pressura Extended fonts.

ESM-2 (Fig 5, Fig 7 pLL curve), Foldseek and epitope figures read caches
written by the sibling scripts (esm2_compute.py, foldseek_cluster.py,
epitope_cluster.py); those builders are skipped gracefully if a cache is
absent. Fig 3 (BLI sensorgrams) is excluded by request.

Outputs: figures/paper/{fig1_ipsae_distribution, fig2_expression_binding,
fig4_design_methods, fig6_sequence_length, fig7_metric_roc,
fig_identity_heatmap, fig_foldseek_clustering, fig_epitope_regions,
fig5_esm2_umap}.svg ; analyses/helen/{report.md, summary.json}.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import auc, roc_curve

from scripts.plotting import apply_theme
from scripts.utils import load_designs, repo_root

OUT = repo_root() / "figures" / "paper"
HERE = Path(__file__).resolve().parent

# Blog brand tokens (see figures/blog/*.html CSS custom properties).
AGENT = "#30C5F5"
HUMAN = "#142933"
AGENT_SOFT = "#9EDFFF"
HUMAN_SOFT = "#E5E7EB"
INK = "#0F1419"
GRID = "#D9DEE3"

# 5-category design taxonomy (Cotet et al. 2025, as listed in the Notion draft).
# This crosswalk is an interpretation of design_method_normalized; documented
# here and in report.md so it can be challenged.
METHOD_TO_CATEGORY = {
    "AF+PXDesign+Boltz2": "De novo",
    "BindCraft": "De novo",            # Notion explicitly buckets BindCraft as De novo
    "BoltzGen": "De novo",
    "BoltzGen+Boltz2": "De novo",
    "BoltzGen+Boltz2+ipSAE": "De novo",
    "Foundry": "De novo",              # Baker-lab RFdiffusion3 pipeline
    "PPIFLOW+MPNN+FAMPNN": "De novo",
    "PXDesign": "De novo",
    "PXDesign+MPNN": "De novo",
    "PXDesign+Protenix+AF2": "De novo",
    "Protpardelle": "De novo",
    "RFDiffusion+LigandMPNN+Boltz": "De novo",
    "RFDiffusion+MPNN+Boltz2": "De novo",
    "RFDiffusion2+MPNN+Boltz2": "De novo",
    "RFDiffusion3+MPNN+RosettaFold": "De novo",
    "RFPeptides": "De novo",
    "AFHall+MPNN+PyR+Boltz2": "Hallucination",
    "Mosaic": "Optimized",             # composite-objective optimization wrapper
    "Struct_Evo": "Optimized",         # structure-guided evolutionary optimization
    "MPNN+Boltz2": "Diversified",      # ProteinMPNN inverse-folding redesign
    "evo+ESM": "Diversified",          # pLM-guided sequence diversification
}
CATEGORY_ORDER = ["De novo", "Hallucination", "Optimized", "Diversified", "Rational/Hybrid"]
# Solid mid-tone palette taken from the Adaptyv "Design Class Distribution"
# blog pie (image #4): no dark outline, thin white separators between stacks.
CATEGORY_COLOR = {
    "De novo": "#8FD3EA",        # light sky blue (Miniprotein)
    "Hallucination": "#8B86CE",  # lavender (Other)
    "Optimized": "#3C9BD9",      # medium blue (Peptide)
    "Diversified": "#6FA995",    # teal-green (Nanobody)
    "Rational/Hybrid": "#D3D3D3",  # light grey (scFv)
}

# Bars/histograms: light fill + dark outline, same two cohort hues as the
# sequence-length figure (Human navy, Agent cyan).
COHORT_FILL = {"Human": HUMAN_SOFT, "Agent": AGENT_SOFT}
COHORT_EDGE = {"Human": HUMAN, "Agent": AGENT}


def _style() -> None:
    apply_theme()
    mpl.rcParams["font.family"] = "sans-serif"
    mpl.rcParams["font.sans-serif"] = ["Geist", "Helvetica", "Arial", "DejaVu Sans"]
    mpl.rcParams["svg.fonttype"] = "none"  # keep text selectable + brand-font-aware
    mpl.rcParams["axes.edgecolor"] = INK
    mpl.rcParams["text.color"] = INK
    mpl.rcParams["axes.labelcolor"] = INK
    mpl.rcParams["xtick.color"] = INK
    mpl.rcParams["ytick.color"] = INK
    # Airy Plotly-style frame (Nipah blog): faint horizontal grid behind data,
    # no left spine — the grid carries the y-scale.
    mpl.rcParams["axes.grid"] = True
    mpl.rcParams["axes.grid.axis"] = "y"
    mpl.rcParams["axes.axisbelow"] = True
    mpl.rcParams["grid.color"] = "#E5E7EB"
    mpl.rcParams["grid.linewidth"] = 0.6
    mpl.rcParams["axes.spines.left"] = False
    mpl.rcParams["ytick.major.size"] = 0


# Soft neutral title grey (blog uses a medium-weight near-black, not pure ink).
_TITLE_GREY = "#2A2F35"
_SUBTITLE = "#36B7F6"  # brand cyan-deep — blog panel-subtitle accent


def _title(ax, text: str) -> None:
    ax.set_title(text, fontfamily=["GT Pressura Extended", "Geist", "DejaVu Sans"],
                 fontsize=10, fontweight="medium", color=_TITLE_GREY,
                 loc="center", pad=16)


def _subtitle(ax, text: str) -> None:
    """Panel sub-title in the blog's muted-cyan accent."""
    ax.set_title(text, fontfamily=["Geist", "DejaVu Sans"], fontsize=8.5,
                 fontweight="medium", color=_SUBTITLE, loc="center", pad=8)


_FONT_DIR = repo_root() / "figures" / "blog" / "fonts"
# (family, weight, woff2 file) — exact families the SVG text declares.
_FONT_FACES = [
    ("Geist", 400, "Geist-Regular.woff2"),
    ("Geist", 500, "Geist-Medium.woff2"),
    ("Geist", 700, "Geist-Bold.woff2"),
    ("GT Pressura Extended", 400, "GT-Pressura-Extended-Regular.woff2"),
]


def _font_style_block() -> str:
    import base64

    faces = []
    for family, weight, fname in _FONT_FACES:
        b64 = base64.b64encode((_FONT_DIR / fname).read_bytes()).decode("ascii")
        faces.append(
            f"@font-face{{font-family:'{family}';font-style:normal;"
            f"font-weight:{weight};src:url(data:font/woff2;base64,{b64}) "
            f"format('woff2');}}"
        )
    return ('<style type="text/css"><![CDATA[\n' + "\n".join(faces)
            + "\n]]></style>")


def _embed_fonts(svg_path: Path) -> None:
    """Inline the blog woff2 as @font-face so the SVG is self-contained.

    NOTE: GT Pressura Extended is a paid display face; embedding it in a
    distributed artifact is a licensing consideration (user-approved here).
    """
    import re

    svg = svg_path.read_text()
    style = _font_style_block()
    svg = re.sub(r"(<svg\b[^>]*>)", lambda m: m.group(1) + "\n" + style,
                 svg, count=1)
    svg_path.write_text(svg)


def _save(fig, name: str) -> Path:
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / f"{name}.svg"
    fig.savefig(path, format="svg", bbox_inches="tight")
    plt.close(fig)
    _embed_fonts(path)
    return path


def _wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float, float]:
    """Wilson score interval for a binomial proportion. Returns (p, lo, hi)."""
    if n == 0:
        return (math.nan, math.nan, math.nan)
    p = k / n
    denom = 1 + z**2 / n
    centre = (p + z**2 / (2 * n)) / denom
    half = (z * math.sqrt(p * (1 - p) / n + z**2 / (4 * n**2))) / denom
    return p, centre - half, centre + half


def fig1_ipsae(df) -> Path:
    """ipSAE distribution, human vs agent, with the top-100 selection cutoff."""
    fig, ax = plt.subplots(figsize=(4.6, 3.4))
    groups = [("Human", HUMAN), ("Agent", AGENT)]
    data = []
    for label, _ in groups:
        is_h = label == "Human"
        vals = df.loc[df.is_human.fillna(False) == is_h, "submitted_ipsae"].dropna()
        data.append(vals.to_numpy())

    parts = ax.violinplot(data, showextrema=False, widths=0.75)
    for body, (_, colour) in zip(parts["bodies"], groups, strict=False):
        body.set_facecolor(colour)
        body.set_alpha(0.30)
        body.set_edgecolor(colour)
        body.set_linewidth(1.0)
    bp = ax.boxplot(data, widths=0.18, showfliers=False, patch_artist=True)
    # Human box drawn in black (per request); agent keeps the cohort cyan.
    box_colours = ["black", AGENT]
    for i, colour in enumerate(box_colours):
        bp["boxes"][i].set(facecolor="white", edgecolor=colour, linewidth=1.2)
        bp["medians"][i].set(color=colour, linewidth=1.6)
        for j in (2 * i, 2 * i + 1):  # 2 whiskers + 2 caps per box
            bp["whiskers"][j].set_color(colour)
            bp["caps"][j].set_color(colour)

    cutoff = df.loc[df.submitted_to_lab.fillna(False), "submitted_ipsae"].min()
    ax.axhline(cutoff, ls="--", lw=1.0, color="#5C6773")
    ax.text(2.46, cutoff, f"  top-100 cutoff\n  ipSAE = {cutoff:.3f}", va="center",
            ha="left", fontsize=7, color="#5C6773")

    ax.set_xticks([1, 2])
    ax.set_xticklabels([f"Human\n(n={len(data[0])})", f"Agent\n(n={len(data[1])})"])
    ax.set_ylabel("Boltz-2 ipSAE")
    ax.set_xlim(0.4, 3.05)
    _title(ax, "ipSAE distribution by cohort")
    return _save(fig, "fig1_ipsae_distribution")


def fig2_expression_binding(df) -> Path:
    """Expression rate and hit rate by cohort, Wilson 95% CI."""
    scr = df[df.submitted_to_lab.fillna(False)]
    fig, ax = plt.subplots(figsize=(4.6, 3.4))
    metrics = ["Expression rate", "Hit rate"]
    x = np.arange(len(metrics))
    width = 0.34

    for offset, label in [(-0.20, "Human"), (0.20, "Agent")]:
        sub = scr[scr.is_human.fillna(False) == (label == "Human")]
        n = len(sub)
        expr_k = int(sub.expressed.fillna(False).sum())
        hit_k = int(sub.is_hit.fillna(False).sum())
        stats = [_wilson(expr_k, n), _wilson(hit_k, n)]
        heights = [s[0] for s in stats]
        lo = [s[0] - s[1] for s in stats]
        hi = [s[2] - s[0] for s in stats]
        ax.bar(x + offset, heights, width, label=f"{label} (n={n})",
               color=COHORT_FILL[label], edgecolor=COHORT_EDGE[label],
               linewidth=1.0)
        ax.errorbar(x + offset, heights, yerr=[lo, hi], fmt="none",
                    ecolor=INK, elinewidth=0.9, capsize=3)
        for xi, h in zip(x + offset, heights, strict=False):
            ax.text(xi, h + 0.02, f"{h:.0%}", ha="center", va="bottom", fontsize=7)

    ax.set_xticks(x)
    ax.set_xticklabels(metrics)
    ax.set_ylabel("Proportion of lab-tested designs")
    ax.set_ylim(0, 1.05)
    ax.legend(frameon=False, fontsize=7, loc="upper right")
    _title(ax, "Expression and binding by cohort")
    return _save(fig, "fig2_expression_binding")


def fig4_design_methods(df) -> Path:
    """Stacked design-method taxonomy by cohort (all 141 submissions)."""
    d = df.copy()
    d["category"] = d.design_method_normalized.map(METHOD_TO_CATEGORY).fillna("Rational/Hybrid")
    fig, ax = plt.subplots(figsize=(4.6, 3.4))
    cohorts = [("Human", df.is_human.fillna(False)), ("Agent", ~df.is_human.fillna(False))]
    x = np.arange(len(cohorts))

    bottoms = np.zeros(len(cohorts))
    for cat in CATEGORY_ORDER:
        vals = []
        for _, mask in cohorts:
            vals.append(int((d[mask].category == cat).sum()))
        vals = np.array(vals)
        ax.bar(x, vals, 0.55, bottom=bottoms, label=cat,
               color=CATEGORY_COLOR[cat], edgecolor="none")
        for xi, v, b in zip(x, vals, bottoms, strict=False):
            if v > 0:
                ax.text(xi, b + v / 2, str(v), ha="center", va="center",
                        fontsize=7, color=INK)
        bottoms += vals

    ax.set_xticks(x)
    ax.set_xticklabels([f"Human\n(n={int(c[1].sum())})" for c in cohorts[:1]]
                       + [f"Agent\n(n={int(cohorts[1][1].sum())})"])
    ax.set_ylabel("Designs submitted")
    ax.legend(frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    _title(ax, "Design-method taxonomy by cohort")
    return _save(fig, "fig4_design_methods")


def fig6_sequence_length(df) -> Path:
    """Overlaid sequence-length distributions, human vs agent (all 141)."""
    from scipy.stats import gaussian_kde

    fig, ax = plt.subplots(figsize=(4.6, 3.4))
    bins = np.linspace(df.sequence_length.min(), df.sequence_length.max(), 26)
    xs = np.linspace(df.sequence_length.min(), df.sequence_length.max(), 200)
    for label in ("Human", "Agent"):
        v = df.loc[df.is_human.fillna(False) == (label == "Human"), "sequence_length"].dropna()
        ax.hist(v, bins=bins, density=True, color=COHORT_FILL[label],
                alpha=0.55, edgecolor=COHORT_EDGE[label], linewidth=1.0)
        kde = gaussian_kde(v)
        ax.plot(xs, kde(xs), color=COHORT_EDGE[label], lw=1.8,
                label=f"{label} (n={len(v)}, med={v.median():.0f} aa)")
        ax.axvline(v.median(), color=COHORT_EDGE[label], ls="--", lw=1.0,
                   alpha=0.8)

    ax.set_xlabel("Sequence length (aa)")
    ax.set_ylabel("Density")
    ax.legend(frameon=False, fontsize=7, loc="upper right")
    _title(ax, "Sequence length by cohort")
    return _save(fig, "fig6_sequence_length")


def _load_esm2():
    """Return (design_id->emb, design_id->pll) from the ESM-2 cache, or (None, None)."""
    cache = HERE / "esm2_cache.npz"
    if not cache.exists():
        return None, None
    z = np.load(cache)
    ids = z["design_id"].astype(int)
    emb = {int(i): e for i, e in zip(ids, z["emb"], strict=False)}
    pll = {int(i): float(p) for i, p in zip(ids, z["pll"], strict=False)}
    return emb, pll


def fig7_metric_roc(df) -> Path:
    """ROC of in-silico metrics as binder predictors (screened + expressed)."""
    s = df[df.submitted_to_lab.fillna(False) & df.expressed.fillna(False)].copy()
    s = s[s.binding_label.isin(["binder", "weak", "non_binder"])]
    y = s.binding_label.isin(["binder", "weak"]).astype(int).to_numpy()

    _, pll = _load_esm2()
    if pll is not None:
        s["esm2_pll"] = s.design_id.map(pll)

    # (column, display label, higher_is_better)
    specs = [
        ("submitted_ipsae", "ipSAE", True),
        ("boltz2_iptm", "ipTM", True),
        ("boltz2_plddt", "pLDDT (binder)", True),
        ("pb_boltz2_complex_pde", "iPAE proxy (Boltz-2 PDE)", False),
    ]
    colours = [AGENT, "#36B7F6", HUMAN, "#5C6773", "#1FE48F"]
    if pll is not None:
        specs.append(("esm2_pll", "ESM-2 pLL (650M)", True))

    fig, ax = plt.subplots(figsize=(4.4, 4.0))
    ax.plot([0, 1], [0, 1], ls="--", lw=0.9, color=GRID)
    for (col, lab, hib), colour in zip(specs, colours, strict=False):
        m = s[col].notna().to_numpy()
        if m.sum() < 5:
            continue
        score = s.loc[m, col].to_numpy().astype(float)
        if not hib:
            score = -score
        fpr, tpr, _ = roc_curve(y[m], score)
        ax.plot(fpr, tpr, lw=1.8, color=colour, label=f"{lab}  AUC={auc(fpr, tpr):.2f}")

    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_aspect("equal")
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    _title(ax, f"Metric ROC (binder vs non-binder, n={len(s)})")
    return _save(fig, "fig7_metric_roc")


def fig5_esm2_umap(df) -> Path | None:
    """UMAP of ESM-2 650M embeddings; one SVG, two panels (cohort, team)."""
    emb, _ = _load_esm2()
    if emb is None:
        return None
    import umap

    d = df[["design_id", "is_human", "team"]].copy()
    d = d[d.design_id.isin(emb)].sort_values("design_id").reset_index(drop=True)
    X = np.vstack([emb[i] for i in d.design_id])
    reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, metric="cosine", random_state=0)
    xy = reducer.fit_transform(X)

    fig, (axc, axt) = plt.subplots(1, 2, figsize=(8.4, 4.0))
    for a in (axc, axt):  # UMAP coords are arbitrary — no grid
        a.grid(False)

    for label, colour in [("Human", HUMAN), ("Agent", AGENT)]:
        m = (d.is_human.fillna(False) == (label == "Human")).to_numpy()
        axc.scatter(xy[m, 0], xy[m, 1], s=22, c=colour, alpha=0.8,
                    edgecolors="white", linewidths=0.4, label=f"{label} (n={m.sum()})")
    axc.legend(frameon=False, fontsize=7, loc="best")
    _subtitle(axc, "By cohort")
    axc.set_xlabel("UMAP-1")
    axc.set_ylabel("UMAP-2")

    teams = sorted(d.team.dropna().unique())
    cmap = plt.get_cmap("tab20")
    for k, t in enumerate(teams):
        m = (d.team == t).to_numpy()
        axt.scatter(xy[m, 0], xy[m, 1], s=22, color=cmap(k % 20), alpha=0.85,
                    edgecolors="white", linewidths=0.4, label=t)
    # Team legend below the panel, multi-column, so the 16 names aren't clipped.
    axt.legend(frameon=False, fontsize=6, loc="center left",
               bbox_to_anchor=(1.02, 0.5), ncol=1, handletextpad=0.4,
               labelspacing=0.6)
    _subtitle(axt, "By team")
    axt.set_xlabel("UMAP-1")
    axt.set_ylabel("UMAP-2")
    fig.suptitle("ESM-2 embedding UMAP", fontsize=11, fontweight="medium",
                 color=_TITLE_GREY,
                 fontfamily=["GT Pressura Extended", "Geist", "DejaVu Sans"],
                 y=1.02)
    return _save(fig, "fig5_esm2_umap")


def fig_identity_heatmap(df) -> Path:
    """Pairwise sequence-identity heatmap, designs ordered by team (within +
    between team structure visible as on/off-diagonal blocks)."""
    import Levenshtein

    d = df[["design_id", "team", "is_human", "sequence"]].copy()
    d = d.sort_values(["is_human", "team", "design_id"], ascending=[False, True, True])
    d = d.reset_index(drop=True)
    seqs = d.sequence.tolist()
    n = len(seqs)
    M = np.eye(n, dtype=np.float32)
    for i in range(n):
        for j in range(i + 1, n):
            L = max(len(seqs[i]), len(seqs[j]))
            ident = 1.0 - Levenshtein.distance(seqs[i], seqs[j]) / L
            M[i, j] = M[j, i] = ident

    fig, ax = plt.subplots(figsize=(5.6, 5.0))
    im = ax.imshow(M, cmap="mako" if "mako" in plt.colormaps() else "viridis",
                   vmin=0, vmax=1, origin="upper", interpolation="nearest")
    # team block boundaries
    bounds, labels, pos = [], [], []
    start = 0
    for team, g in d.groupby("team", sort=False):
        end = start + len(g)
        bounds.append(end)
        labels.append(team)
        pos.append((start + end) / 2)
        start = end
    for b in bounds[:-1]:
        ax.axhline(b - 0.5, color="white", lw=0.6)
        ax.axvline(b - 0.5, color="white", lw=0.6)
    ax.set_xticks(pos)
    ax.set_xticklabels(labels, rotation=90, fontsize=5.5)
    ax.set_yticks(pos)
    ax.set_yticklabels(labels, fontsize=5.5)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Pairwise identity", fontsize=7)
    cb.ax.tick_params(labelsize=6)
    _title(ax, "Pairwise sequence identity (ordered by cohort, team)")
    return _save(fig, "fig_identity_heatmap")


def fig_foldseek_clustering(df) -> Path | None:
    """Foldseek TM-score matrix, hierarchically ordered, cohort-annotated."""
    cache = HERE / "foldseek_cache.npz"
    if not cache.exists():
        return None
    from scipy.cluster.hierarchy import fcluster, leaves_list, linkage
    from scipy.spatial.distance import squareform

    z = np.load(cache)
    ids = z["design_id"].astype(int)
    tm = z["tm"].astype(float)
    tm = np.clip((tm + tm.T) / 2, 0.0, 1.0)
    np.fill_diagonal(tm, 1.0)
    dist = np.clip(1.0 - tm, 0.0, 1.0)
    np.fill_diagonal(dist, 0.0)
    link = linkage(squareform(dist, checks=False), method="average")
    order = leaves_list(link)
    clusters = fcluster(link, t=0.5, criterion="distance")  # TM > 0.5 ~ same fold
    n_clusters = len(np.unique(clusters))

    cohort = (df.set_index("design_id").loc[ids, "is_human"]
              .fillna(False).map({True: HUMAN, False: AGENT}).to_numpy())

    fig, ax = plt.subplots(figsize=(5.6, 5.2))
    M = tm[np.ix_(order, order)]
    im = ax.imshow(M, cmap="viridis", vmin=0, vmax=1, interpolation="nearest")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("Binder (Foldseek-ordered)")
    ax.set_ylabel("Binder (Foldseek-ordered)")
    # cohort strip along the top
    import matplotlib.colors as mcolors

    strip = ax.inset_axes([0, 1.01, 1, 0.03])
    strip_rgb = np.array([mcolors.to_rgb(c) for c in cohort[order]])[None, :, :]
    strip.imshow(strip_rgb, aspect="auto")
    strip.set_xticks([])
    strip.set_yticks([])
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Foldseek TM-score", fontsize=7)
    cb.ax.tick_params(labelsize=6)
    med = float(np.median(tm[~np.eye(len(ids), dtype=bool)]))
    _title(ax, f"Structural clustering — {n_clusters} folds @TM>0.5, median TM {med:.2f}")
    return _save(fig, "fig_foldseek_clustering")


def fig_epitope(df) -> Path | None:
    """TREM2 contact frequency by cohort + distinct epitope-patch count."""
    cache = HERE / "epitope_cache.npz"
    if not cache.exists():
        return None
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import pdist

    z = np.load(cache)
    ids = z["design_id"].astype(int)
    foot = z["footprints"].astype(bool)
    resid = z["trem2_resid"].astype(int)
    is_h = (df.set_index("design_id").loc[ids, "is_human"].fillna(False)).to_numpy()

    # Epitope patches: cluster footprints by Jaccard distance.
    jac = pdist(foot.astype(float), metric="jaccard")
    link = linkage(jac, method="average")
    patch = fcluster(link, t=0.7, criterion="distance")
    n_h = len(np.unique(patch[is_h]))
    n_a = len(np.unique(patch[~is_h]))

    fig, (axf, axb) = plt.subplots(1, 2, figsize=(8.4, 3.4),
                                   gridspec_kw={"width_ratios": [3, 1]})
    x = np.arange(len(resid))
    for label, mask in [("Human", is_h), ("Agent", ~is_h)]:
        freq = foot[mask].mean(0)
        axf.fill_between(x, freq, color=COHORT_FILL[label], alpha=0.55,
                         step="mid")
        axf.plot(x, freq, color=COHORT_EDGE[label], lw=1.3,
                 drawstyle="steps-mid", label=f"{label} (n={mask.sum()})")
    step = max(1, len(resid) // 10)
    axf.set_xticks(x[::step])
    axf.set_xticklabels(resid[::step], fontsize=6)
    axf.set_xlabel("TREM2 residue (author numbering)")
    axf.set_ylabel("Fraction of cohort contacting")
    axf.legend(frameon=False, fontsize=7, loc="upper right")
    _subtitle(axf, "Contact frequency by residue")

    axb.bar(["Human", "Agent"], [n_h, n_a],
            color=[COHORT_FILL["Human"], COHORT_FILL["Agent"]],
            edgecolor=[COHORT_EDGE["Human"], COHORT_EDGE["Agent"]],
            linewidth=1.0, width=0.6)
    for i, v in enumerate([n_h, n_a]):
        axb.text(i, v + 0.1, str(v), ha="center", va="bottom", fontsize=8)
    axb.set_ylabel("Distinct epitope patches")
    axb.set_ylim(0, max(n_h, n_a) + 1.5)
    _subtitle(axb, "Distinct patches (Jaccard < 0.7)")
    fig.suptitle("TREM2 epitope usage by cohort", fontsize=11,
                 fontweight="medium", color=_TITLE_GREY,
                 fontfamily=["GT Pressura Extended", "Geist", "DejaVu Sans"],
                 y=1.02)
    return _save(fig, "fig_epitope_regions")


def fig_cohort_funnel(df) -> Path | None:
    """The selection funnel + the gap that collapses (reads cohort_funnel json)."""
    sp = HERE / "cohort_funnel_summary.json"
    if not sp.exists():
        return None
    j = json.loads(sp.read_text())
    f = j["funnel"]
    sv, ex, bo = (f["survived_ipsae_cut"], f["expressed_given_screened"],
                  f["bound_given_screened"])
    # cohort cascade counts: submitted -> survived -> expressed -> bound
    human = [sv["human"][1], sv["human"][0], ex["human"][0], bo["human"][0]]
    agent = [sv["agent"][1], sv["agent"][0], ex["agent"][0], bo["agent"][0]]
    stages = ["Submitted", "Survived\nipSAE cut", "Expressed", "Bound"]

    fig, (axf, axg) = plt.subplots(1, 2, figsize=(8.6, 3.8),
                                   gridspec_kw={"width_ratios": [1.25, 1]})
    xs = np.arange(len(stages))
    for label, vals, colour in [("Human", human, HUMAN),
                                ("Agent", agent, AGENT)]:
        axf.plot(xs, vals, "-o", lw=2.0, ms=6, color=colour, label=label)
        for xi, v in zip(xs, vals, strict=True):
            axf.text(xi, v + 2.5, str(v), ha="center", va="bottom",
                     fontsize=7, color=colour)
    axf.set_xticks(xs)
    axf.set_xticklabels(stages, fontsize=7.5)
    axf.set_ylabel("Designs remaining")
    axf.set_ylim(0, max(human) + 14)
    axf.legend(frameon=False, fontsize=7.5, loc="upper right")
    _subtitle(axf, "Selection funnel (count remaining)")

    # The decisive contrast: a big, significant survival gap vs a tied hit rate
    pairs = [
        ("Survived\nipSAE cut", sv["rate_human"], sv["rate_agent"],
         sv["fisher_p_two_sided"]),
        ("Hit rate\n| tested", bo["rate_human"], bo["rate_agent"],
         bo["fisher_p_two_sided"]),
    ]
    gx = np.arange(len(pairs))
    w = 0.36
    axg.bar(gx - w / 2, [p[1] for p in pairs], w, color=COHORT_FILL["Human"],
            edgecolor=COHORT_EDGE["Human"], linewidth=1.0, label="Human")
    axg.bar(gx + w / 2, [p[2] for p in pairs], w, color=COHORT_FILL["Agent"],
            edgecolor=COHORT_EDGE["Agent"], linewidth=1.0, label="Agent")
    for i, (_, rh, ra, pv) in enumerate(pairs):
        axg.text(i - w / 2, rh + 0.02, f"{rh:.0%}", ha="center",
                 va="bottom", fontsize=7)
        axg.text(i + w / 2, ra + 0.02, f"{ra:.0%}", ha="center",
                 va="bottom", fontsize=7)
        tag = f"p={pv:.3f}{'  *' if pv < 0.05 else '  (ns)'}"
        axg.text(i, max(rh, ra) + 0.11, tag, ha="center", va="bottom",
                 fontsize=7, color="#5C6773")
    axg.set_xticks(gx)
    axg.set_xticklabels([p[0] for p in pairs], fontsize=7.5)
    axg.set_ylabel("Rate")
    axg.set_ylim(0, 1.05)
    axg.legend(frameon=False, fontsize=7.5, loc="upper right")
    _subtitle(axg, "The decisive gap collapses")

    fig.suptitle("Humans won the proxy, tied on reality", fontsize=11,
                 fontweight="medium", color=_TITLE_GREY,
                 fontfamily=["GT Pressura Extended", "Geist", "DejaVu Sans"],
                 y=1.03)
    return _save(fig, "fig_cohort_funnel")


def _write_report(df, written: list[Path], esm2_ok: bool) -> None:
    scr = df[df.submitted_to_lab.fillna(False)]
    summary = {
        "n_designs": len(df),
        "n_screened": len(scr),
        "n_figures": len(written),
        "figures": [f"figures/paper/{p.name}" for p in written],
        "esm2_model": "esm2_t33_650M_UR50D" if esm2_ok else None,
        "foldseek_done": (HERE / "foldseek_cache.npz").exists(),
        "epitope_done": (HERE / "epitope_cache.npz").exists(),
        "fig3_status": "excluded by request",
    }
    (HERE / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    md = [
        "# Manuscript figures (analyses/helen)\n",
        "**Headline:** the Notion-spec figures, rendered as single SVGs in "
        "blog styling (Agent `#30C5F5`, Human `#142933`).\n",
        "## Method\n",
        "`load_designs()` → matplotlib (`apply_theme()` + blog tokens), one SVG "
        "per figure. ESM-2 (`esm2_t33_650M_UR50D`), Foldseek TM (local static "
        "binary, ESMFold CIFs), and epitope footprints (gemmi, Boltz-2 CIFs, "
        "5 Å contacts) are precomputed into `*_cache.npz` by the sibling "
        "scripts; figure code is deterministic and reads the caches.\n",
        "## Results\n",
    ]
    md += [f"- `figures/paper/{p.name}`" for p in written]
    md += [
        "\n## Caveats\n",
        "- **Fonts** are embedded per-SVG as base64 @font-face (blog Geist + "
        "GT Pressura Extended woff2); self-contained, overrides the "
        "matplotlibrc Helvetica/DejaVu policy by design. GT Pressura is a "
        "paid face — embedding in distributed artifacts is a licensing call.",
        "- **Fig 3** (BLI sensorgrams): excluded by request.",
        "- **ESM-2 pLL** is a length-normalized single-pass log-likelihood; true "
        "masked-marginal is infeasible on CPU/MPS at 650M. Document if revised.",
        "- **Fig 4 taxonomy** (design_method_normalized → 5-category, Cotet 2025) "
        "is an interpretation; PXDesign and BindCraft bucket as De novo "
        "(BindCraft per the Notion draft's own grouping).",
        "- **Foldseek/epitope** are screened-only (100/141); the 41 non-screened "
        "designs have no ProteinBase CIF.",
        "- Per-cohort epitope/structure counts inherit the top-100 selection "
        "imbalance (65 human / 35 agent).",
        f"- ESM-2 figures present: {esm2_ok} (cache "
        f"{'found' if esm2_ok else 'absent — rerun esm2_compute.py'}).",
    ]
    (HERE / "report.md").write_text("\n".join(md) + "\n")


def main() -> None:
    _style()
    df = load_designs()
    written = [
        fig1_ipsae(df),
        fig2_expression_binding(df),
        fig4_design_methods(df),
        fig6_sequence_length(df),
        fig7_metric_roc(df),
        fig_identity_heatmap(df),
    ]
    for builder in (fig_foldseek_clustering, fig_epitope, fig5_esm2_umap,
                    fig_cohort_funnel):
        out = builder(df)
        if out is not None:
            written.append(out)

    _write_report(df, written, esm2_ok=(HERE / "esm2_cache.npz").exists())
    print(f"[helen] wrote {len(written)} SVGs to figures/paper/")


if __name__ == "__main__":
    main()
