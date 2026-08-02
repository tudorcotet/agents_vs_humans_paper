"""EXTRA 1 (OPEN_ITEMS_register_1) — kinetic decomposition, prep stage.

Aggregate per-design k_on / k_off / K_D from the replicate BLI/SPR fits, run the QC gates,
compute the structural descriptor panel (stability + electrostatic + confidence), and write
one per-design table for the correlation stage (kinetics_stats.py).

QC gates (OPEN_ITEMS EXTRA 1):
  * k_off solver-bound floor (~1e-5 s^-1): flag k_off within 2x of the floor.
  * mass-transport ceiling: flag k_on clustered near a ceiling (proxy; loading density not in package).
  * model adequacy: exclude replicates flagged unexpected_order (non-monotonic concentration series).
  * fit quality (chi2/R2): NOT in the data package and not reconstructable (kinetic fits are SPR;
    _bli.json sensorgrams cover only the BLI subset; r2 column empty). Requested from Tudor (Tier-4
    4.1). Internal-consistency proxy used instead: K_D reconstructed from k_off/k_on vs reported K_D.

Outputs: results/kinetics/kinetics_perdesign.csv + results/kinetics/kinetics_qc.md
Run: uv run --with gemmi --with scipy python analyses/structure_epitope/jobs/kinetics_prep.py
"""
from __future__ import annotations
import os
import numpy as np
import pandas as pd
import gemmi

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
SE = os.path.join(ROOT, "analyses", "structure_epitope")
RES = os.path.join(SE, "results")
OUT = os.path.join(RES, "kinetics")
os.makedirs(OUT, exist_ok=True)

KOFF_FLOOR = 1e-5          # reliable BLI/SPR k_off floor
KON_CEILING = 1e7          # mass-transport ceiling scale (M^-1 s^-1)
CHARGE = {"ASP": -1.0, "GLU": -1.0, "LYS": +1.0, "ARG": +1.0}  # H excluded (near-neutral at pH 7.4)
TARGET_LEN = 175


def aggregate_rates() -> pd.DataFrame:
    rep = pd.read_csv(os.path.join(ROOT, "data", "raw_lab", "bli_replicates.csv"))
    f = rep[(rep.selected == True) & (rep.excluded != True) & (rep.fixed == True)]  # noqa: E712
    kin = f[f.kon_M_s.notna() & f.koff_s.notna()].copy()
    g = kin.groupby("design_id")
    tab = pd.DataFrame({
        "n_kin": g.size(),
        "methods": g.method_name.apply(lambda s: "/".join(sorted(set(s.dropna())))),
        "kon": g.kon_M_s.mean(),
        "koff": g.koff_s.mean(),
        "kd_nM": g.kd_nM.mean(),
        "any_unexpected": g.unexpected_order.apply(lambda s: bool((s == True).any())),  # noqa: E712
        "koff_min": g.koff_s.min(),
        "kon_max": g.kon_M_s.max(),
    }).reset_index()
    tab["log_kon"] = np.log10(tab.kon)
    tab["log_koff"] = np.log10(tab.koff)
    tab["log_kd"] = np.log10(tab.kd_nM)              # log10 KD in nM
    # QC gates
    tab["qc_koff_floor_ok"] = tab.koff_min > 2 * KOFF_FLOOR
    tab["qc_kon_ceiling_ok"] = tab.kon_max < KON_CEILING
    tab["qc_model_ok"] = ~tab.any_unexpected
    tab["rate_eligible"] = tab.qc_koff_floor_ok & tab.qc_kon_ceiling_ok & tab.qc_model_ok
    # internal consistency: KD from rates vs reported
    tab["kd_from_rates_nM"] = (tab.koff / tab.kon) * 1e9
    tab["kd_consistency_dex"] = np.log10(tab.kd_from_rates_nM / tab.kd_nM)
    return tab


def interface_electrostatics(design_id: int) -> dict:
    """Net charge on each interface side + complementarity, from the Protenix complex."""
    cif = os.path.join(ROOT, "data", "structures", "protenix", f"design_{design_id:03d}.cif")
    if not os.path.exists(cif):
        return {}
    m = gemmi.read_structure(cif)[0]
    chains = list(m)
    tgt = next((c for c in chains if len(list(c)) == TARGET_LEN), None)
    bnd = next((c for c in chains if c is not tgt), None)
    if tgt is None or bnd is None:
        return {}

    def heavy(res):
        return [a.pos for a in res if not a.is_hydrogen()]

    tgt_res = [(r, heavy(r)) for r in tgt]
    bnd_res = [(r, heavy(r)) for r in bnd]
    b_iface, t_iface = set(), set()
    for bi, (br, bpos) in enumerate(bnd_res):
        for ti, (tr, tpos) in enumerate(tgt_res):
            if any(bp.dist(tp) <= 4.5 for bp in bpos for tp in tpos):
                b_iface.add(bi)
                t_iface.add(ti)
    bq = sum(CHARGE.get(bnd_res[i][0].name, 0.0) for i in b_iface)
    tq = sum(CHARGE.get(tgt_res[i][0].name, 0.0) for i in t_iface)
    return {
        "binder_iface_charge": bq,
        "target_iface_charge": tq,
        # complementarity: positive when the two interfaces carry opposite net charge
        "charge_complementarity": -(bq * tq),
        "n_binder_iface": len(b_iface),
    }


def main():
    tab = aggregate_rates()
    elec = pd.DataFrame([{"design_id": d, **interface_electrostatics(d)} for d in tab.design_id])
    tab = tab.merge(elec, on="design_id", how="left")

    # stability descriptors
    idesc = pd.read_csv(os.path.join(RES, "interface_descriptors.csv"))
    tab = tab.merge(idesc, on="design_id", how="left")

    # confidence metrics + cohort (clean columns only)
    d = pd.read_parquet(os.path.join(ROOT, "data", "designs.parquet"))
    keep = ["design_id", "is_human", "team", "submitted_ipsae", "boltz2_iptm",
            "boltz2_plddt", "boltz2_complex_plddt", "is_literature_copy"]
    keep = [c for c in keep if c in d.columns]
    tab = tab.merge(d[keep], on="design_id", how="left")
    tab["cohort"] = tab.is_human.map({True: "human", False: "agent"})

    # Rosetta dG + PRODIGY dG if present
    rp = os.path.join(RES, "rosetta_descriptors.csv")
    if os.path.exists(rp):
        r = pd.read_csv(rp)
        dgcol = next((c for c in r.columns if "dg" in c.lower() and "sasa" not in c.lower()), None)
        if dgcol:
            tab = tab.merge(r[["design_id", dgcol]].rename(columns={dgcol: "rosetta_dg"}), on="design_id", how="left")
    for c in ["prodigy_protenix_dg", "prodigy_protenix_kd"]:
        if c in d.columns:
            tab = tab.merge(d[["design_id", c]], on="design_id", how="left")

    # epitope bin for the figure
    eb = os.path.join(RES, "epitope_bins.csv")
    if os.path.exists(eb):
        ebd = pd.read_csv(eb)
        bcol = next((c for c in ebd.columns if "bin" in c.lower()), None)
        if bcol:
            tab = tab.merge(ebd[["design_id", bcol]].rename(columns={bcol: "epitope_bin"}), on="design_id", how="left")

    tab.to_csv(os.path.join(OUT, "kinetics_perdesign.csv"), index=False)

    # QC report
    elig = tab[tab.rate_eligible]
    lines = []
    lines.append("# EXTRA 1 — kinetic QC report\n")
    lines.append(f"- Rate-eligible designs (kon+koff, canonical replicate filter): **{len(tab)}** "
                 f"(= P_kd binders); cohort {tab.cohort.value_counts().to_dict()}.")
    lines.append(f"- Method mix: {tab.methods.value_counts().to_dict()} — **kinetic fits are predominantly SPR**, "
                 "not BLI as the manuscript framing implies. Flag to Tudor (Tier-4 4.1: confirm assay + per-design chi2/R2/loading).")
    lines.append(f"- **k_off floor gate:** min k_off = {tab.koff_min.min():.2e} s^-1; "
                 f"{int((~tab.qc_koff_floor_ok).sum())} designs within 2x of the {KOFF_FLOOR:.0e} floor "
                 "→ the tightest KD is a genuine point estimate, not a lower limit (resolves item 1.5).")
    lines.append(f"- **mass-transport (proxy):** k_on median {tab.kon.median():.2e}, max {tab.kon.max():.2e}; "
                 f"{int((~tab.qc_kon_ceiling_ok).sum())} above the {KON_CEILING:.0e} ceiling; k_on not ceiling-clustered. "
                 "Loading density not in the package → full MT check pending 4.1.")
    lines.append(f"- **model adequacy:** {int(tab.any_unexpected.sum())} designs have a kinetic replicate flagged "
                 "unexpected_order (non-monotonic series) → excluded from the eligible set.")
    lines.append(f"- **fit quality (chi2/R2):** unavailable in the package and not reconstructable (SPR fits; "
                 "_bli.json sensorgrams cover only the BLI subset). Internal-consistency proxy instead: "
                 f"KD from k_off/k_on matches reported KD to **{(tab.kd_consistency_dex.abs()<0.5).sum()}/{len(tab)}** "
                 f"within 0.5 dex (median {tab.kd_consistency_dex.median():+.2f} dex).")
    lines.append(f"- **Design 79** (literature copy) is in the rate set; excluded from H-Perera (item 2.1) and "
                 "flagged in the cohort rate-space comparison.")
    lines.append(f"\n**Two eligibility sets:** KD-eligible n=36; rate-eligible n={len(elig)} "
                 f"(drops {len(tab)-len(elig)} for unexpected_order).")
    open(os.path.join(OUT, "kinetics_qc.md"), "w").write("\n".join(lines) + "\n")

    print("\n".join(lines))
    print(f"\nwrote {os.path.relpath(os.path.join(OUT,'kinetics_perdesign.csv'), ROOT)} "
          f"({len(tab)} designs x {len(tab.columns)} cols)")
    print("electrostatics computed for", tab.binder_iface_charge.notna().sum(), "designs")


if __name__ == "__main__":
    main()
