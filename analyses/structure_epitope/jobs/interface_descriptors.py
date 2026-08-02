"""
Phase 2: interface descriptor panel (Biopython-only core) on the complex predictors.
Per (design, predictor): buried surface area, interface size, interface chemistry, geometric
H-bonds / salt bridges. Then binder-monomer RMSD (ESMFold monomer vs binder-in-complex) for the
Bennett Type I/II split.

Predictors: Protenix and Chai (both 141 designs, clean) plus ESMFold2 (100 screened designs,
team-contributed complex panel). **Protenix is the section's PRIMARY structural model**, so the
primary file (interface_descriptors.csv) is Protenix-only — every main-text consumer (§2 phase2_stats,
§5 phase6_correlations, §6 kinetics, Figs 9/12) inherits Protenix. The per-predictor long file carries
all three (PREDS) for the SI cross-predictor comparison (jobs/predictor_panel.py), which anchors every
comparison to Protenix.

Descriptors labelled by route in Methods: ShrakeRupley SASA (Biopython, BSD; NOT FreeSASA),
geometric H-bond/salt-bridge (direct implementation, brief §5 fallback).

Output: results/interface_descriptors.csv (Protenix-only, primary) +
        results/interface_descriptors_perpredictor.csv (long, all three predictors).
"""
import warnings, numpy as np, pandas as pd
from pathlib import Path
from Bio.PDB import MMCIFParser, PDBParser, NeighborSearch, Superimposer
from Bio.PDB.SASA import ShrakeRupley
warnings.filterwarnings("ignore")
ROOT = Path("/Users/amyhe/Desktop/trem2_2026/agents_vs_humans_paper")
RES = ROOT / "analyses/structure_epitope/results"
PREDS = ["protenix", "chai", "esmfold2", "boltz2", "af2m"]  # per-predictor long (boltz2/af2m = regenerated, PR #6; af2m is PDB)
PRIMARY_PRED = "protenix"                  # section-wide primary; interface_descriptors.csv = this only
TARGET_LEN = 175
HYDROPHOBIC = set("AVLIMFWPC"); POLAR = set("STNQGYH"); POS = set("KR"); NEG = set("DE")
DON_ACC = {"N", "O"}  # crude H-bond: any N/O donor-acceptor pair <=3.5A

def chains(model):
    cs = list(model.get_chains())
    tgt = next((c for c in cs if len([r for r in c if r.id[0]==" "]) == TARGET_LEN), None)
    if tgt is None: tgt = max(cs, key=lambda c: len([r for r in c if r.id[0]==" "]))
    bnd = next(c for c in cs if c is not tgt)
    return tgt, bnd

def resname1(res):
    from Bio.PDB.Polypeptide import three_to_index, index_to_one
    try: return index_to_one(three_to_index(res.resname))
    except Exception: return "X"

def sasa_sum(entity):
    return sum(a.sasa for a in entity.get_atoms())

def complex_model(predictor, design_id):
    """Load a predicted complex; AF2M ships as PDB, the rest as mmCIF."""
    if predictor == "af2m":
        p = ROOT / f"data/structures/{predictor}/design_{design_id:03d}.pdb"
        parser = PDBParser(QUIET=True)
    else:
        p = ROOT / f"data/structures/{predictor}/design_{design_id:03d}.cif"
        parser = MMCIFParser(QUIET=True)
    if not p.exists(): return None
    return parser.get_structure("m", str(p))[0]

def descriptors(design_id, predictor, sr):
    model = complex_model(predictor, design_id)
    if model is None: return None
    tgt, bnd = chains(model)
    tgt_res = [r for r in tgt if r.id[0]==" "]; bnd_res = [r for r in bnd if r.id[0]==" "]
    # --- BSA via dSASA: SASA(A_alone)+SASA(B_alone)-SASA(complex), all over full complex atoms ---
    sr.compute(model, level="A"); complex_sasa = sasa_sum(bnd) + sasa_sum(tgt)
    # target alone
    bnd_copy_id = bnd.id; model.detach_child(bnd_copy_id); sr.compute(model, level="A")
    tgt_alone = sasa_sum(tgt); model.add(bnd)
    # binder alone
    tgt_copy_id = tgt.id; model.detach_child(tgt_copy_id); sr.compute(model, level="A")
    bnd_alone = sasa_sum(bnd); model.add(tgt)
    bsa = tgt_alone + bnd_alone - complex_sasa   # total buried surface area
    # --- interface residues (heavy-atom 4.5A) ---
    tgt_heavy = [a for a in tgt.get_atoms() if a.element != "H"]
    bnd_heavy = [a for a in bnd.get_atoms() if a.element != "H"]
    ns_t = NeighborSearch(tgt_heavy)
    iface_b = set(); iface_t = set()
    hbonds = 0; salt = 0
    for a in bnd_heavy:
        near = ns_t.search(a.coord, 4.5, level="A")
        for x in near:
            iface_b.add(a.get_parent().id[1]); iface_t.add(x.get_parent().id[1])
        # H-bond-ish: N/O pair <=3.5
        for x in ns_t.search(a.coord, 3.5, level="A"):
            if a.element in DON_ACC and x.element in DON_ACC: hbonds += 1
    # salt bridges: binder charged sidechain atom N/O within 4.0 of target opposite charge
    for rb in bnd_res:
        cb = resname1(rb)
        if cb not in POS | NEG: continue
        for ab in rb.get_atoms():
            if ab.element not in ("N", "O"): continue
            for x in ns_t.search(ab.coord, 4.0, level="A"):
                rt = x.get_parent(); ct = resname1(rt)
                if (cb in POS and ct in NEG) or (cb in NEG and ct in POS):
                    salt += 1; break
    # --- binder interface chemistry (by interface binder residue type) ---
    ib_res = [r for r in bnd_res if r.id[1] in iface_b]
    types = [resname1(r) for r in ib_res]
    n = len(types) or 1
    frac_hyd = sum(t in HYDROPHOBIC for t in types)/n
    frac_pol = sum(t in POLAR for t in types)/n
    frac_chg = sum(t in POS|NEG for t in types)/n
    return dict(design_id=design_id, predictor=predictor, bsa=round(bsa,1),
                n_iface_binder=len(iface_b), n_iface_target=len(iface_t),
                hbonds=hbonds, salt_bridges=salt,
                iface_frac_hydrophobic=round(frac_hyd,3), iface_frac_polar=round(frac_pol,3),
                iface_frac_charged=round(frac_chg,3),
                hbond_density_per100=round(100*hbonds/bsa,3) if bsa>0 else np.nan)

def monomer_rmsd(design_id, predictor):
    """RMSD of binder chain (in complex) vs ESMFold monomer, CA, after superposition. Type I signal."""
    mono = ROOT / f"data/structures/esmfold/design_{design_id:03d}.cif"
    if not mono.exists():
        mono = ROOT / f"data/structures/proteintyper/design_{design_id:03d}.cif"
        if not mono.exists(): return np.nan
    try:
        m1 = complex_model(predictor, design_id)
        if m1 is None: return np.nan
        _, bnd = chains(m1)
        m2 = MMCIFParser(QUIET=True).get_structure("mo", str(mono))[0]
        mono_ch = max(m2.get_chains(), key=lambda c: len([r for r in c if r.id[0]==" "]))
        ca1 = {r.id[1]: r["CA"] for r in bnd if r.id[0]==" " and "CA" in r}
        ca2 = {r.id[1]: r["CA"] for r in mono_ch if r.id[0]==" " and "CA" in r}
        common = sorted(set(ca1) & set(ca2))
        if len(common) < 5: return np.nan
        si = Superimposer(); si.set_atoms([ca1[i] for i in common], [ca2[i] for i in common])
        return round(si.rms, 3)
    except Exception:
        return np.nan

def main():
    df = pd.read_parquet(ROOT / "data/designs.parquet")
    sr = ShrakeRupley()
    rows = []
    import time; t0 = time.time()
    for i, did in enumerate(df.design_id):
        for pred in PREDS:
            d = descriptors(did, pred, sr)
            if d: d["monomer_rmsd"] = monomer_rmsd(did, pred); rows.append(d)
        if (i+1) % 20 == 0:
            import sys; print(f"...{i+1}/141 {time.time()-t0:.0f}s", file=sys.stderr)
    long = pd.DataFrame(rows)
    long.to_csv(RES / "interface_descriptors_perpredictor.csv", index=False)
    # predictor-averaged file: Protenix+Chai ONLY (AVG_PREDS) so downstream consumers are unchanged;
    # ESMFold2 lives only in the per-predictor long file and the SI panel (jobs/predictor_panel.py).
    # Primary interface_descriptors.csv = PROTENIX only (section-wide primary structural model,
    # used by §2/§5/§6 and Figs 9/12). Chai + ESMFold2 live in the per-predictor file for the SI panel.
    primary = long[long.predictor == PRIMARY_PRED].drop(columns=["predictor"]).reset_index(drop=True)
    primary.to_csv(RES / "interface_descriptors.csv", index=False)
    print(f"wrote interface_descriptors.csv rows={len(avg)} time={time.time()-t0:.0f}s")
    print(avg.describe().round(2).to_string())

if __name__ == "__main__":
    main()
