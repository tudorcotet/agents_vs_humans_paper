"""Superpose the three crystal TREM2 chains onto the design_017 TREM2 frame (rigid CA fit on shared
UniProt residues), move each binding PARTNER into that frame, and write combined PDBs for the Fig-8
hero panels. Because every panel shares design_017's chain A, TREM2 has the SAME orientation in all four.

Partners:  6YYE -> scFv-2 (protein),  6Y6C -> scFv-4 (protein),  6B8O -> phosphatidylserine (HETATM).
Only the partner copy actually contacting TREM2 (any atom <= 5 A) is kept, so a 2nd crystallographic copy
doesn't float in.  Outputs: pdbs/fig8_{scfv2,scfv4,ps}.pdb  (chain A = TREM2, partner = chains P/Q/L).
Run:  uv run python figures/blender/superpose_crystals.py
"""
import warnings, numpy as np
from pathlib import Path
from Bio.PDB import MMCIFParser, PDBParser, Superimposer, PDBIO, NeighborSearch, Select
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
warnings.filterwarnings("ignore")

HERE = Path(__file__).resolve().parent
REFS = HERE.parents[1] / "refs"
PDBS = HERE / "pdbs"
DESIGN = PDBS / "design_017.pdb"
UOFF = 17                       # UniProt = construct + 17
IGSF = ("HNTTVFQGVAGQSLQVSCPYDSMKHWGRRKAWCRQLGEKGPCQRVVSTHNLWLLSFLRRWNGSTAITDDTLGG"
        "TLTITLRNLQPHDAGLYQCQSLHGSEADTLRKVLVEVLADPL")

def one(res):
    from Bio.PDB.Polypeptide import three_to_index, index_to_one
    try: return index_to_one(three_to_index(res.resname))
    except Exception: return "X"

def chain_seq(ch):
    return "".join(one(r) for r in ch if r.id[0] == " ")

def igsf_score(ch):
    s = chain_seq(ch)
    return sum(1 for i in range(0, len(IGSF) - 15, 5) if IGSF[i:i + 15] in s)

def find_trem2_chain(model):
    best, best_score = None, 0
    for ch in model.get_chains():
        score = igsf_score(ch)
        if score > best_score:
            best, best_score = ch, score
    return best

# --- design_017 TREM2 (chain A): CA by UniProt residue ---
dmodel = PDBParser(QUIET=True).get_structure("d017", str(DESIGN))[0]
dA = dmodel["A"]
design_ca = {r.id[1] + UOFF: r["CA"] for r in dA if r.id[0] == " " and "CA" in r}   # UniProt -> Atom
design_atoms = [a for r in dA if r.id[0] == " " for a in r if a.element != "H"]

def process(pdb, kind, partner_chain_label):
    st = MMCIFParser(QUIET=True).get_structure(pdb, str(REFS / f"{pdb}.cif"))
    model = st[0]
    t2 = find_trem2_chain(model)
    trem2_ids = {ch.id for ch in model.get_chains() if igsf_score(ch) > 0}   # every TREM2 copy
    cryst_ca = {r.id[1]: r["CA"] for r in t2 if r.id[0] == " " and "CA" in r}         # UniProt -> Atom
    shared = sorted(set(design_ca) & set(cryst_ca))
    fixed = [design_ca[u] for u in shared]
    moving = [cryst_ca[u] for u in shared]
    sup = Superimposer()
    sup.set_atoms(fixed, moving)                    # maps crystal -> design frame
    # collect partner atoms (protein chains != TREM2, or HETATM ligand), transform, keep contacting copy
    partner_atoms = []
    for ch in model.get_chains():
        for r in ch:
            if kind == "protein":
                if ch.id in trem2_ids or r.id[0] != " ":   # exclude every TREM2 copy
                    continue
                partner_atoms += [a for a in r if a.element != "H"]
            else:  # hetatm ligand
                if r.id[0].strip().startswith("H") and r.resname != "HOH":
                    partner_atoms += [a for a in r if a.element != "H"]
    sup.apply(partner_atoms)                        # now in design frame
    ns = NeighborSearch(design_atoms)
    from collections import Counter
    if kind == "protein":
        # scFv-2/-4 are SINGLE-CHAIN Fvs; crystals may hold 2 copies -> keep only the best-contacting chain
        contact = Counter()
        for a in partner_atoms:
            if ns.search(a.coord, 5.0, level="A"):
                contact[a.get_parent().get_parent().id] += 1
        best = contact.most_common(1)[0][0]
        keep_res = [r for r in model[best] if r.id[0] == " "]
        print(f"   partner chain contacts: {dict(contact)} -> keeping scFv chain {best}")
    else:  # ligand: keep only the phosphatidylserine (PSF), not crystallization additives
        keep_res = []
        for a in partner_atoms:
            r = a.get_parent()
            if r.resname == "PSF" and r not in keep_res and ns.search(a.coord, 6.0, level="A"):
                keep_res.append(r)
        print(f"   PS ligand: kept {len(keep_res)} PSF residue(s)")
    keep_atoms = [a for r in keep_res for a in r if a.element != "H"]
    # geometry report: partner centroid direction vs domain centre
    dc = np.mean([a.coord for a in design_atoms], axis=0)
    pc = np.mean([a.coord for a in keep_atoms], axis=0)
    binder_dir = (pc - dc) / np.linalg.norm(pc - dc)
    print(f"{pdb} [{kind}]: TREM2 chain {t2.id}, {len(shared)} shared CA, fit RMSD {sup.rms:.2f} A | "
          f"partner {len(keep_res)} res / {len(keep_atoms)} atoms kept | dir {binder_dir.round(2)}")

    # --- write combined PDB: chain A = design TREM2, each partner source-chain kept separate ---
    from Bio.PDB.Chain import Chain
    out_struct = Structure("x"); out_model = Model(0); out_struct.add(out_model)
    ca_chain = Chain("A")
    for r in dA:
        if r.id[0] == " ":
            ca_chain.add(r.copy())
    out_model.add(ca_chain)
    by_src = {}
    for r in keep_res:
        by_src.setdefault(r.get_parent().id, []).append(r)
    labels = iter([c for c in "PQRSTUV" if c != "A"])
    for src, residues in by_src.items():
        lab = partner_chain_label if len(by_src) == 1 else next(labels)
        ch = Chain(lab)
        for i, r in enumerate(sorted(residues, key=lambda r: r.id[1])):
            nr = r.copy(); nr.id = (r.id[0], i + 1, " ")       # renumber; keep hetflag
            try: ch.add(nr)
            except Exception: pass
        out_model.add(ch)
    io = PDBIO(); io.set_structure(out_struct)
    out = PDBS / f"fig8_{pdb.lower()}.pdb"
    io.save(str(out))
    return out

print("=== superposing crystals onto design_017 TREM2 frame ===")
process("6YYE", "protein", "P")   # scFv-2
process("6Y6C", "protein", "P")   # scFv-4
process("6B8O", "hetatm", "L")    # phosphatidylserine
print("done -> pdbs/fig8_6yye.pdb, fig8_6y6c.pdb, fig8_6b8o.pdb")
