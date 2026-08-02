"""
Extract crystallographic TREM2 epitopes from reference PDBs and compute per-design footprint
overlap. Reference TREM2 chains use UniProt numbering; our design footprints are already in UniProt
(construct+17), so overlap is directly comparable.
"""
import warnings, numpy as np, pandas as pd
from pathlib import Path
from Bio.PDB import MMCIFParser, NeighborSearch
warnings.filterwarnings("ignore")
ROOT = Path("/Users/amyhe/Desktop/trem2_2026/agents_vs_humans_paper")
REFS = ROOT / "analyses/structure_epitope/refs"
RES = ROOT / "analyses/structure_epitope/results"
IGSF = "HNTTVFQGVAGQSLQVSCPYDSMKHWGRRKAWCRQLGEKGPCQRVVSTHNLWLLSFLRRWNGSTAITDDTLGGTLTITLRNLQPHDAGLYQCQSLHGSEADTLRKVLVEVLADPL"

def one(res):
    from Bio.PDB.Polypeptide import three_to_index, index_to_one
    try: return index_to_one(three_to_index(res.resname))
    except Exception: return "X"

def chain_seq(ch):
    return "".join(one(r) for r in ch if r.id[0] == " ")

def find_trem2_chain(model):
    """chain whose sequence best contains/overlaps the IgSF sequence."""
    best = None; best_score = 0
    for ch in model.get_chains():
        s = chain_seq(ch)
        # crude overlap: longest common 15-mer presence
        score = sum(1 for i in range(0, len(IGSF)-15, 5) if IGSF[i:i+15] in s)
        if score > best_score: best_score = score; best = ch
    return best, best_score

def epitope(pdb, partner_is_hetatm=False):
    st = MMCIFParser(QUIET=True).get_structure(pdb, str(REFS / f"{pdb}.cif"))
    model = st[0]
    t2, sc = find_trem2_chain(model)
    if t2 is None or sc == 0: return None, f"no TREM2 chain (score {sc})"
    # sample TREM2 numbering
    nums = [r.id[1] for r in t2 if r.id[0] == " "]
    info = f"TREM2 chain {t2.id}, {len(nums)} res, resnum {min(nums)}..{max(nums)}"
    # partner atoms = everything not in the TREM2 chain (protein partners), or HETATM ligands
    partner_atoms = []
    for ch in model.get_chains():
        for r in ch:
            is_het = r.id[0].strip() not in ("", "W") and r.id[0].startswith("H")
            if partner_is_hetatm:
                if r.id[0].strip().startswith("H") and r.resname not in ("HOH",):
                    partner_atoms += [a for a in r if a.element != "H"]
            else:
                if ch.id == t2.id: continue
                if r.id[0] == " ":
                    partner_atoms += [a for a in r if a.element != "H"]
    if not partner_atoms: return None, info + " | no partner atoms"
    ns = NeighborSearch(partner_atoms)
    epi = set()
    for r in t2:
        if r.id[0] != " ": continue
        for a in r:
            if a.element == "H": continue
            if ns.search(a.coord, 4.5, level="A"):
                epi.add(r.id[1]); break
    return sorted(epi), info

# ---- extract reference epitopes (UniProt numbering) ----
refs = {}
print("=== reference epitope extraction ===")
for pdb, het in [("6YYE", False), ("6Y6C", False), ("6B8O", True), ("5ELI", False)]:
    epi, info = epitope(pdb, partner_is_hetatm=het)
    print(f"{pdb}: {info}")
    if epi is not None:
        # keep only IgSF-range residues (19..132) to avoid stalk/partner artifacts
        epi_igsf = [e for e in epi if 19 <= e <= 132]
        refs[pdb] = set(epi_igsf)
        print(f"   epitope (UniProt, IgSF-range, n={len(epi_igsf)}): {epi_igsf}")

# name them
REFSETS = {}
if "6YYE" in refs: REFSETS["scFv2_6YYE"] = refs["6YYE"]
if "6Y6C" in refs: REFSETS["scFv4_6Y6C"] = refs["6Y6C"]
if "6B8O" in refs: REFSETS["PS_site_6B8O"] = refs["6B8O"]

# ---- per-design overlap with reference epitopes ----
cons = pd.read_parquet(RES / "footprints_consensus.parquet")
meta = pd.read_parquet(ROOT / "analyses/structure_epitope/data/designs_canonical.parquet")
byd = cons[(cons.consensus == True) & (cons.region == "IgSF")].groupby("design_id")["uniprot"].apply(set).to_dict()
rows = []
for did in meta.design_id:
    fp = byd.get(did, set()); row = {"design_id": did}
    for name, rs in REFSETS.items():
        inter = len(fp & rs); uni = len(fp | rs)
        row[f"jaccard_{name}"] = round(inter/uni, 3) if uni else 0
        row[f"nover_{name}"] = inter
        row[f"frac_{name}"] = round(inter/len(rs), 3) if rs else 0
    rows.append(row)
ov = pd.DataFrame(rows).merge(meta[["design_id","cohort","is_hit","p_expressed","kd_arith_mean_nM_all"]], on="design_id")
ov.to_csv(RES / "epitope_reference_crystal_overlap.csv", index=False)
# save reference sets
pd.DataFrame([(n, sorted(s)) for n, s in REFSETS.items()], columns=["reference","uniprot_residues"]).to_csv(
    RES / "reference_epitope_sets.csv", index=False)

print("\n=== per-design overlap with crystallographic epitopes (expressed; medians) ===")
exp = ov[ov.p_expressed == True]
for name in REFSETS:
    b = exp[exp.is_hit==True][f"frac_{name}"]; nb = exp[exp.is_hit==False][f"frac_{name}"]
    from scipy import stats
    U,p = stats.mannwhitneyu(b, nb, alternative="two-sided")
    print(f"  {name:16s}: binder frac_of_epitope med={b.median():.2f} vs non {nb.median():.2f}  MWU p={p:.3f}")
print("\nwrote epitope_reference_crystal_overlap.csv, reference_epitope_sets.csv")
