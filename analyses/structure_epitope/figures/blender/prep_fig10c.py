"""Fig 10c prep — build one combined PDB: design_017 TREM2 (chain A) + every hit binder's chain B,
all already superposed on the shared frame. Cohort is baked into the B-factor (human=1, agent=2, TREM2=0)
so the Blender renderer colours the binder fan by cohort without chain-order bookkeeping.
Run:  uv run python figures/blender/prep_fig10c.py  ->  pdbs/fig10c_superposed.pdb
"""
import csv, warnings
from pathlib import Path
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
warnings.filterwarnings("ignore")

HERE = Path(__file__).resolve().parent
PDBS = HERE / "pdbs"
manifest = {int(r["design_id"]): r for r in csv.DictReader(open(HERE / "manifest.csv"))}
hits = [d for d, r in manifest.items() if r["is_hit"] == "True" and (PDBS / f"design_{d:03d}.pdb").exists()]
hits.sort()

parser = PDBParser(QUIET=True)
out = Structure("x"); model = Model(0); out.add(model)

# chain A = design_017 TREM2, B-factor 0
d17 = parser.get_structure("d17", str(PDBS / "design_017.pdb"))[0]
A = Chain("A")
for r in d17["A"]:
    if r.id[0] == " ":
        nr = r.copy()
        for a in nr: a.bfactor = 0.0
        A.add(nr)
model.add(A)

pool = [c for c in "BCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789" if c != "A"]
labels = iter(pool)
nh = na = 0
for did in hits:
    st = parser.get_structure(str(did), str(PDBS / f"design_{did:03d}.pdb"))[0]
    if "B" not in st: continue
    cohort = manifest[did]["cohort"].strip().lower()
    bf = 1.0 if cohort == "human" else 2.0
    nh += cohort == "human"; na += cohort == "agent"
    ch = Chain(next(labels))
    i = 0
    for r in st["B"]:
        if r.id[0] != " ": continue
        i += 1
        nr = r.copy(); nr.id = (" ", i, " ")
        for a in nr: a.bfactor = bf
        try: ch.add(nr)
        except Exception: pass
    model.add(ch)

io = PDBIO(); io.set_structure(out)
outpath = PDBS / "fig10c_superposed.pdb"
io.save(str(outpath))
print(f"wrote {outpath.name}: TREM2 + {len(hits)} binders ({nh} human, {na} agent)")
