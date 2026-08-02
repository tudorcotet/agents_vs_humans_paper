"""Build the three Fig-8 epitope residue sets on the TREM2 construct frame and verify the
'opposite faces' geometry, so the figure's colouring + camera are backed by coordinates.

Sets (UniProt residues -> construct res_id via construct = UniProt - 17; chain A = design_017 TREM2):
  - design consensus epitope : epitope_per_residue_freq.csv, freq_binder >= 0.5
  - scFv antibody epitope     : union of scFv2_6YYE + scFv4_6Y6C (reference_epitope_sets.csv)
  - PS ligand site            : PS_site_6B8O
Emits fig8_epitope_sets.json for the Blender render script and prints centroid opposition angles.
Run:  uv run python figures/blender/prep_fig8_epitope_sets.py
"""
import csv, json, ast
import numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent
SEC = HERE.parents[1]
PDB = HERE / "pdbs" / "design_017.pdb"
REFSETS = SEC / "results" / "reference_epitope_sets.csv"
FREQ = SEC / "results" / "epitope_per_residue_freq.csv"
OUT = HERE / "fig8_epitope_sets.json"
UNIPROT_OFFSET = 17                       # UniProt = construct + 17  (construct 2 = UniProt 19)

# --- chain A (TREM2) CA coordinates by construct res_id, from design_017 ---
ca = {}
for line in open(PDB):
    if line.startswith("ATOM") and line[12:16].strip() == "CA" and line[21] == "A":
        ca[int(line[22:26])] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
present = set(ca)
print(f"design_017 chain A: {len(ca)} CA, construct res {min(present)}..{max(present)}")

def to_construct(uniprot_list):
    return sorted(r for r in (u - UNIPROT_OFFSET for u in uniprot_list) if r in present)

refs = {row["reference"]: ast.literal_eval(row["uniprot_residues"]) for row in csv.DictReader(open(REFSETS))}
scfv = sorted(set(to_construct(refs["scFv2_6YYE"])) | set(to_construct(refs["scFv4_6Y6C"])))
ps = to_construct(refs["PS_site_6B8O"])
design = sorted(int(r["construct_pos"]) for r in csv.DictReader(open(FREQ))
                if float(r["freq_binder"]) >= 0.50 and int(r["construct_pos"]) in present)

sets = {"design_epitope": design, "scfv_epitope": scfv, "ps_site": ps}
for k, v in sets.items():
    print(f"{k:16s} n={len(v):3d}  construct {v}")

# --- geometry: is scFv really on the opposite face from design/PS? ---
center = np.mean(list(ca.values()), axis=0)
def cen(ids): return np.mean([ca[i] for i in ids], axis=0)
def unit(v): return v / np.linalg.norm(v)
vecs = {k: unit(cen(v) - center) for k, v in sets.items()}
def ang(a, b): return float(np.degrees(np.arccos(np.clip(np.dot(vecs[a], vecs[b]), -1, 1))))
print("\nOpposite-face check (angle between epitope-centroid directions from domain centre):")
print(f"  scFv  vs  design : {ang('scfv_epitope','design_epitope'):.0f} deg")
print(f"  scFv  vs  PS     : {ang('scfv_epitope','ps_site'):.0f} deg")
print(f"  design vs PS     : {ang('design_epitope','ps_site'):.0f} deg")

# overlap counts (Jaccard) to sanity-check the caption
def jac(a, b):
    A, B = set(sets[a]), set(sets[b])
    return len(A & B) / len(A | B) if (A | B) else 0.0
print(f"\nJaccard(design, PS)={jac('design_epitope','ps_site'):.2f}  "
      f"Jaccard(design, scFv)={jac('design_epitope','scfv_epitope'):.2f}")

json.dump({"sets": sets,
           "centroids_dir": {k: v.tolist() for k, v in vecs.items()},
           "center": center.tolist()}, open(OUT, "w"), indent=1)
print(f"\nwrote {OUT.name}")
