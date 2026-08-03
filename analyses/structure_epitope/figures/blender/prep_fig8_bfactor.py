"""Prep for Fig 8 — bake TREM2 per-residue binder-contact frequency (and, optionally, the EvoEF2 hotspot
score) into a single TREM2 target chain so Molecular Nodes can colour the surface by ``b_factor`` with a
colour ramp. No per-residue attribute wrangling needed inside Blender.

Inputs
  figures/blender/pdbs/design_017.pdb        chain A = TREM2 IgSF target (pre-superposed common frame)
  results/epitope_per_residue_freq.csv       freq_binder per construct_pos (== PDB residue number: Met=1)

Output
  figures/blender/trem2_target_epitope.pdb   chain A only; B-factor column = 100 * freq_binder  (0-100)

Run
  uv run --with gemmi --with pandas python figures/blender/prep_fig8_bfactor.py

Notes
  - PDB residue number == construct position (Met=1, UniProt 19 -> residue 2); UniProt = construct + 17.
    Verify the offset once (the printed sanity line shows the tip residues L69/L71/F74 lighting up).
  - For panel (c) hotspot labels, the top EvoEF2 target hotspots are in results/alanine_scan.csv.
"""
import gemmi
import pandas as pd
from pathlib import Path

HERE = Path(__file__).resolve().parent
SEC = HERE.parents[1]                       # analyses/structure_epitope
freq = pd.read_csv(SEC / "results/epitope_per_residue_freq.csv").set_index("construct_pos")["freq_binder"]

st = gemmi.read_structure(str(HERE / "pdbs/design_017.pdb"))
model = st[0]
for ch in [c.name for c in model if c.name != "A"]:   # keep TREM2 target (chain A) only
    model.remove_chain(ch)

hot = []
for res in model["A"]:
    f = float(freq.get(res.seqid.num, 0.0))
    for atom in res:
        atom.b_iso = round(100.0 * f, 2)
    if f >= 0.9:
        hot.append(f"{res.name}{res.seqid.num}")

st.setup_entities()
out = HERE / "trem2_target_epitope.pdb"
st.write_pdb(str(out))
print(f"wrote {out.name}: chain A only, B-factor = 100*freq_binder (0-100).")
print(f"sanity — residues with binder contact freq >=0.90 (should be the apical tip): {', '.join(hot)}")
