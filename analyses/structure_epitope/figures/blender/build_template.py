"""Build a STARTER template.blend for Fig 10 / Fig 8 — a Blender + Molecular Nodes scene the operator opens
and refines. Generated headlessly so it is reproducible and reviewable.

Look & quality (publication defaults):
  - engine = CYCLES on GPU (Metal on macOS / OptiX-CUDA on NVIDIA; CPU fallback) with denoising
    (OpenImageDenoise; OptiX denoiser on NVIDIA),
  - TREM2 = a transparent PURPLE surface (apical epitope, freq >= 0.5, highlighted orange) over an OPAQUE
    purple Cartoon of TREM2's own secondary structure, so the fold shows through the surface; the surface
    reads the per-atom Color attribute (alpha = opacity) while the cartoon uses a flat opaque material;
    chain B = human-green Cartoon,
  - a thin dark Freestyle silhouette outline + a shadow-catcher plane (contact shadow) so the light surface
    separates from a white background, exposure trimmed -0.6,
  - camera framing + 3-point area lights, transparent film, high-res RGBA PNG.

Toggles below: COLOUR_BY_EPITOPE, TARGET_ALPHA (look); RES, SAMPLES (quality).

Run:
  /Applications/Blender.app/Contents/MacOS/Blender --background --python build_template.py
Verified on Blender 5.0.1 + Molecular Nodes (bl_ext.blender_org.molecularnodes).
"""
import csv
import bpy
import numpy as np
import mathutils
from pathlib import Path
import molecularnodes as mn

HERE = Path(__file__).resolve().parent
SEC = HERE.parents[1]
PDB = HERE / "pdbs" / "design_017.pdb"
FREQ = SEC / "results" / "epitope_per_residue_freq.csv"
OUT = HERE / "template.blend"

# ---- look ----
COLOUR_BY_EPITOPE = True
TARGET_ALPHA = 0.22             # low, so the TREM2 cartoon underneath shows clearly through the surface
EPITOPE_FREQ_CUTOFF = 0.50
# ---- quality ----
RES = 2000
SAMPLES = 160

ACCENT   = [0.91, 0.38, 0.17]     # epitope orange
PURPLE_S = [0.60, 0.48, 0.78]     # #9A7BC7 TREM2 surface (transparent)
PURPLE_C = [0.36, 0.24, 0.56]     # #5C3D8F TREM2 secondary-structure cartoon (opaque, shows through)
GREEN    = [0.122, 0.894, 0.561]  # #1FE48F human binder


def setup_cycles(scene, samples):
    scene.render.engine = "CYCLES"
    chosen = None
    try:
        prefs = bpy.context.preferences.addons["cycles"].preferences
        for backend in ("METAL", "OPTIX", "CUDA", "HIP", "ONEAPI"):
            try:
                prefs.compute_device_type = backend
                prefs.get_devices()
                if any(d.type != "CPU" for d in prefs.devices):
                    for d in prefs.devices:
                        d.use = True
                    chosen = backend
                    break
            except Exception:
                continue
    except Exception:
        pass
    scene.cycles.device = "GPU" if chosen else "CPU"
    scene.cycles.samples = samples
    scene.cycles.use_denoising = True
    for den in ("OPENIMAGEDENOISE", "OPTIX"):   # OIDN works everywhere; OptiX on NVIDIA
        try:
            scene.cycles.denoiser = den
            break
        except Exception:
            continue
    print(f"CYCLES device={scene.cycles.device} backend={chosen or 'CPU'} "
          f"samples={samples} denoiser={scene.cycles.denoiser}")


for o in list(bpy.data.objects):
    bpy.data.objects.remove(o, do_unlink=True)

mol = mn.Molecule.load(str(PDB))
obj = mol.object

# --- per-atom Color: region colour on chain A, green on B; alpha = opacity ---
chain = np.array(mol.named_attribute("chain_id"))
resid = np.array(mol.named_attribute("res_id"))
Color = np.array(mol.named_attribute("Color"))
is_A = chain == 0
epi = {int(r["construct_pos"]) for r in csv.DictReader(open(FREQ))
       if float(r["freq_binder"]) >= EPITOPE_FREQ_CUTOFF}
epi_mask = is_A & np.isin(resid, list(epi))
Color[is_A] = PURPLE_S + [TARGET_ALPHA]
if COLOUR_BY_EPITOPE:
    Color[epi_mask] = ACCENT + [min(TARGET_ALPHA + 0.40, 1.0)]
Color[~is_A] = GREEN + [1.0]
mol.store_named_attribute(Color, "Color")

# OPAQUE material for the TREM2 cartoon — a flat material overrides the Color attribute, so the cartoon
# keeps alpha=1 and stays visible under the transparent surface (which reads the Color attribute's alpha).
mat_cart = bpy.data.materials.new("MAT_trem2_cartoon"); mat_cart.use_nodes = True
_bsdf = next((n for n in mat_cart.node_tree.nodes if n.type == "BSDF_PRINCIPLED"), None)
if _bsdf:
    _bsdf.inputs["Base Color"].default_value = PURPLE_C + [1.0]
    _bsdf.inputs["Roughness"].default_value = 0.45

selA = mn.MoleculeSelector(mol); selA.chain_id("A"); selA.store_selection("sel_target")
selB = mn.MoleculeSelector(mol); selB.chain_id("B"); selB.store_selection("sel_binder")
# TREM2: transparent purple Surface (reads Color) + opaque purple Cartoon underneath (flat material).
mol.add_style(style=mn.StyleSurface(quality=5), color=None, selection="sel_target", name="target_surface")
# peptide_dssp=False (NOT True): MN's on-the-fly DSSP recompute silently truncates this ribbon to ~1/3 of
# the domain (verified: eval mesh 1044 verts / span 0.15,0.37,0.19 vs CA span 0.30,0.42,0.28). The structure
# already carries a complete sheet/loop annotation from load, so =False draws the FULL beta-sandwich w/ arrows
# (eval mesh 5530 verts / span 0.31,0.43,0.29, matching the CA extent).
mol.add_style(style=mn.StyleCartoon(peptide_dssp=False), selection="sel_target", material=mat_cart, name="target_cartoon")
mol.add_style(style=mn.StyleCartoon(), color=None, selection="sel_binder", name="binder_cartoon")

# Cycles honours BSDF alpha natively; also set EEVEE's BLENDED in case the operator switches engines.
for mat in bpy.data.materials:
    if hasattr(mat, "surface_render_method"):
        mat.surface_render_method = "BLENDED"

# --- centroid + radius ---
bb = [obj.matrix_world @ mathutils.Vector(c) for c in obj.bound_box]
center = sum(bb, mathutils.Vector()) / 8.0
radius = max((v - center).length for v in bb) or 10.0
floor_z = min(v.z for v in bb)

target = bpy.data.objects.new("CAM_TARGET", None)
bpy.context.collection.objects.link(target); target.location = center

cam_data = bpy.data.cameras.new("Camera"); cam_data.lens = 55
cam = bpy.data.objects.new("Camera", cam_data)
bpy.context.collection.objects.link(cam)
cam.location = center + mathutils.Vector((0.35, -1.0, 0.28)).normalized() * radius * 3.0
cam.constraints.new("TRACK_TO").target = target
bpy.context.scene.camera = cam


def add_area(name, direction, energy, size):
    ld = bpy.data.lights.new(name, "AREA"); ld.energy, ld.size = energy, size
    lo = bpy.data.objects.new(name, ld); bpy.context.collection.objects.link(lo)
    lo.location = center + mathutils.Vector(direction).normalized() * radius * 2.5
    lo.constraints.new("TRACK_TO").target = target

add_area("KEY",  (1.2, -1.0, 1.2), 26.0, radius)
add_area("FILL", (-1.4, -0.6, 0.4), 9.0, radius * 1.6)
add_area("RIM",  (0.0, 1.4, 1.0), 16.0, radius)

world = bpy.data.worlds.new("World"); bpy.context.scene.world = world
world.use_nodes = True
bg = world.node_tree.nodes.get("Background")
if bg:
    bg.inputs["Color"].default_value = (0.95, 0.96, 0.97, 1.0)
    bg.inputs["Strength"].default_value = 0.4

sc = bpy.context.scene
setup_cycles(sc, SAMPLES)
sc.render.film_transparent = True
sc.render.resolution_x = sc.render.resolution_y = RES
sc.render.image_settings.file_format = "PNG"
sc.render.image_settings.color_mode = "RGBA"
sc.view_settings.exposure = -0.6

# thin dark OUTER outline so the light surface separates from a white page.
# external_contour = outer silhouette only (no internal folds -> no scratchy lines);
# material_boundary traces where the binder meets the surface.
sc.render.use_freestyle = True
vl = sc.view_layers[0]; vl.use_freestyle = True
fs = vl.freestyle_settings
ls = fs.linesets[0] if len(fs.linesets) else fs.linesets.new("outline")
ls.select_silhouette = False; ls.select_border = False
ls.select_crease = False; ls.select_edge_mark = False
ls.select_external_contour = True
ls.select_material_boundary = True
ls.linestyle.color = (0.10, 0.11, 0.13); ls.linestyle.thickness = 1.3

OUT.parent.mkdir(parents=True, exist_ok=True)
bpy.ops.wm.save_as_mainfile(filepath=str(OUT))
print(f"SAVED {OUT.name} | epitope={COLOUR_BY_EPITOPE} alpha={TARGET_ALPHA} res={RES} "
      f"| {len(epi)} epitope residues | r={radius:.1f}")
