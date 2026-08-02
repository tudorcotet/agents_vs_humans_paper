"""Fig 8 (hero-preset version) — four same-orientation TREM2 panels, one per binding partner.

Every panel uses design_017's chain A as TREM2 (the shared frame), so TREM2 is IDENTICALLY oriented and
never rotated; each partner was superposed into that frame (superpose_crystals.py). Hero preset: TREM2 =
transparent purple surface over an opaque purple secondary-structure cartoon; the partner's interaction
epitope on TREM2 is highlighted orange; the binding partner is drawn on top.
  hero   design_017.pdb   our binder (chain B)         green cartoon
  scfv2  fig8_6yye.pdb    scFv-2 (Szykowska 2021)      slate cartoon
  scfv4  fig8_6y6c.pdb    scFv-4 (Szykowska 2021)      slate cartoon
  ps     fig8_6b8o.pdb    phosphatidylserine (Sudom)   yellow ball-and-stick
Because the surface is transparent, partners on the far (distal) face read through it.
Run:  /Applications/Blender.app/Contents/MacOS/Blender --background --python render_fig8_panels.py
"""
import csv, json, ast
import bpy
import numpy as np
import mathutils
from pathlib import Path
import molecularnodes as mn

HERE = Path(__file__).resolve().parent
SEC = HERE.parents[1]
PDBS = HERE / "pdbs"
OUT = HERE / "sample_renders"; OUT.mkdir(exist_ok=True)
REFSETS = SEC / "results" / "reference_epitope_sets.csv"
SETS = json.load(open(HERE / "fig8_epitope_sets.json"))
RES, SAMPLES = 1500, 140
UOFF = 17

TARGET_ALPHA = 0.22
PURPLE_S = [0.60, 0.48, 0.78]     # transparent TREM2 surface
PURPLE_C = [0.36, 0.24, 0.56]     # opaque TREM2 cartoon
ACCENT   = [0.91, 0.38, 0.17]     # design / PS epitope (orange)
BLUE_EPI = [0.13, 0.44, 0.90]     # antibody epitope (blue, contrasts the pink scFv)
GREEN    = [0.122, 0.894, 0.561]  # our binder
PINK     = [0.91, 0.42, 0.60]     # scFv antibody
YELLOW   = [0.97, 0.80, 0.15]     # PS lipid

# --- epitope residue sets on the construct frame (UniProt-17), per partner ---
def to_construct(u_list):
    return [u - UOFF for u in u_list if 1 <= (u - UOFF) <= 115]
refs = {r["reference"]: ast.literal_eval(r["uniprot_residues"]) for r in csv.DictReader(open(REFSETS))}
EPI = {
    "hero":  SETS["sets"]["design_epitope"],
    "scfv2": to_construct(refs["scFv2_6YYE"]),
    "scfv4": to_construct(refs["scFv4_6Y6C"]),
    "ps":    to_construct(refs["PS_site_6B8O"]),
}
PANELS = [
    ("hero",  "design_017.pdb", "B", "cartoon", GREEN,  ACCENT),
    ("scfv2", "fig8_6yye.pdb",  "P", "cartoon", PINK,   BLUE_EPI),
    ("scfv4", "fig8_6y6c.pdb",  "P", "cartoon", PINK,   BLUE_EPI),
    ("ps",    "fig8_6b8o.pdb",  "L", "ligand",  YELLOW, ACCENT),
]


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
                    chosen = backend; break
            except Exception:
                continue
    except Exception:
        pass
    scene.cycles.device = "GPU" if chosen else "CPU"
    scene.cycles.samples = samples
    scene.cycles.use_denoising = True
    for den in ("OPENIMAGEDENOISE", "OPTIX"):
        try:
            scene.cycles.denoiser = den; break
        except Exception:
            continue


def outline(scene):
    scene.render.use_freestyle = True
    vl = scene.view_layers[0]; vl.use_freestyle = True
    ls = vl.freestyle_settings.linesets[0] if len(vl.freestyle_settings.linesets) \
        else vl.freestyle_settings.linesets.new("outline")
    ls.select_silhouette = False; ls.select_border = False; ls.select_crease = False
    ls.select_external_contour = True; ls.select_material_boundary = True
    ls.linestyle.color = (0.10, 0.11, 0.13); ls.linestyle.thickness = 1.2


def render_panel(name, pdb, partner_chain, partner_style, partner_col, epi_col):
    for o in list(bpy.data.objects):
        bpy.data.objects.remove(o, do_unlink=True)
    mol = mn.Molecule.load(str(PDBS / pdb))
    obj = mol.object
    chain = np.array(mol.named_attribute("chain_id"))
    resid = np.array(mol.named_attribute("res_id"))
    Color = np.array(mol.named_attribute("Color"))
    is_A = chain == 0
    epi = is_A & np.isin(resid, EPI[name])
    Color[is_A] = PURPLE_S + [TARGET_ALPHA]
    Color[epi] = epi_col + [min(TARGET_ALPHA + 0.42, 1.0)]
    Color[~is_A] = partner_col + [1.0]
    mol.store_named_attribute(Color, "Color")

    mat_cart = bpy.data.materials.get("MAT_trem2_cartoon") or bpy.data.materials.new("MAT_trem2_cartoon")
    mat_cart.use_nodes = True
    b = next((n for n in mat_cart.node_tree.nodes if n.type == "BSDF_PRINCIPLED"), None)
    if b:
        b.inputs["Base Color"].default_value = PURPLE_C + [1.0]; b.inputs["Roughness"].default_value = 0.45

    selA = mn.MoleculeSelector(mol); selA.chain_id("A"); selA.store_selection("sel_t")
    selP = mn.MoleculeSelector(mol); selP.chain_id(partner_chain); selP.store_selection("sel_p")
    mol.add_style(style=mn.StyleSurface(quality=5), color=None, selection="sel_t", name="t_surf")
    mol.add_style(style=mn.StyleCartoon(peptide_dssp=False), selection="sel_t", material=mat_cart, name="t_cart")
    if partner_style == "cartoon":
        mol.add_style(style=mn.StyleCartoon(peptide_dssp=False), color=None, selection="sel_p", name="p_cart")
    else:
        mol.add_style(style=mn.StyleSpheres(), color=None, selection="sel_p", name="p_lig")   # PS lipid, chunky
    for mat in bpy.data.materials:
        if hasattr(mat, "surface_render_method"):
            mat.surface_render_method = "BLENDED"

    # --- fixed camera + framing from chain A world coords (identical across panels) ---
    pos = np.array(mol.named_attribute("position"))
    M = obj.matrix_world
    world = pos @ np.array(M.to_3x3()).T + np.array(M.translation)
    wA = world[is_A]
    center = mathutils.Vector(wA.mean(axis=0))
    radius = float(np.linalg.norm(wA - wA.mean(axis=0), axis=1).max()) or 0.4

    target = bpy.data.objects.new("T", None); bpy.context.collection.objects.link(target); target.location = center
    cam = bpy.data.objects.new("C", bpy.data.cameras.new("C")); cam.data.lens = 52
    bpy.context.collection.objects.link(cam)
    cam.location = center + mathutils.Vector((0.35, -1.0, 0.28)).normalized() * radius * 3.7
    cam.constraints.new("TRACK_TO").target = target
    bpy.context.scene.camera = cam
    for d, e, s in [((1.2, -1, 1.2), 28, 1.0), ((-1.4, -.6, .4), 9, 1.6), ((0, 1.4, 1), 16, 1.0)]:
        ld = bpy.data.lights.new("L", "AREA"); ld.energy = e; ld.size = radius * s
        lo = bpy.data.objects.new("L", ld); bpy.context.collection.objects.link(lo)
        lo.location = center + mathutils.Vector(d).normalized() * radius * 2.5
        lo.constraints.new("TRACK_TO").target = target
    w = bpy.data.worlds.new("W"); bpy.context.scene.world = w; w.use_nodes = True
    bg = w.node_tree.nodes.get("Background")
    bg.inputs["Color"].default_value = (0.95, 0.96, 0.97, 1.0); bg.inputs["Strength"].default_value = 0.4

    sc = bpy.context.scene
    setup_cycles(sc, SAMPLES); outline(sc)
    sc.render.film_transparent = True
    sc.render.resolution_x = sc.render.resolution_y = RES
    sc.view_settings.exposure = -0.6
    sc.render.filepath = str(OUT / f"fig8p_{name}")
    bpy.ops.render.render(write_still=True)
    print("SAVED", name)


for spec in PANELS:
    render_panel(*spec)
print("DONE fig8 panels")
