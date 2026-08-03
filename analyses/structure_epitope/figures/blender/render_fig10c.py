"""Fig 10c — the 37 hit binders superposed on one TREM2, coloured by cohort (human green / agent blue).
Hero preset: TREM2 = transparent purple surface + opaque purple fold; binders = cohort-coloured cartoons.
Cohort comes from the B-factor baked by prep_fig10c.py (TREM2=0, human=1, agent=2). Same TREM2 camera as
the gallery.  Run:  /Applications/Blender.app/Contents/MacOS/Blender --background --python render_fig10c.py
"""
import bpy
import numpy as np
import mathutils
from pathlib import Path
import molecularnodes as mn

HERE = Path(__file__).resolve().parent
PDB = HERE / "pdbs" / "fig10c_superposed.pdb"
OUT = HERE / "sample_renders"; OUT.mkdir(exist_ok=True)
RES, SAMPLES = 1600, 150

TARGET_ALPHA = 0.20
PURPLE_S = [0.60, 0.48, 0.78]
PURPLE_C = [0.36, 0.24, 0.56]
HUMAN    = [0.122, 0.894, 0.561]
AGENT    = [0.188, 0.773, 0.961]


def setup_cycles(scene, samples):
    scene.render.engine = "CYCLES"; chosen = None
    try:
        prefs = bpy.context.preferences.addons["cycles"].preferences
        for backend in ("METAL", "OPTIX", "CUDA", "HIP", "ONEAPI"):
            try:
                prefs.compute_device_type = backend; prefs.get_devices()
                if any(d.type != "CPU" for d in prefs.devices):
                    for d in prefs.devices: d.use = True
                    chosen = backend; break
            except Exception:
                continue
    except Exception:
        pass
    scene.cycles.device = "GPU" if chosen else "CPU"
    scene.cycles.samples = samples; scene.cycles.use_denoising = True
    for den in ("OPENIMAGEDENOISE", "OPTIX"):
        try: scene.cycles.denoiser = den; break
        except Exception: continue


def outline(scene):
    scene.render.use_freestyle = True
    vl = scene.view_layers[0]; vl.use_freestyle = True
    ls = vl.freestyle_settings.linesets[0] if len(vl.freestyle_settings.linesets) \
        else vl.freestyle_settings.linesets.new("outline")
    ls.select_silhouette = False; ls.select_border = False; ls.select_crease = False
    ls.select_external_contour = True; ls.select_material_boundary = True
    ls.linestyle.color = (0.10, 0.11, 0.13); ls.linestyle.thickness = 1.1


for o in list(bpy.data.objects):
    bpy.data.objects.remove(o, do_unlink=True)

mol = mn.Molecule.load(str(PDB))
obj = mol.object
chain = np.array(mol.named_attribute("chain_id"))
bf = np.array(mol.named_attribute("b_factor"))
Color = np.array(mol.named_attribute("Color"))
is_t = chain == 0
is_h = (~is_t) & (bf < 1.5)
is_a = (~is_t) & (bf >= 1.5)
Color[is_t] = PURPLE_S + [TARGET_ALPHA]
Color[is_h] = HUMAN + [1.0]
Color[is_a] = AGENT + [1.0]
mol.store_named_attribute(Color, "Color")
mol.store_named_attribute(is_t, "sel_t")
mol.store_named_attribute(~is_t, "sel_b")

mat_cart = bpy.data.materials.new("MAT_trem2_cartoon"); mat_cart.use_nodes = True
b = next((n for n in mat_cart.node_tree.nodes if n.type == "BSDF_PRINCIPLED"), None)
if b:
    b.inputs["Base Color"].default_value = PURPLE_C + [1.0]; b.inputs["Roughness"].default_value = 0.45

mol.add_style(style=mn.StyleSurface(quality=4), color=None, selection="sel_t", name="t_surf")
mol.add_style(style=mn.StyleCartoon(peptide_dssp=False), selection="sel_t", material=mat_cart, name="t_cart")
mol.add_style(style=mn.StyleCartoon(peptide_dssp=False), color=None, selection="sel_b", name="binders")
for mat in bpy.data.materials:
    if hasattr(mat, "surface_render_method"):
        mat.surface_render_method = "BLENDED"

pos = np.array(mol.named_attribute("position"))
M = obj.matrix_world
world = pos @ np.array(M.to_3x3()).T + np.array(M.translation)
wT = world[is_t]
center = mathutils.Vector(wT.mean(axis=0))
radius = float(np.linalg.norm(wT - wT.mean(axis=0), axis=1).max()) or 0.4

target = bpy.data.objects.new("T", None); bpy.context.collection.objects.link(target); target.location = center
cam = bpy.data.objects.new("C", bpy.data.cameras.new("C")); cam.data.lens = 50
bpy.context.collection.objects.link(cam)
cam.location = center + mathutils.Vector((0.35, -1.0, 0.28)).normalized() * radius * 4.0
cam.constraints.new("TRACK_TO").target = target
bpy.context.scene.camera = cam
for d, e, s in [((1.2, -1, 1.2), 30, 1.0), ((-1.4, -.6, .4), 10, 1.6), ((0, 1.4, 1), 16, 1.0)]:
    ld = bpy.data.lights.new("L", "AREA"); ld.energy = e; ld.size = radius * s
    lo = bpy.data.objects.new("L", ld); bpy.context.collection.objects.link(lo)
    lo.location = center + mathutils.Vector(d).normalized() * radius * 2.6
    lo.constraints.new("TRACK_TO").target = target
w = bpy.data.worlds.new("W"); bpy.context.scene.world = w; w.use_nodes = True
bgn = w.node_tree.nodes.get("Background")
bgn.inputs["Color"].default_value = (0.95, 0.96, 0.97, 1.0); bgn.inputs["Strength"].default_value = 0.4

sc = bpy.context.scene
setup_cycles(sc, SAMPLES); outline(sc)
sc.render.film_transparent = True
sc.render.resolution_x = sc.render.resolution_y = RES
sc.view_settings.exposure = -0.6
sc.render.filepath = str(OUT / "fig10c_superposed")
bpy.ops.render.render(write_still=True)
print("SAVED fig10c_superposed")
