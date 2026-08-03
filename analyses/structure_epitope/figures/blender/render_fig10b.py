"""Fig 10b — close-up of the best design (design_017, 1.11 nM) interface, with the top EvoEF2 hotspot
TREM2 side chains as red sticks. Hero preset (transparent purple TREM2 surface + fold, green binder).
Also writes fig10b_labels.json = projected 2-D pixel positions of each hotspot, so compose_fig10b.py can
draw labels.  Run:  /Applications/Blender.app/Contents/MacOS/Blender --background --python render_fig10b.py
"""
import json
import bpy
import numpy as np
import mathutils
from pathlib import Path
from bpy_extras.object_utils import world_to_camera_view
import molecularnodes as mn

HERE = Path(__file__).resolve().parent
PDB = HERE / "pdbs" / "design_017.pdb"
OUT = HERE / "sample_renders"; OUT.mkdir(exist_ok=True)
RES, SAMPLES = 1600, 160

TARGET_ALPHA = 0.20
PURPLE_S = [0.60, 0.48, 0.78]
PURPLE_C = [0.36, 0.24, 0.56]
GREEN    = [0.122, 0.894, 0.561]
HOT      = [0.92, 0.22, 0.16]     # hotspot sticks (red)
# top design-17 TREM2 hotspots by EvoEF2 ddG (construct res -> UniProt label)
HOTSPOTS = {27: "W44", 55: "L72", 54: "L71", 52: "L69", 58: "L75", 57: "F74"}


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
resid = np.array(mol.named_attribute("res_id"))
Color = np.array(mol.named_attribute("Color"))
is_A = chain == 0
hot_ids = list(HOTSPOTS)
is_hot = is_A & np.isin(resid, hot_ids)
Color[is_A] = PURPLE_S + [TARGET_ALPHA]
Color[~is_A] = GREEN + [1.0]
mol.store_named_attribute(Color, "Color")
mol.store_named_attribute(is_hot, "sel_hot")

def flat_mat(name, rgb, rough=0.4):
    m = bpy.data.materials.new(name); m.use_nodes = True
    b = next((n for n in m.node_tree.nodes if n.type == "BSDF_PRINCIPLED"), None)
    if b: b.inputs["Base Color"].default_value = rgb + [1.0]; b.inputs["Roughness"].default_value = rough
    return m
mat_cart = flat_mat("MAT_trem2_cartoon", PURPLE_C, 0.45)
mat_hot = flat_mat("MAT_hot", HOT, 0.3)

selA = mn.MoleculeSelector(mol); selA.chain_id("A"); selA.store_selection("sel_t")
selB = mn.MoleculeSelector(mol); selB.chain_id("B"); selB.store_selection("sel_b")
mol.add_style(style=mn.StyleSurface(quality=6), color=None, selection="sel_t", name="t_surf")
mol.add_style(style=mn.StyleCartoon(peptide_dssp=False), selection="sel_t", material=mat_cart, name="t_cart")
mol.add_style(style=mn.StyleCartoon(peptide_dssp=False), color=None, selection="sel_b", name="binder")
mol.add_style(style=mn.StyleSticks(radius=0.42, quality=3), selection="sel_hot", material=mat_hot, name="hotspots")
for mat in bpy.data.materials:
    if hasattr(mat, "surface_render_method"):
        mat.surface_render_method = "BLENDED"

# --- world coords; frame the interface (target = hotspot centroid) and project hotspot labels ---
pos = np.array(mol.named_attribute("position"))
M = obj.matrix_world
world = pos @ np.array(M.to_3x3()).T + np.array(M.translation)
wA = world[is_A]
domain_center = wA.mean(axis=0)
radius = float(np.linalg.norm(wA - wA.mean(axis=0), axis=1).max()) or 0.4
hot_center = world[is_hot].mean(axis=0)
# aim at the interface (weighted toward the hotspots) and zoom in — this subpanel only
target_pt = mathutils.Vector(0.30 * domain_center + 0.70 * hot_center)

target = bpy.data.objects.new("T", None); bpy.context.collection.objects.link(target); target.location = target_pt
cam = bpy.data.objects.new("C", bpy.data.cameras.new("C")); cam.data.lens = 52
bpy.context.collection.objects.link(cam)
cam.location = target_pt + mathutils.Vector((0.35, -1.0, 0.28)).normalized() * radius * 2.15
cam.constraints.new("TRACK_TO").target = target
bpy.context.scene.camera = cam
for d, e, s in [((1.2, -1, 1.2), 30, 1.0), ((-1.4, -.6, .4), 10, 1.6), ((0, 1.4, 1), 16, 1.0)]:
    ld = bpy.data.lights.new("L", "AREA"); ld.energy = e; ld.size = radius * s
    lo = bpy.data.objects.new("L", ld); bpy.context.collection.objects.link(lo)
    lo.location = mathutils.Vector(domain_center) + mathutils.Vector(d).normalized() * radius * 2.6
    lo.constraints.new("TRACK_TO").target = target
w = bpy.data.worlds.new("W"); bpy.context.scene.world = w; w.use_nodes = True
bgn = w.node_tree.nodes.get("Background")
bgn.inputs["Color"].default_value = (0.95, 0.96, 0.97, 1.0); bgn.inputs["Strength"].default_value = 0.4

sc = bpy.context.scene
setup_cycles(sc, SAMPLES); outline(sc)
sc.render.film_transparent = True
sc.render.resolution_x = sc.render.resolution_y = RES
sc.view_settings.exposure = -0.6
bpy.context.view_layer.update()

# project each hotspot to pixel coords (for compose) AND store its LOCAL centroid (for re-projection after
# the operator rotates the camera in the saved .blend — see rerender_fig10b.py)
labels, hot_local = {}, {}
for cid, name in HOTSPOTS.items():
    m = is_A & (resid == cid)
    if not m.any(): continue
    co = mathutils.Vector(world[m].mean(axis=0))
    ndc = world_to_camera_view(sc, cam, co)
    labels[name] = [round(ndc.x * RES, 1), round((1 - ndc.y) * RES, 1), round(ndc.z, 3)]
    hot_local[name] = [round(float(v), 4) for v in pos[m].mean(axis=0)]
json.dump({"res": RES, "labels": labels}, open(HERE / "fig10b_labels.json", "w"), indent=1)
json.dump({"object": obj.name, "hotspots": hot_local}, open(HERE / "fig10b_hotspots.json", "w"), indent=1)

sc.render.filepath = str(OUT / "fig10b_design017")
bpy.ops.render.render(write_still=True)

# Save an operator-editable scene: bake the camera's current orientation and DROP the TRACK_TO constraint so
# the operator can freely rotate the camera (a track-to would just re-aim it at the target and undo the edit).
cam_mw = cam.matrix_world.copy()
for c in list(cam.constraints):
    cam.constraints.remove(c)
cam.matrix_world = cam_mw
bpy.ops.wm.save_as_mainfile(filepath=str(HERE / "fig10b.blend"))
print("SAVED fig10b_design017 + fig10b.blend | labels", labels)
