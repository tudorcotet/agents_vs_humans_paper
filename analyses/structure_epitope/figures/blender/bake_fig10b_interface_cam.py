"""Bake a 'straight-into-the-interface' camera into fig10b.blend (looking from the TREM2 side along the
TREM2->binder axis, hotspots centred). Two modes:

  # 1) preview several candidate angles (fast, small):
  BL --background --python bake_fig10b_interface_cam.py
  # 2) bake the chosen one into fig10b.blend (+ full preview + hotspots json):
  BL --background --python bake_fig10b_interface_cam.py -- --bake B_down35

Scene is identical to render_fig10b.py (transparent purple TREM2 + fold, green binder, red hotspot sticks);
only the camera differs. The saved camera is a FREE object (no track-to) so the operator can fine-tune.
"""
import sys, json
import bpy
import numpy as np
import mathutils
from pathlib import Path
from bpy_extras.object_utils import world_to_camera_view
import molecularnodes as mn

HERE = Path(__file__).resolve().parent
PDB = HERE / "pdbs" / "design_017.pdb"
OUT = HERE / "sample_renders"; OUT.mkdir(exist_ok=True)
CAND_DIR = Path("/private/tmp/claude-501/-Users-amyhe-Desktop-trem2-2026/"
                "d57e2781-0f96-4e4d-848e-74a070e203f5/scratchpad/fig10b_cand")
CAND_DIR.mkdir(parents=True, exist_ok=True)
argv = sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else []
BAKE = argv[argv.index("--bake") + 1] if "--bake" in argv else None

TARGET_ALPHA = 0.20
PURPLE_S = [0.60, 0.48, 0.78]; PURPLE_C = [0.36, 0.24, 0.56]
GREEN = [0.122, 0.894, 0.561]; HOT = [0.92, 0.22, 0.16]
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


def flat_mat(name, rgb, rough=0.4):
    m = bpy.data.materials.new(name); m.use_nodes = True
    b = next((n for n in m.node_tree.nodes if n.type == "BSDF_PRINCIPLED"), None)
    if b: b.inputs["Base Color"].default_value = rgb + [1.0]; b.inputs["Roughness"].default_value = rough
    return m


def unit(v): return np.array(v) / np.linalg.norm(v)

# ---- build scene ----
for o in list(bpy.data.objects):
    bpy.data.objects.remove(o, do_unlink=True)
mol = mn.Molecule.load(str(PDB)); obj = mol.object
chain = np.array(mol.named_attribute("chain_id"))
resid = np.array(mol.named_attribute("res_id"))
Color = np.array(mol.named_attribute("Color"))
is_A = chain == 0
is_hot = is_A & np.isin(resid, list(HOTSPOTS))
Color[is_A] = PURPLE_S + [TARGET_ALPHA]; Color[~is_A] = GREEN + [1.0]
mol.store_named_attribute(Color, "Color")
mol.store_named_attribute(is_hot, "sel_hot")
mat_cart = flat_mat("MAT_trem2_cartoon", PURPLE_C, 0.45); mat_hot = flat_mat("MAT_hot", HOT, 0.3)
selA = mn.MoleculeSelector(mol); selA.chain_id("A"); selA.store_selection("sel_t")
selB = mn.MoleculeSelector(mol); selB.chain_id("B"); selB.store_selection("sel_b")
mol.add_style(style=mn.StyleSurface(quality=6), color=None, selection="sel_t", name="t_surf")
mol.add_style(style=mn.StyleCartoon(peptide_dssp=False), selection="sel_t", material=mat_cart, name="t_cart")
mol.add_style(style=mn.StyleCartoon(peptide_dssp=False), color=None, selection="sel_b", name="binder")
mol.add_style(style=mn.StyleSticks(radius=0.42, quality=3), selection="sel_hot", material=mat_hot, name="hotspots")
for mat in bpy.data.materials:
    if hasattr(mat, "surface_render_method"):
        mat.surface_render_method = "BLENDED"

# ---- geometry (Blender world space) ----
pos = np.array(mol.named_attribute("position"))
M = obj.matrix_world
world = pos @ np.array(M.to_3x3()).T + np.array(M.translation)
domain_c = world[is_A].mean(0)
binder_c = world[~is_A].mean(0)
hot_c = world[is_hot].mean(0)
radius = float(np.linalg.norm(world[is_A] - domain_c, axis=1).max()) or 0.4
axis = unit(binder_c - domain_c)                       # TREM2 -> binder
up0 = np.array([0, 0, 1.0]) if abs(axis @ [0, 0, 1]) < 0.9 else np.array([0, 1.0, 0])
side = unit(np.cross(axis, up0)); vert = unit(np.cross(side, axis))
if vert[2] < 0: vert = -vert                            # vert points up-ish
focus = mathutils.Vector(0.55 * hot_c + 0.45 * binder_c)   # interface, weighted to the hotspots

# SIDE-ON views: camera roughly PERPENDICULAR to the TREM2->binder axis (only a small axial tilt K), so it
# looks broadside at the interface band. Sweep the azimuth around the axis to find the angle that best reveals
# the hotspot sticks between binder and TREM2.
K = 0.22                                                # small axial tilt (~12 deg off pure perpendicular)
def perp(a_deg):
    t = np.radians(a_deg)
    return np.cos(t) * side + np.sin(t) * vert
CANDS = {f"S{a:03d}": unit(K * axis + perp(a)) for a in (0, 45, 90, 135, 180, 225, 270, 315)}
DFAC = 2.2

# common lights/world/render
target = bpy.data.objects.new("CAM_TARGET", None); bpy.context.collection.objects.link(target)
target.location = focus
cam = bpy.data.objects.new("Camera", bpy.data.cameras.new("Camera")); cam.data.lens = 46
bpy.context.collection.objects.link(cam); bpy.context.scene.camera = cam
trk = cam.constraints.new("TRACK_TO"); trk.target = target
for d, e, s in [((1.2, -1, 1.2), 30, 1.0), ((-1.4, -.6, .4), 10, 1.6), ((0, 1.4, 1), 16, 1.0)]:
    ld = bpy.data.lights.new("L", "AREA"); ld.energy = e; ld.size = radius * s
    lo = bpy.data.objects.new("L", ld); bpy.context.collection.objects.link(lo)
    lo.location = mathutils.Vector(domain_c) + mathutils.Vector(d).normalized() * radius * 2.6
    lo.constraints.new("TRACK_TO").target = target
w = bpy.data.worlds.new("W"); bpy.context.scene.world = w; w.use_nodes = True
bgn = w.node_tree.nodes.get("Background")
bgn.inputs["Color"].default_value = (0.95, 0.96, 0.97, 1.0); bgn.inputs["Strength"].default_value = 0.4
sc = bpy.context.scene
sc.render.film_transparent = True
sc.view_settings.exposure = -0.6

def place(view_dir):
    cam.location = focus - mathutils.Vector(view_dir) * radius * DFAC

if BAKE:
    setup_cycles(sc, 160); outline(sc)
    sc.render.resolution_x = sc.render.resolution_y = 1600
    place(CANDS[BAKE])
    bpy.context.view_layer.update()
    # project hotspot labels + store local centroids
    labels, hot_local = {}, {}
    RES = sc.render.resolution_x
    for cid, name in HOTSPOTS.items():
        m = is_A & (resid == cid)
        if not m.any(): continue
        co = mathutils.Vector(world[m].mean(0))
        ndc = world_to_camera_view(sc, cam, co)
        labels[name] = [round(ndc.x * RES, 1), round((1 - ndc.y) * RES, 1), round(ndc.z, 3)]
        hot_local[name] = [round(float(v), 4) for v in pos[m].mean(0)]
    json.dump({"res": RES, "labels": labels}, open(HERE / "fig10b_labels.json", "w"), indent=1)
    json.dump({"object": obj.name, "hotspots": hot_local}, open(HERE / "fig10b_hotspots.json", "w"), indent=1)
    sc.render.filepath = str(OUT / "fig10b_design017")
    bpy.ops.render.render(write_still=True)
    # free the camera (bake matrix, drop track-to) and save
    mw = cam.matrix_world.copy()
    for c in list(cam.constraints): cam.constraints.remove(c)
    cam.matrix_world = mw
    bpy.ops.wm.save_as_mainfile(filepath=str(HERE / "fig10b.blend"))
    print(f"BAKED {BAKE} -> fig10b.blend + fig10b_design017.png | labels {labels}")
else:
    setup_cycles(sc, 48); outline(sc)
    sc.render.resolution_x = sc.render.resolution_y = 760
    for label, vd in CANDS.items():
        place(vd)
        sc.render.filepath = str(CAND_DIR / f"cand_{label}")
        bpy.ops.render.render(write_still=True)
        print("CAND", label)
    print("CANDS_DONE")
