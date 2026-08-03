"""Headless batch-render for the Fig 10a leaderboard gallery — the WORKING counterpart to the
render_headless.py skeleton (verified Blender 5.0.1 + Molecular Nodes).

For the top-N designs (by leaderboard_rank in manifest.csv) it renders each complex in the NEUTRAL gallery
look — a TRANSPARENT slate-grey TREM2 surface over an OPAQUE grey secondary-structure Cartoon (so the fold
and the buried binder read through the surface) + binder Cartoon in cohort colour (human #1FE48F, agent
#30C5F5) — on the shared, pre-superposed TREM2 frame, transparent film. Quality: CYCLES on GPU (Metal/OptiX)
with OpenImageDenoise, a thin outer Freestyle outline + exposure trim so the surface reads on a white page.
Tiles land in sample_renders/tile_<id>.png; compose with the PIL tiler in FIGURE_GUIDE.md. A committed
sample grid is sample_renders/fig10a_sample_grid.png.

Run:
  /Applications/Blender.app/Contents/MacOS/Blender --background --python render_gallery.py -- --n 16
Knobs: --n N (default 16); RES / SAMPLES below. Cycles + N tiles is minutes, not seconds — expected.
"""
import sys
import csv
import bpy
import numpy as np
import mathutils
from pathlib import Path
import molecularnodes as mn

HERE = Path(__file__).resolve().parent
argv = sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else []
N = int(argv[argv.index("--n") + 1]) if "--n" in argv else 16
# --ids 17,74,102,34 renders an explicit set (e.g. the Kd-ranked grid) instead of the top-N by rank
IDS = [int(x) for x in argv[argv.index("--ids") + 1].split(",")] if "--ids" in argv else None
RES, SAMPLES = 1200, 110

HUMAN = [0.122, 0.894, 0.561]
AGENT = [0.188, 0.773, 0.961]
GREY  = [0.45, 0.49, 0.54]     # #737D8A mid slate-grey (reads on white)
TARGET_ALPHA = 0.25            # transparent TREM2 surface; fold + buried binder show through

rows = list(csv.DictReader(open(HERE / "manifest.csv")))    # already sorted by leaderboard_rank
out = HERE / "sample_renders"; out.mkdir(exist_ok=True)


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
    ls.linestyle.color = (0.10, 0.11, 0.13); ls.linestyle.thickness = 1.3


def render(row):
    did = int(row["design_id"])
    pdb = HERE / "pdbs" / f"design_{did:03d}.pdb"
    if not pdb.exists():
        print(f"skip {did}: no pdb"); return
    for o in list(bpy.data.objects):
        bpy.data.objects.remove(o, do_unlink=True)
    mol = mn.Molecule.load(str(pdb)); obj = mol.object
    chain = np.array(mol.named_attribute("chain_id"))
    Color = np.array(mol.named_attribute("Color"))
    binder = HUMAN if row["cohort"].strip().lower() == "human" else AGENT
    Color[chain == 0] = GREY + [TARGET_ALPHA]     # transparent surface reads this alpha
    Color[chain == 1] = binder + [1.0]
    mol.store_named_attribute(Color, "Color")
    # opaque grey material for the TREM2 cartoon so the fold stays visible under the transparent surface
    mat_cart = bpy.data.materials.get("MAT_cart") or bpy.data.materials.new("MAT_cart")
    mat_cart.use_nodes = True
    _b = next((n for n in mat_cart.node_tree.nodes if n.type == "BSDF_PRINCIPLED"), None)
    if _b:
        _b.inputs["Base Color"].default_value = GREY + [1.0]
        _b.inputs["Roughness"].default_value = 0.45
    for ch, nm in (("A", "sa"), ("B", "sb")):
        s = mn.MoleculeSelector(mol); s.chain_id(ch); s.store_selection(nm)
    mol.add_style(mn.StyleSurface(quality=5), color=None, selection="sa")
    # peptide_dssp=False: MN's DSSP recompute truncates the ribbon (see build_template.py); the loaded
    # sheet/loop annotation is complete, so =False draws the full fold with arrows.
    mol.add_style(mn.StyleCartoon(peptide_dssp=False), selection="sa", material=mat_cart)
    mol.add_style(mn.StyleCartoon(), color=None, selection="sb")
    # Cycles honours BSDF alpha natively; set EEVEE's BLENDED too in case the operator switches engines.
    for mat in bpy.data.materials:
        if hasattr(mat, "surface_render_method"):
            mat.surface_render_method = "BLENDED"

    bb = [obj.matrix_world @ mathutils.Vector(c) for c in obj.bound_box]
    ctr = sum(bb, mathutils.Vector()) / 8.0
    rad = max((v - ctr).length for v in bb) or 0.4
    tgt = bpy.data.objects.new("T", None); bpy.context.collection.objects.link(tgt); tgt.location = ctr
    cam = bpy.data.objects.new("C", bpy.data.cameras.new("C")); bpy.context.collection.objects.link(cam)
    cam.location = ctr + mathutils.Vector((0.35, -1, 0.28)).normalized() * rad * 3
    cam.constraints.new("TRACK_TO").target = tgt; bpy.context.scene.camera = cam
    for d, e in [((1.2, -1, 1.2), 26), ((-1.4, -.6, .4), 9), ((0, 1.4, 1), 16)]:
        ld = bpy.data.lights.new("L", "AREA"); ld.energy = e; ld.size = rad
        lo = bpy.data.objects.new("L", ld); bpy.context.collection.objects.link(lo)
        lo.location = ctr + mathutils.Vector(d).normalized() * rad * 2.5
        lo.constraints.new("TRACK_TO").target = tgt
    w = bpy.data.worlds.new("W"); bpy.context.scene.world = w; w.use_nodes = True
    w.node_tree.nodes["Background"].inputs["Strength"].default_value = 0.4
    sc = bpy.context.scene
    setup_cycles(sc, SAMPLES); outline(sc)
    sc.render.film_transparent = True
    sc.render.resolution_x = sc.render.resolution_y = RES
    sc.view_settings.exposure = -0.6
    sc.render.filepath = str(out / f"tile_{did:03d}")
    bpy.ops.render.render(write_still=True)
    print(f"rendered design {did} ({row['cohort']})")


if IDS is not None:
    byid = {int(r["design_id"]): r for r in rows}
    todo = [byid[i] for i in IDS if i in byid]
else:
    todo = rows[:N]
for row in todo:
    render(row)
print(f"DONE — {len(todo)} tiles in {out}; compose with the PIL tiler in FIGURE_GUIDE.md")
