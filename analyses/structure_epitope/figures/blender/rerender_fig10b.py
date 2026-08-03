"""Re-render Fig 10b from an operator-rotated fig10b.blend, re-projecting the hotspot labels for the NEW
camera angle. The camera in fig10b.blend is a free object (no track-to), so the operator can rotate it in the
GUI, save, and hand it back.

    /Applications/Blender.app/Contents/MacOS/Blender --background fig10b.blend --python rerender_fig10b.py

Writes sample_renders/fig10b_design017.png + fig10b_labels.json (matched to the rotated camera). Then rerun
compose_fig10.py to rebuild fig10_gallery.png. Reads fig10b_hotspots.json (local hotspot centroids + MN
object name) written by render_fig10b.py — keep it next to this script.
"""
import json
import bpy
import numpy as np
import mathutils
from pathlib import Path
from bpy_extras.object_utils import world_to_camera_view

HERE = Path(__file__).resolve().parent
OUT = HERE / "sample_renders"; OUT.mkdir(exist_ok=True)
meta = json.load(open(HERE / "fig10b_hotspots.json"))


def setup_cycles(scene):
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
    scene.cycles.use_denoising = True
    for den in ("OPENIMAGEDENOISE", "OPTIX"):
        try: scene.cycles.denoiser = den; break
        except Exception: continue


sc = bpy.context.scene
setup_cycles(sc)                                   # re-enable GPU/denoise on this machine
cam = sc.camera
obj = bpy.data.objects.get(meta["object"]) or next(
    (o for o in bpy.data.objects if o.type == "MESH" and o.data and "res_id" in o.data.attributes), None)
M3 = np.array(obj.matrix_world.to_3x3()); T = np.array(obj.matrix_world.translation)
RES = sc.render.resolution_x
bpy.context.view_layer.update()

labels = {}
for name, loc in meta["hotspots"].items():
    world = np.array(loc) @ M3.T + T                # local centroid -> current world position
    ndc = world_to_camera_view(sc, cam, mathutils.Vector(world))
    labels[name] = [round(ndc.x * RES, 1), round((1 - ndc.y) * RES, 1), round(ndc.z, 3)]
json.dump({"res": RES, "labels": labels}, open(HERE / "fig10b_labels.json", "w"), indent=1)

sc.render.filepath = str(OUT / "fig10b_design017")
bpy.ops.render.render(write_still=True)
print("re-rendered fig10b from rotated blend | labels", labels)
