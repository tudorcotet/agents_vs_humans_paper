"""Compose the top-10-by-Kd hit-binder gallery (SI Fig S7), reusing the Fig 10 tiler preset.

Pure PIL (no Blender): lays out the ten tightest hit binders (by measured K_D, tightest first) as a 5x2 grid
of the render_gallery.py tiles (sample_renders/tile_<id>.png), with the same cohort-coloured label boxes as
the main Fig 10 gallery (human green / agent blue). Design 79 — a near-copy of anti-TREM2 scFv-4, i.e. the
inadvertent positive control discussed in the Methods rules-compliance note, not a de novo design — is flagged
with a dagger + amber box + footnote so the figure never reads it as a clean de novo hit.

The ranking is computed from manifest.csv (is_hit & kd_nM), so it stays in sync with the data.
Requires the ten tiles to exist (render first):
  BL=/Applications/Blender.app/Contents/MacOS/Blender
  $BL --background --python figures/blender/render_gallery.py -- --ids 17,13,74,132,102,34,79,35,45,103
  uv run python figures/blender/compose_fig10_kd.py   # -> sample_renders/figS7_kd_gallery_top10.png
"""
import csv
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

H = Path(__file__).resolve().parent
SR = H / "sample_renders"
rows = {int(r["design_id"]): r for r in csv.DictReader(open(H / "manifest.csv"))}

# ten tightest hit binders by measured Kd, data-driven from the manifest
binders = [(did, float(r["kd_nM"])) for did, r in rows.items()
           if r["kd_nM"].strip() and r.get("is_hit", "").strip().lower() == "true"]
TOP = [did for did, _ in sorted(binders, key=lambda t: t[1])[:10]]
FLAG = {79: "near-copy of anti-TREM2 scFv-4 (Methods) — inadvertent positive control, not a de novo design"}

W, pad, cols = 600, 12, 5
INK = (20, 20, 22); HUM = (19, 139, 87); AGT = (29, 120, 149); WARN = (176, 96, 20)
ft = lambda s, b: ImageFont.truetype(f"/System/Library/Fonts/Supplemental/Arial{' Bold' if b else ''}.ttf", s)
F1, F2, FT, FF = ft(34, True), ft(28, False), ft(34, True), ft(24, False)
pin, gap = 16, 8
h1, h2 = F1.getbbox("Ag")[3], F2.getbbox("Ag")[3]
lh = pin + h1 + gap + h2 + pin
kdf = lambda v: (lambda k: f"{k:.0f}" if k >= 100 else f"{k:.1f}" if k >= 10 else f"{k:.2f}")(float(v))

cellH = lh + W
rN = (len(TOP) + cols - 1) // cols
top_h, foot_h = 62, 66
Wc = cols * W + (cols + 1) * pad
Hc = top_h + rN * cellH + (rN + 1) * pad + foot_h
cv = Image.new("RGB", (Wc, Hc), (255, 255, 255))
d = ImageDraw.Draw(cv)

d.text((pad, 14), "Fig S7. Top 10 hit binders by measured K_D (tightest first) — cohort-coloured "
                  "(human green / agent blue), on the shared TREM2 frame.", fill=INK, font=FT)

for i, did in enumerate(TOP):
    r = rows[did]; v = r["kd_nM"].strip()
    x = pad + (i % cols) * (W + pad)
    y = top_h + pad + (i // cols) * (cellH + pad)
    col = HUM if r["cohort"].strip().lower() == "human" else AGT
    dagger = " †" if did in FLAG else ""
    l1 = f"Kd rank #{i+1}{dagger}    Kd = {kdf(v)} nM"
    d.rectangle([x, y, x + W, y + lh], fill=(255, 255, 255))
    d.rectangle([x, y, x + W - 1, y + lh - 1], outline=(WARN if did in FLAG else col), width=4)
    d.text((x + pin + 4, y + pin), l1, fill=INK, font=F1)
    d.text((x + pin + 4, y + pin + h1 + gap), f"{r['cohort']} · {r['team']}", fill=col, font=F2)
    im = Image.open(SR / f"tile_{did:03d}.png").convert("RGBA").resize((W, W))
    bg = Image.new("RGBA", im.size, (255, 255, 255, 255)); bg.alpha_composite(im)
    cv.paste(bg.convert("RGB"), (x, y + lh))

fy = top_h + rN * cellH + (rN + 1) * pad + 10
d.text((pad, fy), f"† design 79 — {FLAG[79]}.", fill=WARN, font=FF)

out = SR / "figS7_kd_gallery_top10.png"
cv.save(out)
print("wrote", out.name, cv.size, "| ids:", TOP)
