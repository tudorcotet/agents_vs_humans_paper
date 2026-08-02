"""Assemble the full multi-panel Fig 10:
  a  leaderboard gallery (fig10a_sample_grid.png)         — top, full width
  b  best-design (17) interface + EvoEF2 hotspot labels   — bottom-left
  c  all 37 hit binders superposed, by cohort             — bottom-right
Hotspot labels for b are placed from projected pixel coords (fig10b_labels.json), fanned out from the
cluster centroid so they don't overlap.  Run: uv run python figures/blender/compose_fig10.py
  -> sample_renders/fig10_gallery.png
"""
import json, math
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

H = Path(__file__).resolve().parent
SR = H / "sample_renders"
INK, SUB = (22, 22, 26), (95, 97, 105)
HOTCOL = (200, 40, 30)
DDG = {"W44": 9.7, "L72": 5.7, "L71": 5.6, "L69": 4.6, "L75": 4.3, "F74": 4.2}

def font(sz, bold=False):
    p = "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else \
        "/System/Library/Fonts/Supplemental/Arial.ttf"
    return ImageFont.truetype(p, sz)

FT_T, FT_S, FT_HL = font(38, True), font(25, False), font(26, True)

def white(p, w):
    im = Image.open(p).convert("RGBA")
    im = im.resize((w, int(im.height * w / im.width)))
    bg = Image.new("RGBA", im.size, (255, 255, 255, 255)); bg.alpha_composite(im)
    return bg.convert("RGB")

def halo(d, xy, text, fnt, fill):
    x, y = xy
    for dx in (-2, -1, 0, 1, 2):
        for dy in (-2, -1, 0, 1, 2):
            d.text((x + dx, y + dy), text, font=fnt, fill=(255, 255, 255))
    d.text((x, y), text, font=fnt, fill=fill)

CW = 1860
pad = 14
grid = white(SR / "fig10a_sample_grid.png", CW)
Wbc = (CW - pad) // 2
pb = white(SR / "fig10b_design017.png", Wbc)
pc = white(SR / "fig10c_superposed.png", Wbc)

head = 58
y_a = head
y_bc_title = y_a + grid.height + 26
y_bc = y_bc_title + head
cap_h = 96
CH = y_bc + Wbc + cap_h
cv = Image.new("RGB", (CW, CH), (255, 255, 255))
d = ImageDraw.Draw(cv)

# a
d.text((0, 6), "a   Structural leaderboard — top 6 by ipSAE", fill=INK, font=FT_T)
cv.paste(grid, (0, y_a))
# b / c titles
d.text((0, y_bc_title - 4), "b   Best design (17): interface + hotspots", fill=INK, font=FT_T)
d.text((Wbc + pad, y_bc_title - 4), "c   All 37 hit binders on TREM2, by cohort", fill=INK, font=FT_T)
cv.paste(pb, (0, y_bc))
cv.paste(pc, (Wbc + pad, y_bc))

# hotspot labels on panel b — place on an even ring around the cluster centroid (guaranteed no overlap)
lab = json.load(open(H / "fig10b_labels.json"))
sc = Wbc / lab["res"]
pts = {n: (xy[0] * sc, y_bc + xy[1] * sc) for n, xy in lab["labels"].items()}
cx = sum(p[0] for p in pts.values()) / len(pts)
cy = sum(p[1] for p in pts.values()) / len(pts)
order = sorted(pts, key=lambda n: math.atan2(pts[n][1] - cy, pts[n][0] - cx))
base = math.atan2(pts[order[0]][1] - cy, pts[order[0]][0] - cx)
ring = 168
for i, n in enumerate(order):
    px, py = pts[n]
    ang = base + i * (2 * math.pi / len(order))
    ax, ay = cx + ring * math.cos(ang), cy + ring * math.sin(ang)
    d.line([(px, py), (ax, ay)], fill=(90, 90, 94), width=2)
    d.ellipse([px - 6, py - 6, px + 6, py + 6], fill=HOTCOL, outline=(255, 255, 255), width=2)
    tw = d.textlength(n, font=FT_HL)
    tx = ax - tw / 2 if abs(math.cos(ang)) < 0.5 else (ax + 4 if math.cos(ang) > 0 else ax - tw - 4)
    halo(d, (tx, ay - 14), n, FT_HL, HOTCOL)

# captions (panel b wraps to two lines so it stays within its half-width)
cy0 = y_bc + Wbc + 8
hs = " · ".join(f"{k} {v}" for k, v in DDG.items())
d.text((0, cy0), "Red side chains = top TREM2 hotspots (EvoEF2 ΔΔG, kcal/mol):", fill=SUB, font=FT_S)
d.text((0, cy0 + 30), hs, fill=SUB, font=FT_S)
d.text((Wbc + pad, cy0), "37 hits (25 human · 12 agent)", fill=SUB, font=FT_S)
d.text((Wbc + pad, cy0 + 30), "converge on the apical ligand face.", fill=SUB, font=FT_S)

out = SR / "fig10_gallery.png"
cv.save(out)
print("wrote", out.name, cv.size)
