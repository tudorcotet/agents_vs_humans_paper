"""Compose the 4-panel hero-preset Fig 8 (2x2), partner titles + shared legend + orientation caption.
Pure PIL.  Run:  uv run python figures/blender/compose_fig8_panels.py  ->  sample_renders/fig08_partners.png
"""
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

H = Path(__file__).resolve().parent
SR = H / "sample_renders"
W, pad = 820, 26
INK, SUB = (22, 22, 26), (95, 97, 105)
PURPLE, ORANGE, BLUE = (153, 122, 199), (232, 97, 43), (33, 112, 230)
GREEN, PINK, YELLOW = (31, 178, 118), (232, 107, 153), (233, 193, 36)

def font(sz, bold=False):
    p = "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else \
        "/System/Library/Fonts/Supplemental/Arial.ttf"
    return ImageFont.truetype(p, sz)

FT_T, FT_S, FT_L, FT_C = font(34, True), font(24, False), font(25, True), font(23, False)

def on_white(p):
    im = Image.open(p).convert("RGBA").resize((W, W))
    bg = Image.new("RGBA", im.size, (255, 255, 255, 255)); bg.alpha_composite(im)
    return bg.convert("RGB")

panels = [
    ("fig8p_hero.png",  "a   Our binder — design 17",         "human · MRAZS · Kd 1.11 nM · apical face"),
    ("fig8p_scfv2.png", "b   scFv-2 — PDB 6YYE",              "Szykowska 2021 · distal antibody epitope"),
    ("fig8p_scfv4.png", "c   scFv-4 — PDB 6Y6C",              "Szykowska 2021 · distal antibody epitope"),
    ("fig8p_ps.png",    "d   Phosphatidylserine — PDB 6B8O",  "Sudom 2018 · apical ligand site"),
]
h_head = 90
cw = 2 * W + 3 * pad
band_h = 150
ch = pad + 2 * (h_head + W) + pad + band_h
cv = Image.new("RGB", (cw, ch), (255, 255, 255))
d = ImageDraw.Draw(cv)

for i, (img, title, sub) in enumerate(panels):
    col, row = i % 2, i // 2
    x = pad + col * (W + pad)
    y = pad + row * (h_head + W + pad)
    d.text((x, y), title, fill=INK, font=FT_T)
    d.text((x + 42, y + 46), sub, fill=SUB, font=FT_S)
    cv.paste(on_white(SR / img), (x, y + h_head))

# legend
items = [(PURPLE, "TREM2 (surface + fold)"), (ORANGE, "design / PS epitope"), (BLUE, "antibody epitope"),
         (GREEN, "design binder"), (PINK, "scFv antibody"), (YELLOW, "PS lipid")]
chip, cgap, igap = 28, 10, 40
widths = [chip + cgap + d.textlength(t, font=FT_L) for _, t in items]
x = (cw - (sum(widths) + igap * (len(items) - 1))) / 2
ly = pad + 2 * (h_head + W) + pad + 24
for (col, label), wdt in zip(items, widths):
    d.rounded_rectangle([x, ly + 2, x + chip, ly + 2 + chip], radius=6, fill=col)
    d.text((x + chip + cgap, ly), label, fill=INK, font=FT_L)
    x += wdt + igap

cap1 = "Same TREM2 orientation in every panel — partners superposed onto design_017 (fit RMSD ≤ 0.75 Å)."
cap2 = "scFv-2/-4 bind the distal face; our binder and PS occupy the apical ligand face — opposite sides."
d.text(((cw - d.textlength(cap1, font=FT_C)) / 2, ly + 56), cap1, fill=SUB, font=FT_C)
d.text(((cw - d.textlength(cap2, font=FT_C)) / 2, ly + 88), cap2, fill=SUB, font=FT_C)

out = SR / "fig08_partners.png"
cv.save(out)
print("wrote", out.name, cv.size)
