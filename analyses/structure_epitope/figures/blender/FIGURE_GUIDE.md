# Operator guide — Fig 8 (epitope surfaces) & Fig 10 (complex gallery)

Two figures need a human operator + **Blender + Molecular Nodes** (both GPL-3, commercial-use OK;
licence-clean — ChimeraX excluded). This guide gives a **suggested layout** for each, then the
**build steps**. All inputs already exist in this folder (`pdbs/`, `manifest.csv`) and in
`../../results/`. Structures are **predicted (Protenix), unrelaxed — state this in every caption.**

Shared conventions:
- **Cohort colours:** human `#1FE48F`, agent `#30C5F5`. **Target (TREM2):** mid slate-grey `#737D8A`
  (dark enough to read on a white page — not near-white).
- **Epitope ramp** (Fig 8): white `#F2F2F2` → hot `#E8622C` (low→high binder-contact frequency).
- **One common TREM2 frame:** every `pdbs/*.pdb` chain A is pre-superposed (≤0.2 Å target-Cα spread), so a
  single camera gives the same TREM2 orientation everywhere. Keep Fig 8 and Fig 10 in the *same*
  orientation for cross-figure consistency.
- Render with a **transparent film** (Render → Film → Transparent), 300 dpi equivalent, PNG.
- **Engine, quality & how to preview.** `build_template.py` / `render_gallery.py` default to **Cycles on
  GPU** (Metal on macOS / OptiX-CUDA on NVIDIA; CPU fallback) with **OpenImageDenoise**, res 2000 (hero) /
  1200 (tiles), `StyleSurface(quality=5)`, exposure −0.6. Preview in the 3D viewport's **Rendered** shading
  (`Z` → **Rendered**) or a full **F12** render — not **Solid** / **Material Preview**. Cycles honours the
  surface alpha natively, so transparency "just works"; if you switch to **EEVEE Next** for speed, surface
  transparency then needs the material **render method = BLENDED** (the scripts set it — the EEVEE default
  DITHERED renders alpha as stochastic noise that looks broken in the viewport).
- **Contrast on a white page:** target surface = mid slate-grey `#737D8A` + a thin dark outer Freestyle
  outline (external-contour) + exposure −0.6, so the surface reads clearly against white.

---

# Fig 8 — Epitope on the TREM2 IgSF domain

**Story:** binders converge on the **apical** ligand-binding face (the hydrophobic tip + CDR-like loops),
which is the phosphatidylserine site (6B8O) — *not* the distal antibody face (6YYE/6Y6C). Panel **d is
already produced** (`fig08e_epitope_landscape.png`, the per-residue contact-frequency plot); you are making
the three 3-D surface panels **a–c** that give it structural context.

> **Fig 8 = `sample_renders/fig8_partners.png` (the GO version).** Four same-orientation hero panels (TREM2
> never rotated): our binder (design 17), scFv-2 (6YYE), scFv-4 (6Y6C), and phosphatidylserine (6B8O). Each
> crystal partner is superposed onto design_017's TREM2 (rigid CA fit, RMSD ≤ 0.75 Å) and drawn on the
> transparent purple TREM2 surface + fold, with that partner's interaction epitope highlighted — **orange for
> the design/PS epitope, blue for the antibody epitope**; **our binder green, scFv antibodies pink, PS yellow
> spheres**. scFv-2/-4 sit on the distal face, our binder + PS on the apical face — opposite sides, shown
> through the transparent surface. Reproduce:
> ```bash
> uv run python figures/blender/prep_fig8_epitope_sets.py   # -> fig8_epitope_sets.json (design-epitope set)
> uv run python figures/blender/superpose_crystals.py       # -> pdbs/fig8_{6yye,6y6c,6b8o}.pdb (needs refs/*.cif)
> /Applications/Blender.app/Contents/MacOS/Blender --background --python figures/blender/render_fig8_panels.py
> uv run python figures/blender/compose_fig8_panels.py      # -> sample_renders/fig8_partners.png
> ```
> Colours are tuned at the top of `render_fig8_panels.py` (PINK/BLUE_EPI/ACCENT/GREEN/YELLOW) and the legend
> in `compose_fig8_panels.py`. *(The earlier surface epitope-map version — `render_fig8.py`/`compose_fig8.py`
> — was dropped; `prep_fig8_epitope_sets.py` is retained because the partners renderer reads its
> design-epitope set and it still verifies the opposite-face geometry: centroids 153° apart, Jaccard
> design∩PS 0.31 vs design∩scFv 0.00.)*

### Suggested layout
```
┌──────────────┬──────────────┬──────────────┐
│    (a)       │    (b)       │    (c)       │   three TREM2 surfaces,
│  APICAL      │  DISTAL      │  HOTSPOTS    │   IDENTICAL orientation base
│  epitope     │  (rotate     │  labelled    │   (a: apical up;
│  heat-map    │   ~180°)     │  W44 L69 L71 │    b: same object +180°;
│  on surface  │  scFv face   │  F74 L75 L89 │    c: apical up, cartoon)
│  white→hot   │  ≈ cold      │              │
├──────────────┴──────────────┴──────────────┤
│  (d)  per-residue binder-contact frequency  │   ← ALREADY PRODUCED (fig08e), spans full width
│       (Met+19–132; tip peaks at ~L69–F74)   │
└─────────────────────────────────────────────┘
```
- **(a) Apical epitope, coloured by contact frequency.** TREM2 surface, apical face to camera, coloured by
  **binder contact frequency** (white→hot). The apical tip (W44, L69/W70/L71, F74, L89) + CDR-like loops
  (~67–75) light up; everything else stays pale. This is the 3-D twin of panel (d).
- **(b) Distal face, for contrast.** The *same* object rotated ~180° to the β-A/F/G + C–C′ face (the
  scFv-2/scFv-4 epitope). It should read **cold** (near-zero frequency) — designs avoid the antibody face.
  Optionally outline the scFv contact patch (residue set in `../../results/reference_epitope_sets.csv`).
- **(c) Energetic hotspots, labelled.** TREM2 as cartoon (or the same surface, semi-transparent), apical
  view, with the EvoEF2 hotspot residues **W44, L69, L71, F74, L75, L89** labelled — the nominated
  mutagenesis targets. (Top target hotspots by ΔΔG: `../../results/alanine_scan.csv`.)

### Build steps (Blender + Molecular Nodes)

**Prep (once, outside Blender)** — bake the contact frequency into B-factors so MN can colour by a ramp:
```bash
uv run --with gemmi --with pandas python prep_fig8_bfactor.py
# -> writes trem2_target_epitope.pdb (chain A only; B-factor = 100*freq_binder)
#    and prints the tip residues that should light up — sanity-check the offset.
```

1. **Enable Molecular Nodes** — Edit → Preferences → Add-ons → enable *Molecular Nodes* (install its
   dependencies if prompted).
2. **Import** `trem2_target_epitope.pdb` — File → Import → Molecular Nodes. In the node tree add a
   **Style Surface** node (a smooth molecular surface). Set world to a neutral studio HDRI or a 3-point
   light rig; Render → Film → **Transparent**.
3. **Colour by frequency (panel a).** In the MN node tree, after the Style node:
   `Named Attribute (b_factor)` → **Map Range** (0 → max B, i.e. 0→~100, to 0→1) → **Color Ramp**
   (white `#F2F2F2` at 0 → hot `#E8622C` at 1) → feed into the **Set Color** node. Preview in Material/
   Rendered view; the apical tip should glow.
4. **Camera + orientation.** Add a camera, frame the **apical** face (tip toward camera). **Do not move
   the camera again** — b/c reuse the same frame. Render → `renders/fig8a_apical.png`.
5. **Panel (b).** Duplicate the object (or add a keyframe/second camera); rotate the *object* 180° about
   the vertical axis to show the distal face. Same ramp — it should render cold. Optionally add a thin
   coloured outline over the scFv residue set (select-by-residue → separate Style → flat colour).
   Render → `renders/fig8b_distal.png`.
6. **Panel (c).** Switch the Style to **Cartoon** (or make the surface ~40 % transparent), apical view.
   Add **text labels** for W44, L69, L71, F74, L75, L89: for each, place a small Blender *Text* object at
   that residue's Cα (read the Cα xyz from `trem2_target_epitope.pdb`) with a leader line.
   **Numbering:** the PDB is **construct-numbered** (UniProt = construct + 17), so the manuscript's UniProt
   labels sit at PDB residues **W44→27, L69→52, L71→54, F74→57, L75→58, L89→72** (confirmed by
   `prep_fig8_bfactor.py`'s sanity line: the same residues carry the highest baked B-factor).
   Render → `renders/fig8c_hotspots.png`.
7. **Assemble** a + b + c with the produced **d** (`../fig08e_epitope_landscape.png`) in Inkscape/
   Illustrator (or a small PIL script), three surfaces in a row over the full-width heat-map, panel
   letters a–d. Save `../fig08_epitope.png`.

---

# Fig 10 — Structural complex gallery

**Story:** what the designed binder–TREM2 complexes look like, and that the competition's in-silico
(ipSAE) ranking does **not** track affinity. Full render mechanics are in `README.md`; this section adds
the **layout** and consolidates the steps.

### Suggested layout
```
┌───────────────────────────────┬───────────────────┐
│                               │      (b)          │
│   (a)  LEADERBOARD GRID        │  close-up:        │
│        top-16 by ipSAE, 4×4    │  design 017       │
│        (== Fig 1 order)        │  (1.11 nM)        │
│   tiles: TREM2 grey + binder   │  interface +      │
│   in cohort colour; tile       │  EvoEF2 hotspots  │
│   border/label = is_hit & K_D  ├───────────────────┤
│   (design 13 = #1;             │      (c)          │
│    design 017 = only #24 →     │  37 binders       │
│    the "ipSAE ≠ affinity"      │  superposed on    │
│    point)                      │  one TREM2, by    │
│                               │  cohort           │
└───────────────────────────────┴───────────────────┘
```
- **(a) Leaderboard grid** — top-N in `leaderboard_rank` (ipSAE) order from `manifest.csv` (already
  sorted). Pick N for a clean tile (top-16 → 4×4). Colour each binder by cohort; encode **is_hit** and
  **K_D** in the tile border/label so the grid shows ipSAE interleaving binders and non-binders.
- **(b) Best-design close-up** — `design_017` (1.11 nM), zoomed to the interface, with the EvoEF2 hotspot
  side-chains shown/labelled (`../../results/alanine_scan.csv`). Pairs with (a): the tightest binder ranks
  only #24.
- **(c) All 37 binders superposed** — every binder's chain B imported over one `design_017` target (all
  share the common frame), coloured by cohort — the "one apical epitope" picture in 3-D.

### Build steps (Blender + Molecular Nodes)

1. **Template — already built (purple, transparent, fold-through).** `template.blend` is pre-generated by
   `build_template.py` (verified on Blender 5.0.1 + Molecular Nodes; see `template_preview.png`). It imports
   `pdbs/design_017.pdb` and renders **TREM2 as a transparent purple surface** (apical epitope, freq ≥ 0.5,
   highlighted orange) **over an opaque purple cartoon of TREM2's own secondary structure**, so the β-sandwich
   fold shows *through* the surface; **chain B is a human-green Cartoon**; camera + 3-point lights +
   transparent film. Key trick: the surface reads a **per-atom `Color` attribute** (RGBA, A = opacity =
   `TARGET_ALPHA` 0.22) while the cartoon uses a **flat opaque material** (they can't share one alpha). Knobs
   at the top of `build_template.py` — **`COLOUR_BY_EPITOPE`**, **`TARGET_ALPHA`** — toggle the epitope
   highlight / opacity. **Open it and refine:** the exact **apical** camera angle, colours (edit the `Color`
   attribute / `PURPLE_*` constants), then **lock the camera.** Regenerate any time:
   `/Applications/Blender.app/Contents/MacOS/Blender --background --python build_template.py`.
   - **Two pitfalls the template already handles — keep them if you edit the styles:**
     - **`StyleCartoon(peptide_dssp=False)`, not `True`.** MN's on-the-fly DSSP recompute silently
       **truncates the target ribbon to ~1/3 of the domain** (verified: 1044 vs 5530 mesh verts). The
       loaded structure already carries a complete sheet/loop annotation, so `=False` draws the full
       β-sandwich *with* arrows. `=True` is the cause of a "missing / incomplete secondary structure" look.
     - **Transparency only shows in a real render.** View in **Rendered** shading (`Z` → Rendered) or **F12**
       with the **Cycles** engine — **Solid / Material-Preview shading shows the surface opaque.**
   *(This transparent-surface-over-cartoon is the Fig 8a look. For a continuous epitope heat-map instead of
   the ≥0.5 cutoff, import `trem2_target_epitope.pdb` and colour by `b_factor`. The Fig 10 gallery uses the
   **same transparent-surface-over-fold** look in neutral grey — `render_gallery.py` — so each binder reads
   through the surface while cohort colour stays on the binder.)*
2. **Panel (a) renders.** Use the **working, verified** batch renderer `render_gallery.py` (renders in the
   neutral cohort look — **transparent** grey surface over an opaque grey fold cartoon + green/blue binder —
   to `sample_renders/tile_<id>.png`). Two selection modes:
   ```bash
   # top-N by ipSAE leaderboard_rank
   /Applications/Blender.app/Contents/MacOS/Blender --background --python render_gallery.py -- --n 16
   # an explicit set (e.g. the tightest-Kd designs, which are scattered down the ipSAE ranking)
   /Applications/Blender.app/Contents/MacOS/Blender --background --python render_gallery.py -- --ids 17,74,102,34
   ```
   Then tile with PIL (no Blender). **Two committed sample grids** (both 3×2, same label style):
   - `sample_renders/fig10a_sample_grid.png` — **ipSAE** leaderboard, ranks #1–6 (`ids=[13,14,16,1,35,132]`).
   - `sample_renders/fig10_kd_grid.png` — **Kd** leaderboard, tightest 6 binders (`ids=[17,13,74,132,102,34]`).
     Tells §5's story: the tightest binder `design_017` (1.11 nM) sits at ipSAE rank **24**; `design_074`
     (1.91 nM) at rank **78**.

   Label sits ABOVE each tile: **white box, cohort-coloured outline** (human green / agent blue), two lines —
   near-black `Rank #N   Kd = X nM` over a **cohort-coloured** `human|agent · team` (same green/blue,
   darkened from the brand hues `#1FE48F`/`#30C5F5` so they read on white). **ipSAE non-binders show
   `non-binder`, not a Kd** — high-ipSAE designs that never bound (verified `is_hit=False`, `n_with_kd=0` in
   `data/designs.csv`); do NOT invent a Kd. The Kd grid is all binders, so every tile has a Kd.
   ```python
   from PIL import Image, ImageDraw, ImageFont; import csv; from pathlib import Path
   H=Path("figures/blender"); rows={int(r["design_id"]):r for r in csv.DictReader(open(H/"manifest.csv"))}
   W,pad,cols=600,12,3; INK=(20,20,22); HUM=(19,139,87); AGT=(29,120,149)  # outline + 2nd line = cohort colour
   ft=lambda s,b: ImageFont.truetype(f"/System/Library/Fonts/Supplemental/Arial{' Bold' if b else ''}.ttf",s)
   F1,F2=ft(34,True),ft(28,False); pin,gap=16,8; h1,h2=F1.getbbox("Ag")[3],F2.getbbox("Ag")[3]; lh=pin+h1+gap+h2+pin
   kdf=lambda v:(lambda k: f"{k:.0f}" if k>=100 else f"{k:.1f}" if k>=10 else f"{k:.2f}")(float(v))
   def build(ids,out,mode):
       cellH=lh+W; rN=(len(ids)+cols-1)//cols
       cv=Image.new("RGB",(cols*W+(cols+1)*pad, rN*cellH+(rN+1)*pad),(255,255,255)); d=ImageDraw.Draw(cv)
       for i,did in enumerate(ids):
           r=rows[did]; x=pad+(i%cols)*(W+pad); y=pad+(i//cols)*(cellH+pad); v=r["kd_nM"].strip()
           col=HUM if r["cohort"].strip().lower()=="human" else AGT
           if mode=="ipsae": rk=int(float(r["leaderboard_rank"])); l1=f"Rank #{rk}    Kd = {kdf(v)} nM" if v else f"Rank #{rk}    non-binder"
           else:             l1=f"Kd rank #{i+1}    Kd = {kdf(v)} nM"
           d.rectangle([x,y,x+W,y+lh],fill=(255,255,255)); d.rectangle([x,y,x+W-1,y+lh-1],outline=col,width=4)
           d.text((x+pin+4,y+pin),l1,fill=INK,font=F1); d.text((x+pin+4,y+pin+h1+gap),f"{r['cohort']} · {r['team']}",fill=col,font=F2)
           im=Image.open(H/f"sample_renders/tile_{did:03d}.png").convert("RGBA").resize((W,W))
           bg=Image.new("RGBA",im.size,(255,255,255,255)); bg.alpha_composite(im); cv.paste(bg.convert("RGB"),(x,y+lh))
       cv.save(out)
   build([13,14,16,1,35,132], H/"sample_renders/fig10a_sample_grid.png", "ipsae")
   build([17,13,74,132,102,34], H/"sample_renders/fig10_kd_grid.png", "kd")
   ```
   (`render_headless.py`/`make_grid.py` remain the original single-design skeletons; `render_gallery.py`
   supersedes them for the grid.)
3. **Panel (b) — BUILT.** `render_fig10b.py` renders the design_017 interface close-up (hero preset) with
   the top EvoEF2 hotspot TREM2 side chains as red sticks, and writes `fig10b_labels.json` (projected 2-D
   label positions). Hotspots shown = top 6 by ΔΔG (`alanine_scan.csv`): W44 9.7 · L72 5.7 · L71 5.6 ·
   L69 4.6 · L75 4.3 · F74 4.2 kcal/mol. *(Note: the draft's nominated set swaps L72→L89; L89 is a named
   tip residue but only 1.06 kcal/mol here — decide which set the caption should show.)*
4. **Panel (c) — BUILT.** `prep_fig10c.py` builds one combined PDB (design_017 TREM2 + all 37 hit binders'
   chain B, cohort baked into B-factor); `render_fig10c.py` renders the fan coloured by cohort
   (human green / agent blue) on the transparent purple TREM2.
5. **Assemble — BUILT.** `compose_fig10.py` stacks a (leaderboard grid, full width) over b + c, draws the
   panel-b hotspot labels on a fanned ring, and writes `sample_renders/fig10_gallery.png`. Full pipeline:
   ```bash
   uv run python figures/blender/prep_fig10c.py
   BL=/Applications/Blender.app/Contents/MacOS/Blender
   $BL --background --python figures/blender/render_fig10b.py
   $BL --background --python figures/blender/render_fig10c.py
   uv run python figures/blender/compose_fig10.py     # -> sample_renders/fig10_gallery.png
   ```

### Two reconciliations before locking Fig 10 (from `README.md`)
- **Confirm panel (a) order == the manuscript Fig 1 leaderboard.** Our `selection_rank` and the draft
  disagree on **design 79** (rank 36 here vs "#7" in the rules-compliance note) — resolve which ranking
  Fig 1 uses. If the grid extends past ~rank 35, flag design 79 as the scFv-4 literature copy.
- `README.md` still calls the affinity-ranking point "**§5**"; after the review-session merge it is now
  **§3** (the *"learned scores don't rank affinity"* result).

---

## Caption stubs (drop into the manuscript)
- **Fig 8.** *Designs converge on the apical ligand-binding epitope of the TREM2 IgSF domain.* (a) TREM2
  surface coloured by binder contact frequency (predicted Protenix complexes, n=… screened); (b) the
  distal antibody (scFv) face for contrast; (c) EvoEF2 energetic hotspots (W44, L69, L71, F74, L75, L89);
  (d) per-residue contact frequency.
- **Fig 10.** *Structural gallery of designed TREM2 binders.* (a) top-16 by ipSAE (competition ranking;
  binders and non-binders interleaved — design 017, the tightest at 1.11 nM, ranks only #24); (b) design
  017 interface with energetic hotspots; (c) the 37 experimental binders superposed on TREM2, by cohort.
  All complexes predicted (Protenix), unrelaxed.
