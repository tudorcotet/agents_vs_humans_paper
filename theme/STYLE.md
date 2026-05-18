# Visual language — agents-vs-humans (May 2026)

The house style for figures in this paper. Light canvas, navy ink, single
cyan accent. The tokens below match the seven hand-authored HTML figures
in `figures/blog/`.

## Editorial canon (light canvas)

For paper and blog figures (graphical abstract, social cards, hero), use
these tokens. They're locked across every figure so the paper reads as
one piece of work.

```css
:root {
  --text:        #111827;
  --text-mute:   #4B5563;
  --text-dim:    #6B7280;
  --line:        #D1D5DB;
  --line-soft:   #E5E7EB;
  --brand-black: #142933;   /* navy ink, NOT pure black */
  --brand-cyan:  #30C5F5;   /* the only accent that ships */
  --cyan-soft:   #9EDFFF;
}
```

Body text and labels: `--brand-black` on white. One accent word per hero
in `--brand-cyan`. Numerals tabular (Geist Mono / IBM Plex Mono),
centered.

**Card pattern (chassis for every editorial figure):**
- `background: #FFFFFF`
- `border: 1px solid rgba(15,20,25,0.08)`
- `border-radius: 16px`
- `box-shadow: 0 24px 56px rgba(15,20,25,0.07)`
- `padding: 44px 56px`

**Card-idx row (the signature pattern):** `01  PANEL NAME ────────────`
where `01` is in cyan, `PANEL NAME` is in brand-black, and the hairline
rule on the right is cyan at 45% opacity, 2px tall. Both labels in GT
Pressura Extended, 28px, letter-spacing 0.04em, uppercase.

**Big-number stat block:** GT Pressura Extended, 80px (lead) / 44px
(secondary), `letter-spacing: -0.04em`, `line-height: 0.92`. Lead number
cyan, secondary numbers brand-black. Always paired with a `.lab` (12px
Pressura, 0.18em tracking, uppercase, `--text-mute`) above or beside.

**Title gradient (one accent word per hero, max):**
`linear-gradient(90deg, #9EDFFF 0%, #36B7F6 50%, #142933 100%)` clipped
to text. Light → dark, cyan family only. Solid `--brand-cyan` is also
acceptable for the accent word.

## Cohort glyphs (load-bearing across every figure)

| Cohort | Color | Hex | Marker glyph (matplotlib) | Notes |
|---|---|---|---|---|
| Human (10 teams, 81 designs) | Cyan  | `#30C5F5` | hexagon (`h`) | Primary brand colour. Stays cyan everywhere. |
| Agent (6 teams, 60 designs)  | Amber | `#FFB547` | diamond (`D`) | Warn-amber held in reserve elsewhere; reused here because cohort is the load-bearing categorical of the study. |

Locked in `theme/palettes.json` under `agent_vs_human`.

## Chart conventions (editorial light canvas)

- Axis spines: 0.5pt hairline, brand-black, top + right hidden.
- Tick marks: 2.5pt length, 0.5pt width, outside, brand-black.
- Bars: solid fill, no edge, opacity 0.85–1.0. No gradients.
- Markers: solid fill, opacity 0.80, no outline. Add a brand-black
  outline (0.5pt) only when the marker sits on a same-colour fill.
- Callouts: 0.5pt brand-black leader line + small mono label
  (Geist Mono 7.5pt). No bounding box. Pull right unless that side is
  occupied.
- KD readouts: Geist Mono, 11–12px, tnum + lnum on. Render
  `K_D = 1.11 nM` with `=` (not `:`) and a non-breaking space before
  the unit.

## Typography

Three faces, three jobs:

1. **GT Pressura Extended Regular** (weight 400 only) — display /
   headings. Wide, condensed-feel. Tracking near zero. Sizes calibrated
   low: title-xl 60px, title-lg 50px, title-md 42px (line-height
   1.02–1.06, tracking -1.0 to -0.6px).
2. **Geist (Regular / Medium / Bold)** — body, UI, chart labels. The
   default. Set `font-feature-settings: "ss01", "tnum"`. Thesis
   paragraph: 20px / 1.45.
3. **Geist Mono / IBM Plex Mono** — K_D readouts, sequence IDs, axis
   labels, small uppercase technical metadata. Mono KD value: 22px /
   1.1 / -0.4px / weight 500. **Never** set thesis in Pressura.

Font files live in `theme/fonts/` and `figures/blog/fonts/`. See
`theme/fonts/README.md`.

## Layout idioms

Three patterns. Reuse them.

1. **Text-left / data-right split.** Title + lockup left ~60%, data
   panel right ~40% on a soft cyan backlight. Grid
   `grid-template-columns: 1fr 480px` (or `1fr 600px` for 16:9).
2. **Card on dot grid + directional cyan glow.** Dot grid 32px at 5%
   opacity. Cyan radial wash always top-right, never centred. Card
   solid white, 16px radius, 1px hairline border, `0 24px 56px
   rgba(15,20,25,0.07)` shadow.
3. **Sensorgram as visual proof.** 3 concentrations (200/100/50 nM),
   baseline / association / dissociation phases, dashed phase
   boundaries, mono axis labels. Cyan stroke gradient `#00D9FF` to
   `#33C4FF`, faint cyan fill under the top trace. The data motif —
   never a glowing helix or floating molecule.

## Forbidden moves

AI-biotech defaults. Don't.

- Glowing DNA helices, floating molecules, orbiting atoms, stock-photo
  futuristic lab imagery.
- Gradient text on the hero title. Solid ink with one accent word in
  solid `--brand-cyan`. The cyan gradient is reserved for one specific
  accent word per hero, cyan family only.
- Glassmorphism / `backdrop-filter: blur()`.
- Lavender / purple-pink AI palette.
- Inter, system-ui, or "sans-serif". Body is Geist; display is Pressura.
- "AI" sparkle badge above the title.
- Centred radial glow blob behind the title. Glow is directional,
  top-right.
- Em-dashes inside chart labels or captions. Use commas, periods, or
  spaced hyphens.
