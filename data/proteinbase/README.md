# `data/proteinbase/` — public data package

Mirrored from ProteinBase's public collection for this hackathon. Every
file here is also live at `https://proteinbase-pub.t3.storage.dev/`. We
keep a local snapshot so the analyses run offline and don't break when a
third-party CDN goes down.

## Layout

```
data/proteinbase/
├── boltz2/             # 100 × Boltz-2 complex predictions (.cif)
├── esmfold/            # 100 × ESMFold single-chain predictions (.cif)
├── pae/                # 100 × PAE matrices (.json, residue × residue)
├── images/             # 99 × ESMFold stylised renders (.png)
└── sensorgrams/        # 215 × raw kinetic curves (.json)
                        #   filename: design_NNN_repNN_{bli,spr}.json
                        #   193 SPR + 22 BLI traces across the 100 screened designs
```

## Coverage

The 100 screened designs (top-100 by ipSAE) are the universe here. The
41 non-screened designs are absent — they didn't make the wet lab and
have no kinetic data. Use `submitted_to_lab=True` to filter.

| artifact | files | total size |
|---|---|---|
| Boltz-2 CIF       | 100 | 14 MB |
| ESMFold CIF       | 100 | 5 MB  |
| PAE JSON          | 100 | 88 MB |
| Stylised PNG      |  99 | 27 MB |
| Kinetic curves    | 215 | 27 MB |
| **total**         | 614 | **161 MB** |

## How to access from Python

The local paths are columns in `data/designs.csv`:

```python
from scripts.utils import load_designs

df = load_designs(only_screened=True)
print(df[['design_id', 'name', 'pb_boltz2_cif', 'pb_pae_json']].head())
#  design_id          name                pb_boltz2_cif           pb_pae_json
#          1  1_NovoFy_...   data/proteinbase/boltz2/design_001.cif   data/proteinbase/pae/design_001.json
```

For sensorgrams (multiple per design), glob the folder:

```python
from pathlib import Path
import json

rep_files = sorted(Path("data/proteinbase/sensorgrams").glob("design_017_rep*.json"))
for f in rep_files:
    curve = json.loads(f.read_text())
    # curve has timeseries data per concentration step
```

## What's in each artifact

### `boltz2/design_NNN.cif`

mmCIF complex prediction from Boltz-2. Chain A is TREM2 (residues 19–174
of UniProt Q9NZC2 plus the linker + 10×His), chain B is the designed
binder. Per-atom pLDDT in the B-factor column.

### `esmfold/design_NNN.cif`

mmCIF single-chain prediction from ESMFold. The binder only — no target.
Use this for binder-only stability / pLDDT analyses.

### `pae/design_NNN.json`

Predicted Aligned Error matrix from Boltz-2, residue × residue. Used to
compute ipSAE, ipTM, and interface-confidence statistics. Each value is
a predicted error in Angstroms.

### `images/design_NNN.png`

A pre-rendered cartoon of the ESMFold structure. Used in the leaderboard
figures and the graphical abstract.

### `sensorgrams/design_NNN_repNN_{bli,spr}.json`

Raw kinetic-curve traces from the wet-lab assay. Each file is one
replicate of one design on one instrument (`bli` for Gator BLI, `spr`
for Carterra SPR). Schema (typical):

```json
{
  "trace_01": {
    "raw": {"t": [...], "y": [...]},
    "concentration": 0.0,
    "control": false,
    "aggregated": false,
    "virtual": false,
    "fit": {}
  },
  "trace_02": { ... },
  ...
}
```

Each top-level `trace_NN` key is one concentration step in the
concentration series for that replicate. `raw.t` / `raw.y` are time
(s) and response (nm) arrays. `fit` carries per-trace fit metadata
where available.

Use these when you need the raw response curves — for example, to
re-fit, to plot a multi-concentration overlay, or to inspect a
suspicious binder.

## Provenance

Uploaded to ProteinBase by the hackathon organisers when the wet-lab
data closed out. Same numbers and curves as our internal LIMS, with
internal identifiers stripped.

If a file is missing or looks wrong, regenerate from the public
collection at
`https://proteinbase.bio/collections/adaptyv-x-muni-hackathon-ai-agents-vs-humans`
and open a PR.

## Git LFS

Every binary artifact in this directory — and the top-level
`data/designs.parquet` and `data/raw_lab/*.csv` — lives behind
[Git LFS](https://git-lfs.com). The repo's `.gitattributes` is the
source of truth for the tracked patterns.

A fresh clone needs LFS installed and the blobs pulled:

```bash
git clone <repo>
cd agents_vs_humans_paper
git lfs install
git lfs pull        # ~162 MB across ~620 files
```

Without LFS, every file in `data/proteinbase/`, every parquet, every
font, and the raw-lab CSVs are **pointer stubs of <200 bytes** that
look like this:

```
version https://git-lfs.github.com/spec/v1
oid sha256:abc123…
size 145678
```

Any code that tries to load one of these stubs as a CIF / JSON / PNG
will fail with a parse error. If you see "unexpected end of file" or
"not a PNG" errors after cloning, run `git lfs pull` and try again.
