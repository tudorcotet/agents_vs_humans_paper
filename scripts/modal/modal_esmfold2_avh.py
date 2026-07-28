# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "modal>=1.0",
#     "pandas",
# ]
# ///
"""ESMFold2 complex prediction + ipSAE scoring for the TREM2 (agents-vs-humans)
competition designs.

Same Modal-app shell as the sibling Boltz-2 / Protenix / Chai-1 / AF2-M
scripts in this folder (inline Dunbrack ipSAE + pDockQ + LIS scorer, results
volume, pb_id slug), folding with **local open-weight ESMFold2** on a GPU —
``ESMFold2Model.from_pretrained("biohub/ESMFold2")`` (MIT-licensed weights on
HuggingFace, ~14 GB on GPU). This is the same engine class as the other
folders (local GPU inference, no hosted API / credits).

Target / binder layout matches the other complex folders:

* **Target**: TREM2 ectodomain construct (175 aa) — chain ``A``.
* **Binder**: each row in ``data/designs.csv`` — chain ``B``.

ESMFold2 co-folds both chains and returns the full PAE matrix, so the
*identical* ``compute_ipsae`` used for the other four folders applies unchanged.

Outputs (Modal Volume ``avh-rerun-results``):

* ``esmfold2/{pb_id}.json`` — ipSAE + native ESMFold2 iptm / ptm / plddt
* ``structures/esmfold2/{pb_id}.cif``
* ``raw_data/esmfold2/{pb_id}.npz``

Usage::

    cd <repo_root>
    modal run --detach scripts/modal/modal_esmfold2_avh.py
    modal run scripts/modal/modal_esmfold2_avh.py --download
    modal run scripts/modal/modal_esmfold2_avh.py --limit 1   # smoke test
"""

from __future__ import annotations

import contextlib
import json
import os
from pathlib import Path

import modal

GPU = os.environ.get("GPU", "A100")
TIMEOUT_MIN = int(os.environ.get("TIMEOUT", 120))
CONCURRENCY = int(os.environ.get("CONCURRENCY", 10))

# Local ESMFold2 fold settings (≈ the "fast" preset: few recycles/diffusion steps).
NUM_LOOPS = int(os.environ.get("ESMFOLD2_LOOPS", 3))
NUM_SAMPLING_STEPS = int(os.environ.get("ESMFOLD2_STEPS", 32))
SEED = int(os.environ.get("ESMFOLD2_SEED", 42))
HF_MODEL = os.environ.get("ESMFOLD2_HF_MODEL", "biohub/ESMFold2")

# Modal app / volume names — overridable for external users.
APP_NAME = os.environ.get("MODAL_APP_NAME", "avh-esmfold2")
RESULTS_VOLUME_NAME = os.environ.get("MODAL_RESULTS_VOLUME", "avh-rerun-results")

# TREM2 ectodomain construct (175 aa) — same target the sibling avh folders used.
TREM2_TARGET_SEQ = (
    "MHNTTVFQGVAGQSLQVSCPYDSMKHWGRRKAWCRQLGEKGPCQRVVSTHNLWLLSFLRRWNGSTAITDDTLGGTLTITLRNLQPHDAGLYQCQSLHGSEADTLRKVLVEVLADPLDHRDAGDLWFPGESESFEDAHVEHSISRSLLEGEIPFPPTSGGGSGGGSHHHHHHHHHH"
)

PAE_CUTOFF = 15.0
DIST_CUTOFF = 8.0

PREDICTOR = "esmfold2"

# ---------------------------------------------------------------------------
# IPSAE scoring (identical to the other Modal apps in this folder)
# ---------------------------------------------------------------------------


def _calc_d0(n_res: float, min_value: float = 1.0) -> float:
    n_res = max(27.0, float(n_res))
    return max(min_value, 1.24 * (n_res - 15.0) ** (1.0 / 3.0) - 1.8)


def _ptm_func(pae_values, d0: float):
    return 1.0 / (1.0 + (pae_values / d0) ** 2.0)


def compute_ipsae(
    pae_matrix,
    structure_path: str,
    target_len: int,
    binder_len: int,
    target_chain: str = "A",
    binder_chain: str = "B",
    pae_cutoff: float = PAE_CUTOFF,
    dist_cutoff: float = DIST_CUTOFF,
) -> dict[str, float]:
    import numpy as np

    try:
        import gemmi
    except ImportError:
        return _empty_metrics()

    pae_matrix = np.asarray(pae_matrix, dtype=np.float64)
    total_len = target_len + binder_len

    try:
        st = gemmi.read_structure(str(structure_path))
    except Exception as e:
        print(f"Failed to read structure: {e}")
        return _empty_metrics()

    model = st[0]
    chain_names = [c.name for c in model]
    if target_chain not in chain_names or binder_chain not in chain_names:
        print(f"Chains {target_chain}/{binder_chain} not in {chain_names}")
        return _empty_metrics()

    def _extract(chain_id: str) -> list[list[float]]:
        out: list[list[float]] = []
        for res in model[chain_id]:
            atom = res.find_atom("CB", "*") if res.name != "GLY" else res.find_atom("CA", "*")
            if not atom:
                atom = res.find_atom("CA", "*")
            if atom:
                out.append([atom.pos.x, atom.pos.y, atom.pos.z])
        return out

    target_coords = _extract(target_chain)
    binder_coords = _extract(binder_chain)
    struct_target_len = len(target_coords)
    struct_binder_len = len(binder_coords)
    struct_total = struct_target_len + struct_binder_len

    if pae_matrix.shape[0] != struct_total:
        if pae_matrix.shape == (total_len, total_len) and total_len != struct_total:
            idx = list(range(struct_target_len)) + list(
                range(target_len, target_len + struct_binder_len)
            )
            if max(idx) < pae_matrix.shape[0]:
                pae_matrix = pae_matrix[np.ix_(idx, idx)]
            else:
                return _empty_metrics()
        else:
            return _empty_metrics()

    target_len = struct_target_len
    binder_len = struct_binder_len
    total_len = struct_total

    cb_coords = np.array(target_coords + binder_coords)
    distances = np.sqrt(((cb_coords[:, None, :] - cb_coords[None, :, :]) ** 2).sum(axis=2))
    target_mask = np.arange(total_len) < target_len
    binder_mask = ~target_mask
    pae_bt = pae_matrix[binder_mask][:, target_mask]
    pae_tb = pae_matrix[target_mask][:, binder_mask]

    n0chn = total_len
    d0chn = _calc_d0(n0chn)
    valid_bt = pae_bt < pae_cutoff
    valid_tb = pae_tb < pae_cutoff

    ipsae_d0chn_bt = float(_ptm_func(pae_bt[valid_bt], d0chn).mean()) if valid_bt.any() else 0.0
    ipsae_d0chn_tb = float(_ptm_func(pae_tb[valid_tb], d0chn).mean()) if valid_tb.any() else 0.0
    n0dom_bt = int(valid_bt.sum())
    ipsae_d0dom_bt = (
        float(_ptm_func(pae_bt[valid_bt], _calc_d0(n0dom_bt)).mean()) if n0dom_bt > 0 else 0.0
    )
    n0dom_tb = int(valid_tb.sum())
    ipsae_d0dom_tb = (
        float(_ptm_func(pae_tb[valid_tb], _calc_d0(n0dom_tb)).mean()) if n0dom_tb > 0 else 0.0
    )

    def _d0res_scores(pae_block, n_rows: int) -> float:
        vals: list[float] = []
        for i in range(n_rows):
            row = pae_block[i]
            good = row[row < pae_cutoff]
            if good.size > 0:
                vals.append(float(_ptm_func(good, _calc_d0(good.size)).mean()))
        return max(vals) if vals else 0.0

    ipsae_d0res_bt = _d0res_scores(pae_bt, binder_len)
    ipsae_d0res_tb = _d0res_scores(pae_tb, target_len)

    iptm_d0chn_bt = float(_ptm_func(pae_bt, d0chn).mean())
    iptm_d0chn_tb = float(_ptm_func(pae_tb, d0chn).mean())
    iptm_af_bt = float(_ptm_func(pae_bt, 10.0).mean())
    iptm_af_tb = float(_ptm_func(pae_tb, 10.0).mean())

    dist_bt = distances[binder_mask][:, target_mask]
    interface_mask = dist_bt <= dist_cutoff
    n_interface = int(interface_mask.sum())

    pdockq = pdockq2 = 0.0
    if n_interface > 0:
        binder_iface = np.where(interface_mask.any(axis=1))[0] + target_len
        target_iface = np.where(interface_mask.any(axis=0))[0]
        iface_idx = np.concatenate([binder_iface, target_iface])
        plddt_list = []
        for chain_id in [target_chain, binder_chain]:
            for res in model[chain_id]:
                ca = res.find_atom("CA", "*")
                if ca:
                    plddt_list.append(ca.b_iso)
        plddt = np.array(plddt_list)
        if len(plddt) == total_len:
            mean_iface_plddt = float(plddt[iface_idx].mean())
            x = mean_iface_plddt * np.log10(n_interface)
            pdockq = 0.724 / (1 + np.exp(-0.052 * (x - 152.611))) + 0.018
            ptm_iface = float(_ptm_func(pae_bt[interface_mask], 10.0).mean())
            x2 = mean_iface_plddt * ptm_iface
            pdockq2 = 1.31 / (1 + np.exp(-0.075 * (x2 - 84.733))) + 0.005

    pae_lis = pae_bt[pae_bt < 12.0]
    lis = float(((12.0 - pae_lis) / 12.0).mean()) if pae_lis.size > 0 else 0.0

    return {
        "ipsae_d0res_min": min(ipsae_d0res_bt, ipsae_d0res_tb),
        "ipsae_d0res_max": max(ipsae_d0res_bt, ipsae_d0res_tb),
        "ipsae_d0chn_min": min(ipsae_d0chn_bt, ipsae_d0chn_tb),
        "ipsae_d0chn_max": max(ipsae_d0chn_bt, ipsae_d0chn_tb),
        "ipsae_d0dom_min": min(ipsae_d0dom_bt, ipsae_d0dom_tb),
        "ipsae_d0dom_max": max(ipsae_d0dom_bt, ipsae_d0dom_tb),
        "iptm_d0chn_min": min(iptm_d0chn_bt, iptm_d0chn_tb),
        "iptm_d0chn_max": max(iptm_d0chn_bt, iptm_d0chn_tb),
        "iptm_af_min": min(iptm_af_bt, iptm_af_tb),
        "iptm_af_max": max(iptm_af_bt, iptm_af_tb),
        "pdockq": pdockq,
        "pdockq2": pdockq2,
        "lis": lis,
        "n_interface": n_interface,
    }


def _empty_metrics() -> dict[str, float]:
    keys = [
        "iptm",
        "ptm",
        "mean_plddt",
        "ipsae_d0res_min",
        "ipsae_d0res_max",
        "ipsae_d0chn_min",
        "ipsae_d0chn_max",
        "ipsae_d0dom_min",
        "ipsae_d0dom_max",
        "iptm_d0chn_min",
        "iptm_d0chn_max",
        "iptm_af_min",
        "iptm_af_max",
        "pdockq",
        "pdockq2",
        "lis",
        "n_interface",
    ]
    return dict.fromkeys(keys, 0.0)


# ---------------------------------------------------------------------------
# ESMFold2 helpers (shared across the four competition apps)
# ---------------------------------------------------------------------------


def _to_np(x):
    import numpy as np

    if x is None:
        return None
    if hasattr(x, "detach"):
        x = x.detach().cpu().numpy()
    return np.asarray(x, dtype=np.float64)


def _normalize_plddt(plddt):
    import numpy as np

    arr = _to_np(plddt)
    if arr is None or arr.size == 0:
        return arr
    return arr * 100.0 if float(np.nanmax(arr)) <= 1.0 else arr


def _write_cif_with_plddt(cif_text: str, plddt_per_res, out_path: Path) -> None:
    import gemmi
    import numpy as np  # noqa: F401

    out_path.write_text(cif_text)
    plddt = _normalize_plddt(plddt_per_res)
    if plddt is None or plddt.size == 0:
        return
    try:
        st = gemmi.read_structure(str(out_path))
        model = st[0]
        residues = [res for chain in model for res in chain]
        if len(residues) == len(plddt):
            for res, val in zip(residues, plddt, strict=True):
                for atom in res:
                    atom.b_iso = float(val)
        else:
            bvals = [atom.b_iso for chain in model for res in chain for atom in res]
            if bvals and max(bvals) <= 1.0:
                for chain in model:
                    for res in chain:
                        for atom in res:
                            atom.b_iso = float(atom.b_iso) * 100.0
        st.setup_entities()
        st.make_mmcif_document().write_file(str(out_path))
    except Exception as e:
        print(f"  warn: could not normalize B-factors ({e}); keeping raw CIF")


# ---------------------------------------------------------------------------
# Modal image — local ESMFold2 (GPU). Weights cached in a HF-cache volume.
# ---------------------------------------------------------------------------

esmfold2_image = (
    modal.Image.debian_slim(python_version="3.12")
    .apt_install("git")
    .pip_install(
        "esm @ git+https://github.com/Biohub/esm.git@main",
        "gemmi",
        "numpy<2.0",
    )
    .env({"HF_HUB_ENABLE_HF_TRANSFER": "0"})
)

app = modal.App(APP_NAME)
RESULTS_VOLUME = modal.Volume.from_name(RESULTS_VOLUME_NAME, create_if_missing=True)
RESULTS_DIR = f"/{RESULTS_VOLUME_NAME}"
HF_CACHE = modal.Volume.from_name("esmfold2-hf-cache", create_if_missing=True)
HF_CACHE_DIR = "/root/.cache/huggingface"


@app.cls(
    image=esmfold2_image,
    gpu=GPU,
    timeout=TIMEOUT_MIN * 60,
    max_containers=CONCURRENCY,
    volumes={RESULTS_DIR: RESULTS_VOLUME, HF_CACHE_DIR: HF_CACHE},
)
class ESMFold2Runner:
    @modal.enter()
    def load(self):
        import torch
        from esm.models.esmfold2 import ESMFold2InputBuilder
        from transformers.models.esmfold2.modeling_esmfold2 import ESMFold2Model

        self.torch = torch
        self.model = ESMFold2Model.from_pretrained(HF_MODEL).cuda().eval()
        self.builder = ESMFold2InputBuilder()
        HF_CACHE.commit()

    @modal.method()
    def predict(self, pb_id: str, binder_seq: str, target_seq: str) -> dict:
        safe = _sanitize(pb_id)
        result_path = Path(RESULTS_DIR) / PREDICTOR / f"{safe}.json"
        if result_path.exists():
            try:
                cached = json.loads(result_path.read_text())
                if cached.get("status") == "ok":
                    print(f"  {pb_id}: cached, skipping")
                    return cached
            except Exception:
                pass

        try:
            import numpy as np
            from esm.models.esmfold2 import ProteinInput, StructurePredictionInput

            target_len = len(target_seq)
            binder_len = len(binder_seq)

            fold_input = StructurePredictionInput(
                sequences=[
                    ProteinInput(id="A", sequence=target_seq),
                    ProteinInput(id="B", sequence=binder_seq),
                ]
            )
            with self.torch.no_grad():
                result = self.builder.fold(
                    self.model,
                    fold_input,
                    num_loops=NUM_LOOPS,
                    num_sampling_steps=NUM_SAMPLING_STEPS,
                    seed=SEED,
                    complex_id=safe,
                )
            if isinstance(result, list):
                result = result[0]

            pae_matrix = _to_np(result.pae)
            if pae_matrix is None:
                res = {"pb_id": pb_id, "predictor": PREDICTOR, "status": "failed_no_pae",
                       **_empty_metrics()}
                _save_result(res)
                return res

            plddt_arr = _normalize_plddt(result.plddt)
            native_iptm = float(result.iptm) if result.iptm is not None else 0.0
            native_ptm = float(result.ptm) if result.ptm is not None else 0.0
            native_plddt = float(np.nanmean(plddt_arr)) if plddt_arr is not None else 0.0

            out_dir = Path("/tmp/esmfold2_out")
            out_dir.mkdir(parents=True, exist_ok=True)
            cif_path = out_dir / f"{safe}.cif"
            _write_cif_with_plddt(result.complex.to_mmcif(), result.plddt, cif_path)

            scores = compute_ipsae(pae_matrix, str(cif_path), target_len, binder_len, "A", "B")
            scores["iptm"] = native_iptm
            scores["ptm"] = native_ptm
            scores["mean_plddt"] = native_plddt
            res = {"pb_id": pb_id, "predictor": PREDICTOR, "status": "ok", **scores}
            _save_result(res)
            _save_raw(pb_id, pae_matrix, native_iptm, native_ptm, native_plddt,
                      target_len, binder_len)
            _save_structure(cif_path, pb_id)
            print(
                f"  {pb_id}: iptm={native_iptm:.3f} ptm={native_ptm:.3f} "
                f"plddt={native_plddt:.1f} ipsae_d0chn_max={scores['ipsae_d0chn_max']:.3f}"
            )
            return res

        except Exception as e:
            import traceback

            print(f"ESMFold2 failed for {pb_id}: {traceback.format_exc()}")
            res = {"pb_id": pb_id, "predictor": PREDICTOR, "status": f"error: {e}",
                   **_empty_metrics()}
            _save_result(res)
            return res


def _sanitize(name: str) -> str:
    return (
        name.replace("/", "_SLASH_")
        .replace("\\", "_BSLASH_")
        .replace("|", "_PIPE_")
        .replace(" ", "_")
        .replace(",", "_COMMA_")
    )


def _save_result(result: dict) -> None:
    out_dir = Path(RESULTS_DIR) / PREDICTOR
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / f"{_sanitize(result['pb_id'])}.json").write_text(json.dumps(result))
    RESULTS_VOLUME.commit()


def _save_structure(local_path, pb_id: str) -> None:
    import shutil

    struct_dir = Path(RESULTS_DIR) / "structures" / PREDICTOR
    struct_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(str(local_path), str(struct_dir / f"{_sanitize(pb_id)}{Path(local_path).suffix}"))
    RESULTS_VOLUME.commit()


def _save_raw(
    pb_id: str,
    pae_matrix,
    native_iptm: float,
    native_ptm: float,
    native_plddt: float,
    target_len: int,
    binder_len: int,
) -> None:
    import numpy as np

    raw_dir = Path(RESULTS_DIR) / "raw_data" / PREDICTOR
    raw_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        str(raw_dir / f"{_sanitize(pb_id)}.npz"),
        pae=np.asarray(pae_matrix, dtype=np.float32),
        target_len=np.int32(target_len),
        binder_len=np.int32(binder_len),
        native_iptm=np.float32(native_iptm),
        native_ptm=np.float32(native_ptm),
        native_plddt=np.float32(native_plddt),
    )
    RESULTS_VOLUME.commit()


orchestrator_image = modal.Image.debian_slim(python_version="3.12").pip_install("pandas")


@app.function(image=orchestrator_image, timeout=24 * 3600, volumes={RESULTS_DIR: RESULTS_VOLUME})
def run_batch(pb_ids: list[str], seqs: list[str], target_seq: str) -> None:
    import pandas as pd

    RESULTS_VOLUME.reload()
    pred_dir = Path(RESULTS_DIR) / PREDICTOR
    completed: set[str] = set()
    if pred_dir.exists():
        for f in pred_dir.glob("*.json"):
            with contextlib.suppress(Exception):
                d = json.loads(f.read_text())
                if d.get("status") == "ok":
                    completed.add(d["pb_id"])

    pending = [(p, s) for p, s in zip(pb_ids, seqs, strict=True) if p not in completed]
    print(f"ESMFold2: {len(completed)} done, {len(pending)} pending")

    if pending:
        runner = ESMFold2Runner()
        for done, result in enumerate(
            runner.predict.map(
                [p[0] for p in pending],
                [p[1] for p in pending],
                [target_seq] * len(pending),
                return_exceptions=True,
            ),
            start=1,
        ):
            if isinstance(result, Exception):
                print(f"  exception: {result}")
            elif done % 25 == 0:
                print(f"  {done}/{len(pending)} done")

    RESULTS_VOLUME.reload()
    rows = []
    if pred_dir.exists():
        for f in pred_dir.glob("*.json"):
            with contextlib.suppress(Exception):
                rows.append(json.loads(f.read_text()))
    df = pd.DataFrame(rows)
    df.to_csv(Path(RESULTS_DIR) / f"{PREDICTOR}_summary.csv", index=False)
    RESULTS_VOLUME.commit()
    ok = int((df["status"] == "ok").sum()) if "status" in df.columns else 0
    print(f"ESMFold2: {ok}/{len(df)} succeeded.")


@app.local_entrypoint()
def main(
    designs_csv: str = "./data/designs.csv",
    target_fasta: str = "./data/target/trem2_construct.fasta",
    limit: int | None = None,
    download: bool = False,
    retry_failed: bool = False,
) -> None:
    import pandas as pd

    if download:
        _download_to_local(designs_csv)
        return

    if retry_failed:
        _clear_failed()
        return

    target_seq = _read_fasta(Path(target_fasta))
    if target_seq != TREM2_TARGET_SEQ:
        print(
            f"  note: {target_fasta} ({len(target_seq)} aa) differs from the "
            f"embedded TREM2_TARGET_SEQ ({len(TREM2_TARGET_SEQ)} aa); using the file."
        )

    df = pd.read_csv(designs_csv)
    df = df[df["sequence"].notna() & (df["sequence"].str.len() > 0)].copy()
    if limit:
        df = df.head(limit)
    # Key by the `design_NNN` slug the rest of this repo uses (build_grand_metrics
    # reads data/metrics/<model>/design_NNN.json), NOT pb_id.
    pb_ids = [f"design_{int(x):03d}" for x in df["design_id"]]
    seqs = df["sequence"].tolist()
    print(f"Triggering ESMFold2 batch: {len(pb_ids)} designs (target=TREM2, 175 aa)")
    call = run_batch.spawn(pb_ids, seqs, target_seq)
    print(f"Spawned run_batch ({call.object_id}); waiting for completion...")
    call.get()
    print("Done. Pull results with `--download`.")


def _read_fasta(path: Path) -> str:
    seq_parts: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">") or not line.strip():
            if seq_parts:
                break
            continue
        seq_parts.append(line.strip())
    return "".join(seq_parts).upper()


def _download_to_local(designs_csv: str) -> None:
    import pandas as pd

    df = pd.read_csv(designs_csv)
    expected_ids = {f"design_{int(x):03d}" for x in df["design_id"]}

    json_out = Path("./data/metrics") / PREDICTOR
    cif_out = Path("./data/structures") / PREDICTOR
    json_out.mkdir(parents=True, exist_ok=True)
    cif_out.mkdir(parents=True, exist_ok=True)

    n_json = n_cif = 0
    try:
        for entry in RESULTS_VOLUME.iterdir(PREDICTOR):
            if entry.path.endswith(".json"):
                payload = b"".join(RESULTS_VOLUME.read_file(entry.path))
                data = json.loads(payload)
                pb_id = data.get("pb_id") or Path(entry.path).stem
                if pb_id not in expected_ids:
                    continue
                (json_out / f"{pb_id}.json").write_bytes(payload)
                n_json += 1
    except Exception as e:
        print(f"  no metric JSONs: {e}")

    try:
        for entry in RESULTS_VOLUME.iterdir(f"structures/{PREDICTOR}"):
            if entry.path.endswith(".cif"):
                pb_id = Path(entry.path).stem
                if pb_id not in expected_ids:
                    continue
                payload = b"".join(RESULTS_VOLUME.read_file(entry.path))
                (cif_out / f"{pb_id}.cif").write_bytes(payload)
                n_cif += 1
    except Exception as e:
        print(f"  no CIFs: {e}")

    print(f"Downloaded {n_json} JSON + {n_cif} CIF for {PREDICTOR}.")


def _clear_failed() -> None:
    cleared = 0
    try:
        for entry in RESULTS_VOLUME.iterdir(PREDICTOR):
            if not entry.path.endswith(".json"):
                continue
            data = json.loads(b"".join(RESULTS_VOLUME.read_file(entry.path)))
            status = str(data.get("status", ""))
            if status.startswith("error") or status.startswith("failed"):
                RESULTS_VOLUME.remove_file(entry.path)
                cleared += 1
    except Exception:
        pass
    print(f"  cleared {cleared} failed results")
