# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "modal>=1.0",
#     "pandas",
# ]
# ///
"""DE-STRESS structural-stability metrics on pre-computed TREM2 complex structures.

Reads complex structures persisted on the Modal Volume by the complex
predictors (``structures/{boltz2,protenix,chai,af2m}/<slug>.{cif,pdb}``) and
runs the Wells-Wood DE-STRESS panel:

  * PyRosetta REF15 energies (per-residue)
  * EvoEF2 stability
  * BuDEff BUDE force field (internal + interaction)
  * Biophysical properties (gravy, pI, instability index, MW)

Per-model outputs land on the volume:

* ``destress_{model}/{slug}.json`` — full DE-STRESS metric panel

Usage::

    cd <repo_root>
    # All four predictors
    modal run --detach scripts/modal/modal_destress_avh.py

    # Just one predictor
    modal run --detach scripts/modal/modal_destress_avh.py --predictors boltz2

    # Download finished results to disk
    modal run scripts/modal/modal_destress_avh.py --download
"""

from __future__ import annotations

import contextlib
import json
import os
from pathlib import Path

import modal

TIMEOUT_MIN = int(os.environ.get("TIMEOUT", 30))
CONCURRENCY = int(os.environ.get("CONCURRENCY", 10))

ALL_PREDICTORS = ("boltz2", "protenix", "chai", "af2m")


# ---------------------------------------------------------------------------
# Modal image — EvoEF2 + BuDEff + (optional) PyRosetta
# ---------------------------------------------------------------------------


def _build_evoef2() -> None:
    """Clone and compile EvoEF2 from source inside the image."""
    import subprocess

    subprocess.run(
        "git clone https://github.com/tommyhuangthu/EvoEF2.git /opt/evoef2"
        " && cd /opt/evoef2 && g++ -O3 -o EvoEF2 src/*.cpp -I src/",
        shell=True,
        check=True,
    )


destress_image = (
    modal.Image.micromamba(python_version="3.11")
    .apt_install("git", "wget", "gcc", "g++", "dssp")
    .micromamba_install(
        "pyrosetta",
        channels=["https://conda.graylab.jhu.edu", "conda-forge"],
    )
    .pip_install(
        "biopython>=1.80",
        "pandas",
    )
    .pip_install("Cython")
    .run_commands(
        "pip install ampal BUDEFF || echo 'BuDEff install failed, will skip BuDEff scoring'"
    )
    .run_function(_build_evoef2)
)

app = modal.App(os.environ.get("MODAL_APP_NAME", "avh-destress"))

RESULTS_VOLUME_NAME = os.environ.get("MODAL_RESULTS_VOLUME", "avh-rerun-results")
RESULTS_VOLUME = modal.Volume.from_name(RESULTS_VOLUME_NAME, create_if_missing=True)
RESULTS_DIR = f"/{RESULTS_VOLUME_NAME}"


# ---------------------------------------------------------------------------
# Scoring function
# ---------------------------------------------------------------------------


@app.function(
    image=destress_image,
    cpu=4,
    timeout=TIMEOUT_MIN * 60,
    max_containers=CONCURRENCY,
    volumes={RESULTS_DIR: RESULTS_VOLUME},
)
def score_destress(predictor: str, slug: str) -> dict:
    """Compute Rosetta, EvoEF2, BuDEff, and biophysical metrics for one structure."""
    import shutil
    import subprocess
    import tempfile

    RESULTS_VOLUME.reload()

    result_subdir = f"destress_{predictor}"
    result_path = Path(RESULTS_DIR) / result_subdir / f"{slug}.json"
    if result_path.exists():
        try:
            cached = json.loads(result_path.read_text())
            if cached.get("status") == "ok":
                return cached
        except Exception:
            pass

    empty = _empty_destress()

    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            struct_dir = Path(RESULTS_DIR) / "structures" / predictor
            pdb_path: Path | None = None
            for ext in (".pdb", ".cif"):
                src = struct_dir / f"{slug}{ext}"
                if not src.exists():
                    continue
                dest = Path(tmpdir) / f"{slug}{ext}"
                shutil.copy2(str(src), str(dest))
                if ext == ".cif":
                    from Bio.PDB import PDBIO, MMCIFParser

                    parser = MMCIFParser(QUIET=True)
                    structure = parser.get_structure("s", str(dest))
                    pdb_dest = Path(tmpdir) / f"{slug}.pdb"
                    io = PDBIO()
                    io.set_structure(structure)
                    io.save(str(pdb_dest))
                    pdb_path = pdb_dest
                else:
                    pdb_path = dest
                break

            if pdb_path is None or not pdb_path.exists():
                result = {
                    "slug": slug,
                    "predictor": predictor,
                    "status": "no_pdb",
                    **empty,
                }
                _save_result(result, predictor)
                return result

            metrics: dict = {}

            # --- Rosetta REF15 ---
            try:
                import pyrosetta as pr

                pr.init(
                    extra_options=(
                        "-ignore_unrecognized_res -ignore_zero_occupancy -mute all "
                        "-corrections::beta_nov16 true"
                    )
                )
                pose = pr.pose_from_pdb(str(pdb_path))
                sfxn = pr.get_fa_scorefxn()
                total_score = sfxn(pose)
                n_res = pose.total_residue()
                metrics["rosetta_total_per_aa"] = (
                    total_score / n_res if n_res > 0 else 0.0
                )

                energies = pose.energies()
                terms = {
                    "fa_atr": pr.rosetta.core.scoring.ScoreType.fa_atr,
                    "fa_rep": pr.rosetta.core.scoring.ScoreType.fa_rep,
                    "fa_sol": pr.rosetta.core.scoring.ScoreType.fa_sol,
                    "fa_elec": pr.rosetta.core.scoring.ScoreType.fa_elec,
                    "hbond_sr_bb": pr.rosetta.core.scoring.ScoreType.hbond_sr_bb,
                    "hbond_lr_bb": pr.rosetta.core.scoring.ScoreType.hbond_lr_bb,
                    "hbond_bb_sc": pr.rosetta.core.scoring.ScoreType.hbond_bb_sc,
                    "hbond_sc": pr.rosetta.core.scoring.ScoreType.hbond_sc,
                    "rama_prepro": pr.rosetta.core.scoring.ScoreType.rama_prepro,
                    "fa_dun": pr.rosetta.core.scoring.ScoreType.fa_dun,
                }
                for term_name, term_type in terms.items():
                    term_total = sum(
                        energies.residue_total_energies(i)[term_type]
                        for i in range(1, n_res + 1)
                    )
                    metrics[f"rosetta_{term_name}_per_aa"] = (
                        term_total / n_res if n_res > 0 else 0.0
                    )
            except Exception as e:
                print(f"  Rosetta failed for {predictor}/{slug}: {e}")

            # --- EvoEF2 ---
            try:
                evoef2_bin = "/opt/evoef2/EvoEF2"
                if Path(evoef2_bin).exists():
                    proc = subprocess.run(
                        f'{evoef2_bin} --command=ComputeStability --pdb="{pdb_path}"',
                        shell=True,
                        capture_output=True,
                        text=True,
                        timeout=120,
                        cwd=str(pdb_path.parent),
                    )
                    for line in proc.stdout.split("\n"):
                        if "Total" in line and "=" in line:
                            try:
                                val = float(line.split("=")[-1].strip())
                                from Bio.PDB import PDBParser

                                parser = PDBParser(QUIET=True)
                                structure = parser.get_structure("s", str(pdb_path))
                                n_res_bio = sum(1 for _ in structure.get_residues())
                                metrics["evoef2_total_per_aa"] = (
                                    val / n_res_bio if n_res_bio > 0 else 0.0
                                )
                            except Exception:
                                pass
                            break
            except Exception as e:
                print(f"  EvoEF2 failed for {predictor}/{slug}: {e}")

            # --- BuDEff ---
            try:
                import ampal
                import budeff

                assembly = ampal.load_pdb(str(pdb_path))
                if hasattr(assembly, "__iter__") and not isinstance(assembly, ampal.Assembly):
                    assembly = next(iter(assembly), None) if assembly else None
                if assembly is not None and isinstance(assembly, ampal.AmpalContainer):
                    assembly = assembly[0] if len(assembly) > 0 else None
                if assembly is not None and isinstance(assembly, ampal.Assembly):
                    buff_score = budeff.get_internal_energy(assembly)
                    n_res_buff = sum(1 for _ in assembly.get_monomers())
                    if n_res_buff > 0:
                        metrics["budeff_total_per_aa"] = (
                            buff_score.total_energy / n_res_buff
                        )
                        metrics["budeff_steric_per_aa"] = buff_score.steric / n_res_buff
                        metrics["budeff_desolvation_per_aa"] = (
                            buff_score.desolvation / n_res_buff
                        )
                        metrics["budeff_charge_per_aa"] = buff_score.charge / n_res_buff

                    chains = list(assembly)
                    if len(chains) >= 2:
                        inter_score = budeff.get_interaction_energy(chains)
                        metrics["budeff_interaction_total"] = inter_score.total_energy
                        metrics["budeff_interaction_steric"] = inter_score.steric
                        metrics["budeff_interaction_desolvation"] = inter_score.desolvation
                        metrics["budeff_interaction_charge"] = inter_score.charge
            except ImportError:
                print("  BuDEff/AMPAL not available, skipping")
            except Exception as e:
                print(f"  BuDEff failed for {predictor}/{slug}: {e}")

            # --- Biophysical properties (whole complex, A+B) ---
            try:
                from Bio.PDB import PDBParser
                from Bio.PDB.Polypeptide import protein_letters_3to1
                from Bio.SeqUtils.ProtParam import ProteinAnalysis

                parser = PDBParser(QUIET=True)
                structure = parser.get_structure("s", str(pdb_path))

                seq = ""
                for residue in structure.get_residues():
                    resname = residue.get_resname().strip().upper()
                    if resname in protein_letters_3to1:
                        seq += protein_letters_3to1[resname]

                if seq:
                    pa = ProteinAnalysis(seq)
                    metrics["isoelectric_point"] = pa.isoelectric_point()
                    metrics["molecular_weight"] = pa.molecular_weight()
                    metrics["num_residues"] = len(seq)
                    with contextlib.suppress(Exception):
                        metrics["instability_index"] = pa.instability_index()
                    with contextlib.suppress(Exception):
                        metrics["gravy"] = pa.gravy()
            except Exception as e:
                print(f"  Biophysical props failed for {predictor}/{slug}: {e}")

            result = {
                "slug": slug,
                "predictor": predictor,
                "status": "ok",
                **{**empty, **metrics},
            }
            _save_result(result, predictor)
            return result

    except Exception as e:
        import traceback

        print(f"DE-STRESS failed for {predictor}/{slug}: {traceback.format_exc()}")
        result = {
            "slug": slug,
            "predictor": predictor,
            "status": f"error: {e}",
            **empty,
        }
        _save_result(result, predictor)
        return result


def _empty_destress() -> dict:
    return {
        "rosetta_total_per_aa": None,
        "rosetta_fa_atr_per_aa": None,
        "rosetta_fa_rep_per_aa": None,
        "rosetta_fa_sol_per_aa": None,
        "rosetta_fa_elec_per_aa": None,
        "rosetta_hbond_sr_bb_per_aa": None,
        "rosetta_hbond_lr_bb_per_aa": None,
        "rosetta_hbond_bb_sc_per_aa": None,
        "rosetta_hbond_sc_per_aa": None,
        "rosetta_rama_prepro_per_aa": None,
        "rosetta_fa_dun_per_aa": None,
        "evoef2_total_per_aa": None,
        "budeff_total_per_aa": None,
        "budeff_steric_per_aa": None,
        "budeff_desolvation_per_aa": None,
        "budeff_charge_per_aa": None,
        "budeff_interaction_total": None,
        "budeff_interaction_steric": None,
        "budeff_interaction_desolvation": None,
        "budeff_interaction_charge": None,
        "isoelectric_point": None,
        "molecular_weight": None,
        "num_residues": None,
        "instability_index": None,
        "gravy": None,
    }


def _save_result(result: dict, predictor: str) -> None:
    out_dir = Path(RESULTS_DIR) / f"destress_{predictor}"
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / f"{result['slug']}.json").write_text(json.dumps(result))
    RESULTS_VOLUME.commit()


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

orchestrator_image = modal.Image.debian_slim(python_version="3.11").pip_install("pandas")


@app.function(
    image=orchestrator_image,
    timeout=24 * 3600,
    volumes={RESULTS_DIR: RESULTS_VOLUME},
)
def run_batch(predictors: list[str], slugs: list[str]) -> None:
    import pandas as pd

    RESULTS_VOLUME.reload()

    for predictor in predictors:
        struct_dir = Path(RESULTS_DIR) / "structures" / predictor
        if not struct_dir.exists():
            print(f"No structures dir for {predictor}, skipping.")
            continue

        available: list[str] = []
        for f in struct_dir.iterdir():
            if f.suffix in (".pdb", ".cif") and f.stem in slugs:
                available.append(f.stem)
        available = sorted(set(available))
        print(f"DE-STRESS/{predictor}: {len(available)} structures found")

        dest_dir = Path(RESULTS_DIR) / f"destress_{predictor}"
        completed: set[str] = set()
        if dest_dir.exists():
            for f in dest_dir.glob("*.json"):
                with contextlib.suppress(Exception):
                    data = json.loads(f.read_text())
                    if data.get("status") == "ok":
                        completed.add(data["slug"])

        pending = [s for s in available if s not in completed]
        print(f"  {len(completed)} done, {len(pending)} pending")
        if not pending:
            continue

        for done, result in enumerate(
            score_destress.map(
                [predictor] * len(pending),
                pending,
                return_exceptions=True,
                wrap_returned_exceptions=False,
            ),
            start=1,
        ):
            if isinstance(result, Exception):
                print(f"  exception: {result}")
            elif done % 25 == 0:
                print(f"  {done}/{len(pending)} done")

        RESULTS_VOLUME.reload()
        rows: list[dict] = []
        if dest_dir.exists():
            for f in dest_dir.glob("*.json"):
                with contextlib.suppress(Exception):
                    rows.append(json.loads(f.read_text()))
        if rows:
            df = pd.DataFrame(rows)
            csv_path = Path(RESULTS_DIR) / f"destress_{predictor}_summary.csv"
            df.to_csv(csv_path, index=False)
            RESULTS_VOLUME.commit()
            ok = int((df["status"] == "ok").sum())
            print(f"  destress_{predictor}: {ok}/{len(df)} ok.")


@app.local_entrypoint()
def main(
    designs_csv: str = "./data/designs.csv",
    predictors: str = "all",
    limit: int | None = None,
    download: bool = False,
    retry_failed: bool = False,
) -> None:
    """Trigger DE-STRESS scoring on persisted complexes, or pull results to disk."""
    import pandas as pd

    pred_list = (
        list(ALL_PREDICTORS) if predictors == "all" else [p.strip() for p in predictors.split(",")]
    )

    if download:
        _download_to_local(designs_csv, pred_list)
        return

    if retry_failed:
        for p in pred_list:
            _clear_failed(p)
        return

    df = pd.read_csv(designs_csv)
    df = df[df["sequence"].notna() & (df["sequence"].str.len() > 0)].copy()
    if limit:
        df = df.head(limit)

    slugs = [f"design_{int(d):03d}" for d in df["design_id"].tolist()]
    print(f"Triggering DE-STRESS batch: {len(slugs)} slugs across predictors={pred_list}")
    run_batch.remote(pred_list, slugs)
    print("Done. Pull results with `--download`.")


def _download_to_local(designs_csv: str, predictors: list[str]) -> None:
    """Pull all per-design DE-STRESS JSONs into ./data/metrics/destress_<model>/."""
    import pandas as pd

    df = pd.read_csv(designs_csv)
    expected = {f"design_{int(d):03d}" for d in df["design_id"].tolist()}

    total = 0
    for predictor in predictors:
        out_dir = Path("./data/metrics") / f"destress_{predictor}"
        out_dir.mkdir(parents=True, exist_ok=True)
        n = 0
        try:
            for entry in RESULTS_VOLUME.iterdir(f"destress_{predictor}"):
                if not entry.path.endswith(".json"):
                    continue
                payload = b"".join(RESULTS_VOLUME.read_file(entry.path))
                data = json.loads(payload)
                slug = data.get("slug") or Path(entry.path).stem
                if slug not in expected:
                    continue
                (out_dir / f"{slug}.json").write_bytes(payload)
                n += 1
        except Exception as e:
            print(f"  destress_{predictor}: no results yet ({e})")
        print(f"  destress_{predictor}: downloaded {n} JSONs.")
        total += n
    print(f"Downloaded {total} total JSONs.")


def _clear_failed(predictor: str) -> None:
    cleared = 0
    try:
        for entry in RESULTS_VOLUME.iterdir(f"destress_{predictor}"):
            if not entry.path.endswith(".json"):
                continue
            data = json.loads(b"".join(RESULTS_VOLUME.read_file(entry.path)))
            status = str(data.get("status", ""))
            if status.startswith("error") or status.startswith("failed"):
                RESULTS_VOLUME.remove_file(entry.path)
                cleared += 1
    except Exception:
        pass
    print(f"  destress_{predictor}: cleared {cleared} failed results")
