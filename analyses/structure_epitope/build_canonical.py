"""
Build the canonical design table for the Structure & Epitope Analysis section.

Reads the repo's source-of-truth (data/designs.parquet, n=141), adds explicit
population-membership flags and a mutually-exclusive tier label, and freezes a
schema version. EVERYTHING downstream in this section reads the output file and
nothing else.

Populations (Handoff Brief §3.3), computed here, not hard-coded:
  P_all       = all submitted                      (141)
  P_screened  = wet-lab tested (submitted_to_lab)  (100)
  P_expressed = expressed/loaded (expressed==True) ( 89)  <- denominator for binder-vs-nonbinder
  P_kd        = binders with a fitted, uncensored Kd (kd_arith_mean_nM_all notna) (36)

Run: uv run python analyses/structure_epitope/build_canonical.py
"""
import hashlib
import json
from pathlib import Path

import pandas as pd

SCHEMA_VERSION = "1.0.0"
REPO = Path(__file__).resolve().parents[2]
SRC = REPO / "data" / "designs.parquet"
OUTDIR = Path(__file__).resolve().parent / "data"
OUT = OUTDIR / "designs_canonical.parquet"

# TREM2 numbering: BLI construct chain A = Met + UniProt(Q9NZC2) 19-174 + GGGSGGGS + 10xHis,
# authored 1..175 in the shipped Boltz-2 CIFs. Verified: construct pos 27 = W = UniProt W44;
# construct pos 19 = C = UniProt C36.  => UniProt_resnum = construct_pos + 17.
NUMBERING = {
    "construct_len": 175,
    "chainA_is_target": True,
    "chainB_is_binder": True,
    "uniprot": "Q9NZC2",
    "uniprot_from_construct_pos": "uniprot = construct_pos + 17",
    "igsf_domain_uniprot": [19, 132],
    "igsf_domain_construct_pos": [2, 115],
    "stalk_uniprot": [133, 174],
    "stalk_construct_pos": [116, 157],
    "linker_construct_pos": [158, 165],
    "his_tag_construct_pos": [166, 175],
}


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> None:
    df = pd.read_parquet(SRC)
    assert len(df) == 141, f"expected 141 rows, got {len(df)}"

    scr = df["submitted_to_lab"].fillna(False).astype(bool)
    exp = (df["expressed"] == True).fillna(False)  # noqa: E712
    hit = (df["is_hit"] == True).fillna(False)  # noqa: E712
    has_kd = df["kd_arith_mean_nM_all"].notna()

    df["p_all"] = True
    df["p_screened"] = scr
    df["p_expressed"] = exp
    df["p_kd"] = has_kd

    # mutually-exclusive finest tier
    def tier(r):
        if not r["p_screened"]:
            return "non_screened"
        if not r["p_expressed"]:
            return "screened_not_expressed"
        if not (r["is_hit"] == True):  # noqa: E712
            return "expressed_nonbinder"
        if not r["p_kd"]:
            return "hit_unfittable"
        return "hit_fitted"

    df["population_tier"] = df.apply(tier, axis=1)

    # a light Kd-quality flag (NOT censoring — all fitted Kd are point estimates)
    df["kd_quality_flag"] = ""
    df.loc[df["assay_methods_mixed"] == True, "kd_quality_flag"] += "assay_mixed;"  # noqa: E712
    df.loc[df["kd_replicate_cv_pct"] > 50, "kd_quality_flag"] += "cv_gt50;"

    counts = df["population_tier"].value_counts().to_dict()
    pops = {
        "P_all": int(df["p_all"].sum()),
        "P_screened": int(df["p_screened"].sum()),
        "P_expressed": int(df["p_expressed"].sum()),
        "P_kd": int(df["p_kd"].sum()),
    }
    assert pops == {"P_all": 141, "P_screened": 100, "P_expressed": 89, "P_kd": 36}, pops
    assert sum(counts.values()) == 141

    OUTDIR.mkdir(parents=True, exist_ok=True)
    df.to_parquet(OUT, index=False)

    schema = {
        "schema_version": SCHEMA_VERSION,
        "n_rows": int(len(df)),
        "n_cols": int(df.shape[1]),
        "source_parquet": str(SRC.relative_to(REPO)),
        "source_sha256": sha256(SRC),
        "output_sha256": sha256(OUT),
        "populations": pops,
        "population_tier_counts": counts,
        "added_columns": [
            "p_all", "p_screened", "p_expressed", "p_kd",
            "population_tier", "kd_quality_flag",
        ],
        "trem2_numbering": NUMBERING,
        "notes": [
            "P_expressed (expressed==True) is the denominator for binder-vs-non-binder comparisons.",
            "P_kd = 36: 37 hits minus design 5 (NovoFy, binding_label=weak) which has no fittable Kd.",
            "No interval censoring exists; every fitted Kd is a point estimate from >=2 replicate fits.",
        ],
    }
    (OUTDIR / "canonical_schema.json").write_text(json.dumps(schema, indent=2))
    print(json.dumps({"populations": pops, "tiers": counts,
                      "out": str(OUT.relative_to(REPO)),
                      "schema_version": SCHEMA_VERSION}, indent=2))


if __name__ == "__main__":
    main()
