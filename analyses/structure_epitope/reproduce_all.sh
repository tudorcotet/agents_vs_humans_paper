#!/usr/bin/env bash
# One-command reproduction of the Structure & Epitope section's Tier A (permissive-stack)
# pipeline, from the committed source of truth to every committed 2D figure + results CSV.
#
#   source of truth:  data/designs.{csv,parquet}  +  data/structures/{protenix,chai,esmfold}/
#   regenerated here:  designs_canonical.parquet, footprints_{long,consensus}.parquet (all
#                      gitignored intermediates) and every results/*.csv + fig08e/09/11/12 +
#                      figS5/S6/S7/S8.
#
# Deterministic (all seeds fixed); re-running reproduces the committed outputs modulo float repr.
#
# NOT covered here — external tools / Tier B / GPU (see results/reproducibility_audit.md):
#   figS4 + fold_tmalign_clusters  -> jobs/tm_cluster.py  (needs TMalign binary + tm_allvsall.py; Tier B)
#   b2_*/af2m_* recompute, regen   -> jobs/recompute_scalars.py, regen_boltz2_footprints.py (GPU regen structures)
#   rosetta descriptors            -> PyRosetta (Tier B, .rosetta_venv)
#   EvoEF2 alanine scan            -> EvoEF2 binary
#   reference crystal epitopes     -> refs/*.cif from RCSB (run reference_epitopes.py; non-fatal below)
set -uo pipefail
cd "$(dirname "$0")/../.."                       # repo root
J=analyses/structure_epitope/jobs
UV="uv run --quiet --with pydssp --with biopython --with scikit-learn --with gemmi --with scipy"
pass=0; fail=0; failed=()
run() {
  local name="$1" script="$2"
  printf '── %-22s' "$name"
  if $UV python "$script" >"/tmp/repro_$name.log" 2>&1; then echo "ok"; pass=$((pass+1))
  else echo "FAIL (/tmp/repro_$name.log)"; fail=$((fail+1)); failed+=("$name"); fi
}

echo "== Stage 0: canonical table =="
run build_canonical      analyses/structure_epitope/build_canonical.py
echo "== Stage 1: footprints (structures -> long form) =="
run extract_footprints   $J/extract_footprints.py
echo "== Stage 2: epitope core (consensus + per-residue) =="
run build_epitope        $J/build_epitope.py
echo "== Stage 3: epitope downstream + interface descriptors =="
run epitope_binning      $J/epitope_binning.py
run epitope_stalk_refsets $J/epitope_stalk_refsets.py
run interface_descriptors $J/interface_descriptors.py
run rules_compliance     $J/rules_compliance.py
run robustness_HPerera   $J/robustness_HPerera.py
echo "== Stage 4: stats + fold topology =="
run fold_topology        $J/fold_topology.py
run phase2_stats         $J/phase2_stats.py
run phase6_correlations  $J/phase6_correlations.py
echo "== Stage 5: motifs / failure taxonomy / helical bias =="
run motifs               $J/motifs.py
run type_failure         $J/type_failure.py
run helical_bias         $J/helical_bias.py
echo "== Stage 6: figures =="
run figures              $J/figures.py
run supp_figures         $J/supp_figures.py
run supp_figures2        $J/supp_figures2.py
run fig_auroc_predictor  $J/fig_auroc_by_predictor.py   # Fig S8 (reads results/metric_predictor_table.csv)
echo "== (optional) reference crystal epitopes — needs refs/*.cif from RCSB =="
run reference_epitopes   $J/reference_epitopes.py

echo ""
echo "==================== SUMMARY: ${pass} ok, ${fail} failed ===================="
[ "$fail" -gt 0 ] && printf '  failed: %s\n' "${failed[@]}"
echo "Compare regenerated vs committed:  git -C \"$(pwd)\" status --short analyses/structure_epitope/results analyses/structure_epitope/figures"
exit 0
