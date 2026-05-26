#!/usr/bin/env bash
# Poll Modal volume + local disk for new model outputs, regenerate
# data/grand_metrics.csv, commit + push to the current branch.
#
# Designed to run forever inside a Monitor — every iteration emits exactly
# one summary line to stdout (the event the Monitor surfaces).
#
# Usage:
#   bash scripts/data/poll_and_populate.sh         # one iteration
#   while true; do bash scripts/data/poll_and_populate.sh; sleep 300; done
#
# Env knobs:
#   AVH_MODELS                       space-separated list of models to poll
#                                    (default: "boltz2 protenix chai")
#   AVH_SKIP_MODAL=1                 don't even try Modal CLI; only watch local disk
#   AVH_SKIP_PUSH=1                  commit but don't push
set -u
shopt -s nullglob

cd "$(git rev-parse --show-toplevel)"

MODELS="${AVH_MODELS:-boltz2 protenix chai af2m}"
TS="$(date -u +%H:%M:%S)"

script_for() {
    case "$1" in
        boltz2) echo "scripts/modal/modal_boltz2_avh.py" ;;
        protenix) echo "scripts/modal/modal_protenix_avh.py" ;;
        chai) echo "scripts/modal/modal_chai1_avh.py" ;;
        af2m) echo "scripts/modal/modal_af2m_avh.py" ;;
        *) echo "" ;;
    esac
}

cif_count() {
    # Portable file-count under a directory. AF2-M writes .pdb; rest write .cif.
    local dir="$1"
    local n=0
    if [[ -d "$dir" ]]; then
        n=$(find "$dir" -maxdepth 1 \( -name "*.cif" -o -name "*.pdb" \) -type f 2>/dev/null \
            | wc -l | tr -d ' ')
    fi
    echo "$n"
}

# Snapshot file counts BEFORE the pull. One line per model: "model count".
BEFORE_FILE="${TMPDIR:-/tmp}/avh_before_$$"
trap "rm -f '$BEFORE_FILE'" EXIT
: > "$BEFORE_FILE"
for m in $MODELS; do
    echo "$m $(cif_count data/structures/${m})" >> "$BEFORE_FILE"
done
echo "proteintyper $(cif_count data/structures/proteintyper)" >> "$BEFORE_FILE"

count_before() {
    awk -v m="$1" '$1==m {print $2; exit}' "$BEFORE_FILE"
}

# 1. Try to pull from Modal. Modal CLI may be down — that's fine, just note it.
modal_status="skipped"
if [[ "${AVH_SKIP_MODAL:-0}" != "1" ]]; then
    modal_status="ok"
    for m in $MODELS; do
        script="$(script_for "$m")"
        [[ -n "$script" && -f "$script" ]] || continue
        if ! modal run "$script" --download >/dev/null 2>&1; then
            modal_status="modal_unreachable"
        fi
    done
fi

# 2. Count files after the pull, compute deltas.
new_files=0
new_per_model=""
for m in $MODELS proteintyper; do
    after=$(cif_count "data/structures/${m}")
    before=$(count_before "$m")
    delta=$(( after - before ))
    if (( delta > 0 )); then
        new_files=$(( new_files + delta ))
        new_per_model="${new_per_model} ${m}+${delta}"
    fi
done

# 3. If nothing new on disk, emit a quiet heartbeat and return.
if (( new_files == 0 )); then
    echo "[poll ${TS}] no new structures (modal=${modal_status})"
    exit 0
fi

# 4. Rebuild data/grand_metrics.csv.
if ! PYTHONPATH=. python3 scripts/data/build_grand_metrics.py >/dev/null 2>&1; then
    echo "[poll ${TS}] FAILED grand_metrics rebuild — modal=${modal_status} new=${new_per_model}"
    exit 1
fi

# 5. Stage everything new under data/structures + data/metrics + grand csv.
git add data/structures data/metrics data/grand_metrics.csv 2>/dev/null

if git diff --cached --quiet; then
    echo "[poll ${TS}] no staged changes (despite ${new_files} new files? — investigate)"
    exit 0
fi

branch="$(git rev-parse --abbrev-ref HEAD)"
commit_msg="rerun:${branch} +${new_files} structures (${new_per_model# })"
if ! git commit -m "$commit_msg" --no-gpg-sign >/dev/null 2>&1; then
    echo "[poll ${TS}] FAILED commit — ${commit_msg}"
    exit 1
fi
short=$(git rev-parse --short HEAD)

# 6. Push (skip on the explicit knob).
push_status="skipped"
if [[ "${AVH_SKIP_PUSH:-0}" != "1" ]]; then
    if git push origin "$branch" >/dev/null 2>&1; then
        push_status="ok"
    else
        push_status="push_failed"
    fi
fi

echo "[poll ${TS}] +${new_files} structures —${new_per_model} → ${short} push=${push_status}"
