#!/bin/bash
#
# make_rerun.sh - Build a disBatch job that reruns models whose last run ended with
#                 "termination code: min_timestep_limit".
#
# Each model resumes from its latest photo (or starts fresh if it has none),
# with whatever inlist_base is current (min_timestep_limit = 1d-8 since 2026-09-23).
#
# Selection:
#   - only models whose rn.out contains a termination code (a running model never has one,
#     since rn.out is rewritten when a run starts);
#   - by default only failures that ran with the OLD limit (MESA echoes the limit it used,
#     e.g. "min_timestep_limit 1.0000000000000001D-01"); those that already failed with
#     1d-8 need a physics/numerics change, not a rerun. --all includes them too;
#   - models listed in the task file of any rn_rerun job still queued or running are
#     skipped, so nothing is picked up twice.
#
# Every batch gets its own files (rerun_<stamp>_tasks.txt, submit_rerun_<stamp>.sh), so
# building a batch never touches the files of a job that is still running.
#
# Usage (from grids/):
#   ./make_rerun.sh [--all] [--max-nodes N]     # default N = 6 (8 models per node)
#   sbatch submit_rerun_<stamp>.sh              # command is printed at the end
#
set -euo pipefail

GRID_BASE="$(cd "$(dirname "$0")" && pwd)"
cd "$GRID_BASE"

OMP=16
TASK_TIMEOUT=48h   # guard against a model crawling at tiny dt for the whole allocation
MAX_NODES=6        # fewer slots than models: quick failures hand their slot to the next model
INCLUDE_NEW=0

while [ $# -gt 0 ]; do
    case "$1" in
        --all) INCLUDE_NEW=1 ;;
        --max-nodes) MAX_NODES="$2"; shift ;;
        *) echo "unknown option $1"; exit 1 ;;
    esac
    shift
done

# Models belonging to rerun jobs that are still queued/running: read each job's submit
# script, take its task file, collect the model directories in it.
ACTIVE=$(mktemp)
trap 'rm -f "$ACTIVE"' EXIT
for script in $(squeue -u "$USER" -h -n rn_rerun -o %o); do
    [ -f "$script" ] || continue
    tf=$(awk '/^disBatch/ {gsub(/"/, "", $NF); print $NF}' "$script")
    [ -f "$tf" ] && grep -o "${GRID_BASE}/[^ ]*/w[^ ]*/M[^ ]*" "$tf" | sed "s|${GRID_BASE}/||" >> "$ACTIVE"
done
n_active=$(sort -u "$ACTIVE" | wc -l)

STAMP=$(date '+%y%m%d_%H%M')
TASKFILE="$GRID_BASE/rerun_${STAMP}_tasks.txt"
SUBMIT="$GRID_BASE/submit_rerun_${STAMP}.sh"
> "$TASKFILE"

n=0; n_skip_new=0; n_skip_active=0
for dir in $(ls -d */w*/M* | sort); do
    out="$dir/rn.out"
    [ -f "$out" ] || continue
    tailtxt=$(tail -c 200000 "$out")
    code=$(grep -a "termination code" <<< "$tailtxt" | tail -1 | awk '{print $3}' || true)
    [ "$code" = "min_timestep_limit" ] || continue
    if grep -qx "$dir" "$ACTIVE"; then
        n_skip_active=$((n_skip_active + 1)); continue
    fi
    lim=$(grep -a -o "min_timestep_limit *[0-9.]*D[-+][0-9]*" <<< "$tailtxt" | tail -1 | awk '{print $2}' | tr D E)
    if [ "$INCLUDE_NEW" -eq 0 ] && [ -n "$lim" ] && awk -v l="$lim" 'BEGIN {exit !(l < 1e-3)}'; then
        n_skip_new=$((n_skip_new + 1)); continue
    fi
    echo "cd ${GRID_BASE}/${dir} && export OMP_NUM_THREADS=${OMP}"' && photo=$(ls -1t photos/ 2>/dev/null | head -1) && if [ -n "$photo" ]; then timeout '"${TASK_TIMEOUT}"' ./re $photo > rn.out 2>&1; else timeout '"${TASK_TIMEOUT}"' ./rn > rn.out 2>&1; fi' >> "$TASKFILE"
    echo "$(date '+%F %T') selected for rerun ($(basename "$TASKFILE")): $dir (old limit: ${lim:-?}, photo: $(ls -1t $dir/photos/ 2>/dev/null | head -1))" >> "$GRID_BASE/rerun_log.txt"
    n=$((n + 1))
done

echo "skipped: ${n_skip_active} in active rerun jobs (${n_active} models listed), ${n_skip_new} already failed with the new limit (use --all to include)"

if [ "$n" -eq 0 ]; then
    echo "No models to rerun."
    rm -f "$TASKFILE"
    exit 0
fi

# 128 cores/node / 16 threads = 8 models per node
NODES=$(( (n + 7) / 8 ))
[ "$NODES" -gt "$MAX_NODES" ] && NODES=$MAX_NODES

cat > "$SUBMIT" << SUBMIT_EOF
#!/bin/bash
#SBATCH --job-name=rn_rerun
#SBATCH --partition=cca
#SBATCH --account=cca
#SBATCH --qos=gen
#SBATCH --nodes=${NODES}
#SBATCH --exclusive
#SBATCH --time=7-00:00:00
#SBATCH --output=rerun_%j.log
#SBATCH --error=rerun_%j.err
#
# ${n} models that ended with min_timestep_limit (built $(date '+%F %T') by make_rerun.sh)

module load modules/2.4-20250724
module load disBatch

export MESASDK_ROOT=/mnt/home/mcantiello/mesasdk-26.6.1
source "\$MESASDK_ROOT/bin/mesasdk_init.sh"
export MESA_DIR=/mnt/home/mcantiello/mesa-26.04.1
export OMP_NUM_THREADS=${OMP}

cd "${GRID_BASE}"
disBatch -c ${OMP} "${TASKFILE}"
SUBMIT_EOF
chmod +x "$SUBMIT"

echo "${n} models -> $(basename "$TASKFILE") (${NODES} nodes, $((NODES * 8)) concurrent). Submit with:"
echo "  sbatch $SUBMIT"
