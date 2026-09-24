#!/bin/bash
#SBATCH --job-name=rednoise_grid
#SBATCH --partition=cca
#SBATCH --account=cca
#SBATCH --qos=gen
#SBATCH --nodes=6
#SBATCH --exclusive
#SBATCH --time=7-00:00:00
#SBATCH --output=grid_%j.log
#SBATCH --error=grid_%j.err
#
# submit_grid.sh - Run the rednoise MESA grid using disBatch
#
# Each MESA model uses OMP_NUM_THREADS cores via OpenMP.
# disBatch distributes models across the allocated nodes.
#
# Usage:
#   1. Run setup_grid.sh first to create model directories
#   2. sbatch submit_grid.sh
#
# Resource math:
#   128 cores/node, OMP_NUM_THREADS=16 => 8 models per node
#   6 nodes => 48 models running concurrently
#   34 total models => completes in 1 wave

module load modules/2.4-20250724
module load disBatch

GRID_DIR=/mnt/home/mcantiello/work/rednoise/models
OMP=16

# MESA environment setup (will be inherited by tasks)
export MESASDK_ROOT=/mnt/home/mcantiello/mesasdk-26.6.1
source "$MESASDK_ROOT/bin/mesasdk_init.sh"
export MESA_DIR=/mnt/home/mcantiello/mesa-26.04.1
export OMP_NUM_THREADS=$OMP

cd "$GRID_DIR"

# Generate the disBatch task file
TASKFILE="$GRID_DIR/disBatch_tasks.txt"
> "$TASKFILE"

for dir in $(ls -d M* 2>/dev/null | sort -t'M' -k2 -g); do
    if [ -d "$dir" ] && [ -x "$dir/star" ]; then
        # Find the latest photo and restart from it; fall back to fresh run if no photos exist
        echo "cd ${GRID_DIR}/${dir} && export OMP_NUM_THREADS=${OMP} && photo=\$(ls -1 photos/ 2>/dev/null | sort -t'x' -k2 -g | tail -1) && if [ -n \"\$photo\" ]; then ./re \$photo > rn.out 2>&1; else ./rn > rn.out 2>&1; fi" >> "$TASKFILE"
    fi
done

NTASKS=$(wc -l < "$TASKFILE")
echo "Task file: $TASKFILE"
echo "Number of tasks: $NTASKS"
echo "Cores per task: $OMP"

# Run disBatch: -c sets cores per task
disBatch -c "$OMP" "$TASKFILE"
