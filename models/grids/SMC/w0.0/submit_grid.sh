#!/bin/bash
#SBATCH --job-name=rn_SMC_w0.0
#SBATCH --partition=cca
#SBATCH --account=cca
#SBATCH --qos=gen
#SBATCH --nodes=6
#SBATCH --exclusive
#SBATCH --time=7-00:00:00
#SBATCH --output=grid_%j.log
#SBATCH --error=grid_%j.err
#
# Sub-grid: SMC/w0.0
# 34 mass models, 16 OpenMP threads each
# 128 cores/node => 8 models per node => all 34 fit in 1 wave on 6 nodes

module load modules/2.4-20250724
module load disBatch

SUBGRID_DIR=/mnt/home/mcantiello/work/rednoise/models/grids/SMC/w0.0
OMP=16

export MESASDK_ROOT=/mnt/home/mcantiello/mesasdk-26.6.1
source "$MESASDK_ROOT/bin/mesasdk_init.sh"
export MESA_DIR=/mnt/home/mcantiello/mesa-26.04.1
export OMP_NUM_THREADS=$OMP

cd "$SUBGRID_DIR"

# Generate disBatch task file
TASKFILE="$SUBGRID_DIR/disBatch_tasks.txt"
> "$TASKFILE"

for dir in $(ls -d M* 2>/dev/null | sort -t'M' -k2 -g); do
    if [ -d "$dir" ] && [ -e "$dir/star" ]; then
        echo "cd ${SUBGRID_DIR}/${dir} && export OMP_NUM_THREADS=${OMP}"' && photo=$(ls -1t photos/ 2>/dev/null | head -1) && if [ -n "$photo" ]; then ./re $photo > rn.out 2>&1; else ./rn > rn.out 2>&1; fi' >> "$TASKFILE"
    fi
done

NTASKS=$(wc -l < "$TASKFILE")
echo "Sub-grid: SMC/w0.0"
echo "Tasks: $NTASKS, Cores/task: $OMP"

disBatch -c "$OMP" "$TASKFILE"
