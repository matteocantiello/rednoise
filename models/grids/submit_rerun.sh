#!/bin/bash
#SBATCH --job-name=rn_rerun
#SBATCH --partition=cca
#SBATCH --account=cca
#SBATCH --qos=gen
#SBATCH --nodes=2
#SBATCH --exclusive
#SBATCH --time=7-00:00:00
#SBATCH --output=rerun_%j.log
#SBATCH --error=rerun_%j.err
#
# 14 models that ended with min_timestep_limit (built 2026-09-23 18:44:53 by make_rerun.sh)

module load modules/2.4-20250724
module load disBatch

export MESASDK_ROOT=/mnt/home/mcantiello/mesasdk-26.6.1
source "$MESASDK_ROOT/bin/mesasdk_init.sh"
export MESA_DIR=/mnt/home/mcantiello/mesa-26.04.1
export OMP_NUM_THREADS=16

cd "/mnt/home/mcantiello/work/rednoise/models/grids"
disBatch -c 16 "/mnt/home/mcantiello/work/rednoise/models/grids/rerun_tasks.txt"
