#!/bin/bash
# --chdir: logs go here, whatever directory sbatch is run from
#SBATCH --job-name=rn2_MW_nzfix
#SBATCH --chdir=/mnt/home/mcantiello/work/rednoise/models/grids_v2/MW
#SBATCH --partition=cca
#SBATCH --account=cca
#SBATCH --qos=gen
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=7-00:00:00
#SBATCH --output=rerun_nzfix_%j.log
#SBATCH --error=rerun_nzfix_%j.err
#
# v2 MW: M22 and M25 in all four rotation sub-grids crashed at models 1-4 (array bound in
# run_star_extras get_conv_regions_mlt, fixed 2026-09-24). Restart them from scratch with the
# rebuilt template_v2 binary. 8 models x 16 threads = 1 node.

module load modules/2.4-20250724
module load disBatch

export MESASDK_ROOT=/mnt/home/mcantiello/mesasdk-26.6.1
source "$MESASDK_ROOT/bin/mesasdk_init.sh"
export MESA_DIR=/mnt/home/mcantiello/mesa-26.04.1
export OMP_NUM_THREADS=16

disBatch -c 16 "/mnt/home/mcantiello/work/rednoise/models/grids_v2/MW/rerun_nzfix_tasks.txt"
