#!/bin/bash
#
# setup_grid.sh - Create model directories from template_new for the rednoise grid
#
# Usage: ./setup_grid.sh
#
# This script creates a directory for each mass value, copies the template,
# compiles the code, and injects the initial mass into inlist_grid.
# Run this on the gate node BEFORE submitting the Slurm job.

set -euo pipefail

# MESA environment
export MESASDK_ROOT=/mnt/home/mcantiello/mesasdk-26.6.1
source "$MESASDK_ROOT/bin/mesasdk_init.sh"
export MESA_DIR=/mnt/home/mcantiello/mesa-26.04.1

GRID_DIR="$(cd "$(dirname "$0")" && pwd)"
TEMPLATE="template_new"

# Bulk output goes to ceph (large, no backup, appropriate for simulation data)
# Code and config stay on /mnt/home (small, backed up)
CEPH_OUTPUT="/mnt/ceph/users/mcantiello/rednoise/grid"

# Mass grid: 5 to 120 Msun
MASSES=(5.0 5.2 5.4 5.6 5.8 6.0 6.5 7.0 7.5 8.0 9.0 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 30 40 50 60 80 100 120)

cd "$GRID_DIR"

# Check template exists and is compiled
if [ ! -x "$TEMPLATE/star" ]; then
    echo "ERROR: $TEMPLATE/star not found or not executable."
    echo "Compile the template first: cd $TEMPLATE && ./mk"
    exit 1
fi

# Create ceph output directory
mkdir -p "$CEPH_OUTPUT"

echo "Setting up grid with ${#MASSES[@]} models in $GRID_DIR"
echo "Template: $TEMPLATE"
echo "Output:   $CEPH_OUTPUT"
echo "Masses: ${MASSES[*]}"
echo ""

for mass in "${MASSES[@]}"; do
    dir="M${mass}"
    if [ -d "$dir" ]; then
        echo "SKIP: $dir already exists"
        continue
    fi

    echo "Creating $dir (M = ${mass} Msun)..."

    # Create directory and copy template contents
    mkdir -p "$dir"
    cp "$TEMPLATE"/{inlist,inlist_grid,rn,re,mk,clean,star,README_first} "$dir/" 2>/dev/null || true
    cp -r "$TEMPLATE"/{make,src} "$dir/"

    # Create LOGS and photos on ceph, symlink from model directory
    mkdir -p "$CEPH_OUTPUT/$dir/LOGS" "$CEPH_OUTPUT/$dir/photos"
    ln -s "$CEPH_OUTPUT/$dir/LOGS" "$dir/LOGS"
    ln -s "$CEPH_OUTPUT/$dir/photos" "$dir/photos"

    # Inject the mass value
    sed -i "s/initial_mass     =/initial_mass     = ${mass} /" "$dir/inlist_grid"

    echo "  Done: $dir"
done

echo ""
echo "Grid setup complete. ${#MASSES[@]} model directories created."
echo "Output data (LOGS, photos) will be written to: $CEPH_OUTPUT"
echo "Next step: submit the Slurm job with sbatch submit_grid.sh"
