#!/bin/bash
#
# setup_grids.sh - Create the full multi-Z, multi-rotation MESA grid
#
# Creates 12 sub-grids: 3 metallicities (MW, LMC, SMC) x 4 rotation rates (0, 0.2, 0.4, 0.6)
# Each sub-grid has 34 mass models from 5.0 to 120 Msun.
#
# Usage:
#   ./setup_grids.sh              # Set up all 12 sub-grids
#   ./setup_grids.sh MW           # Set up all 4 rotation rates for MW
#   ./setup_grids.sh MW w0.2      # Set up only MW/w0.2
#
# Prerequisites:
#   - template_new/star must be compiled (cd ../template_new && ./mk)
#   - Run from the grids/ directory
#
set -euo pipefail

# MESA environment
export MESASDK_ROOT=/mnt/home/mcantiello/mesasdk-26.6.1
source "$MESASDK_ROOT/bin/mesasdk_init.sh"
export MESA_DIR=/mnt/home/mcantiello/mesa-26.04.1

# Directories
GRID_BASE="$(cd "$(dirname "$0")" && pwd)"
MODELS_DIR="$(dirname "$GRID_BASE")"
TEMPLATE="$MODELS_DIR/template_new"
CEPH_BASE="/mnt/ceph/users/mcantiello/rednoise/grids"

# Verify template
if [ ! -x "$TEMPLATE/star" ]; then
    echo "ERROR: $TEMPLATE/star not found or not executable."
    echo "Compile the template first: cd $TEMPLATE && ./mk"
    exit 1
fi

# Mass grid: 5 to 120 Msun
MASSES=(5.0 5.2 5.4 5.6 5.8 6.0 6.5 7.0 7.5 8.0 9.0 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 30 40 50 60 80 100 120)

# Metallicity labels
Z_LABELS=(MW LMC SMC MW_mltpp)

# Rotation labels and corresponding omega/omega_crit values
W_LABELS=(w0.0 w0.2 w0.4 w0.6)
W_VALUES=(0.0 0.2 0.4 0.6)

# Optional filtering
FILTER_Z="${1:-all}"
FILTER_W="${2:-all}"

OMP=16

# ============================================================
# Function: write the per-model inlist (router file)
# ============================================================
write_inlist() {
    local model_dir="$1"
    cat > "$model_dir/inlist" << 'INLIST_EOF'

&star_job

      read_extra_star_job_inlist(1) = .true.
      extra_star_job_inlist_name(1) = '../../../inlist_base'

      read_extra_star_job_inlist(2) = .true.
      extra_star_job_inlist_name(2) = '../../inlist_common'

      read_extra_star_job_inlist(3) = .true.
      extra_star_job_inlist_name(3) = 'inlist_grid'

/ ! end of star_job namelist

&eos

      read_extra_eos_inlist(1) = .true.
      extra_eos_inlist_name(1) = '../../../inlist_base'

/ ! end of eos namelist

&kap

      read_extra_kap_inlist(1) = .true.
      extra_kap_inlist_name(1) = '../../../inlist_base'

      read_extra_kap_inlist(2) = .true.
      extra_kap_inlist_name(2) = '../../inlist_common'

/ ! end of kap namelist

&controls

      read_extra_controls_inlist(1) = .true.
      extra_controls_inlist_name(1) = '../../../inlist_base'

      read_extra_controls_inlist(2) = .true.
      extra_controls_inlist_name(2) = '../../inlist_common'

      read_extra_controls_inlist(3) = .true.
      extra_controls_inlist_name(3) = 'inlist_grid'

/ ! end of controls namelist

&pgstar

      read_extra_pgstar_inlist(1) = .true.
      extra_pgstar_inlist_name(1) = '../../../inlist_pgstar'

/ ! end of pgstar namelist
INLIST_EOF
}

# ============================================================
# Function: write inlist_grid for a non-rotating model
# ============================================================
write_inlist_grid_norot() {
    local model_dir="$1"
    local mass="$2"
    cat > "$model_dir/inlist_grid" << GRID_EOF

&star_job

      create_pre_main_sequence_model = .true.

/ ! end of star_job namelist

&controls

      initial_mass     = ${mass}
      xa_central_lower_limit_species(1) = 'he4'
      xa_central_lower_limit(1) = 1d-4

/ ! end of controls namelist

&pgstar

/ ! end of pgstar namelist
GRID_EOF
}

# ============================================================
# Function: write inlist_grid for a rotating model
# ============================================================
write_inlist_grid_rot() {
    local model_dir="$1"
    local mass="$2"
    local omega="$3"
    cat > "$model_dir/inlist_grid" << GRID_EOF

&star_job

      create_pre_main_sequence_model = .true.

      ! Rotation setup
      change_rotation_flag = .true.
      new_rotation_flag = .true.
      set_near_zams_omega_div_omega_crit_steps = 10
      near_zams_relax_omega_div_omega_crit = .true.
      num_steps_to_relax_rotation = 50

/ ! end of star_job namelist

&controls

      initial_mass     = ${mass}
      xa_central_lower_limit_species(1) = 'he4'
      xa_central_lower_limit(1) = 1d-4

      ! Rotation rate
      new_omega_div_omega_crit = ${omega}

      ! Angular momentum transport: Spruit-Tayler + hydrodynamic instabilities
      am_nu_visc_factor = 0
      am_D_mix_factor = 0.0333333333d0
      D_DSI_factor = 1
      D_SH_factor = 1
      D_SSI_factor = 1
      D_ES_factor = 1
      D_GSF_factor = 1
      D_ST_factor = 1

/ ! end of controls namelist

&pgstar

/ ! end of pgstar namelist
GRID_EOF
}

# ============================================================
# Function: write submit_grid.sh for a sub-grid
# ============================================================
write_submit_script() {
    local subgrid_dir="$1"
    local z_label="$2"
    local w_label="$3"
    cat > "$subgrid_dir/submit_grid.sh" << SUBMIT_EOF
#!/bin/bash
#SBATCH --job-name=rn_${z_label}_${w_label}
#SBATCH --partition=cca
#SBATCH --account=cca
#SBATCH --qos=gen
#SBATCH --nodes=6
#SBATCH --exclusive
#SBATCH --time=7-00:00:00
#SBATCH --output=grid_%j.log
#SBATCH --error=grid_%j.err
#
# Sub-grid: ${z_label}/${w_label}
# 34 mass models, 16 OpenMP threads each
# 128 cores/node => 8 models per node => all 34 fit in 1 wave on 6 nodes

module load modules/2.4-20250724
module load disBatch

SUBGRID_DIR=${subgrid_dir}
OMP=${OMP}

export MESASDK_ROOT=/mnt/home/mcantiello/mesasdk-26.6.1
source "\$MESASDK_ROOT/bin/mesasdk_init.sh"
export MESA_DIR=/mnt/home/mcantiello/mesa-26.04.1
export OMP_NUM_THREADS=\$OMP

cd "\$SUBGRID_DIR"

# Generate disBatch task file
TASKFILE="\$SUBGRID_DIR/disBatch_tasks.txt"
> "\$TASKFILE"

for dir in \$(ls -d M* 2>/dev/null | sort -t'M' -k2 -g); do
    if [ -d "\$dir" ] && [ -e "\$dir/star" ]; then
        echo "cd \${SUBGRID_DIR}/\${dir} && export OMP_NUM_THREADS=\${OMP}"' && photo=\$(ls -1t photos/ 2>/dev/null | head -1) && if [ -n "\$photo" ]; then ./re \$photo > rn.out 2>&1; else ./rn > rn.out 2>&1; fi' >> "\$TASKFILE"
    fi
done

NTASKS=\$(wc -l < "\$TASKFILE")
echo "Sub-grid: ${z_label}/${w_label}"
echo "Tasks: \$NTASKS, Cores/task: \$OMP"

disBatch -c "\$OMP" "\$TASKFILE"
SUBMIT_EOF
    chmod +x "$subgrid_dir/submit_grid.sh"
}

# ============================================================
# Main loop
# ============================================================

echo "================================================================"
echo "MESA Grid Setup: 3 metallicities x 4 rotation rates = 12 sub-grids"
echo "Template: $TEMPLATE"
echo "Ceph output: $CEPH_BASE"
echo "================================================================"
echo ""

total_models=0

for zi in "${!Z_LABELS[@]}"; do
    z_label="${Z_LABELS[$zi]}"

    if [ "$FILTER_Z" != "all" ] && [ "$FILTER_Z" != "$z_label" ]; then continue; fi

    z_dir="$GRID_BASE/$z_label"

    # Check that inlist_common exists for this metallicity
    if [ ! -f "$z_dir/inlist_common" ]; then
        echo "ERROR: $z_dir/inlist_common not found."
        echo "Create it before running setup."
        exit 1
    fi

    for wi in "${!W_LABELS[@]}"; do
        w_label="${W_LABELS[$wi]}"
        omega="${W_VALUES[$wi]}"

        if [ "$FILTER_W" != "all" ] && [ "$FILTER_W" != "$w_label" ]; then continue; fi

        w_dir="$z_dir/$w_label"
        mkdir -p "$w_dir"

        echo "Setting up ${z_label}/${w_label} (omega/omega_crit = ${omega})"

        # Create model directories
        for mass in "${MASSES[@]}"; do
            dir="$w_dir/M${mass}"

            if [ -d "$dir" ]; then
                echo "  SKIP: M${mass} already exists"
                total_models=$((total_models + 1))
                continue
            fi

            mkdir -p "$dir"

            # Symlink binary, scripts, make, src to template
            for f in star rn re mk clean README_first; do
                ln -s "$TEMPLATE/$f" "$dir/$f"
            done
            ln -s "$TEMPLATE/make" "$dir/make"
            ln -s "$TEMPLATE/src" "$dir/src"

            # Create LOGS and photos on Ceph, symlink back
            ceph_dir="$CEPH_BASE/$z_label/$w_label/M${mass}"
            mkdir -p "$ceph_dir/LOGS" "$ceph_dir/photos"
            ln -s "$ceph_dir/LOGS" "$dir/LOGS"
            ln -s "$ceph_dir/photos" "$dir/photos"

            # Write inlist (router file)
            write_inlist "$dir"

            # Write inlist_grid (mass + rotation + stopping condition)
            if [ "$omega" = "0.0" ]; then
                write_inlist_grid_norot "$dir" "$mass"
            else
                write_inlist_grid_rot "$dir" "$mass" "$omega"
            fi

            total_models=$((total_models + 1))
        done

        # Generate submit script for this sub-grid
        write_submit_script "$w_dir" "$z_label" "$w_label"

        echo "  Created ${#MASSES[@]} model directories + submit_grid.sh"
        echo ""
    done
done

echo "================================================================"
echo "Grid setup complete: $total_models model directories created."
echo ""
echo "To submit, run these commands (one per sub-grid):"
echo ""
for z_label in "${Z_LABELS[@]}"; do
    if [ "$FILTER_Z" != "all" ] && [ "$FILTER_Z" != "$z_label" ]; then continue; fi
    for w_label in "${W_LABELS[@]}"; do
        if [ "$FILTER_W" != "all" ] && [ "$FILTER_W" != "$w_label" ]; then continue; fi
        echo "  sbatch $GRID_BASE/$z_label/$w_label/submit_grid.sh"
    done
done
echo "================================================================"
