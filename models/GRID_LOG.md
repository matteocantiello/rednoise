# MESA grid log

Chronological record of the rednoise MESA grids: what was run, what changed, and why.
Newest entries at the bottom. Paths are relative to `models/` unless absolute.

MESA r26.04.1 (`/mnt/home/mcantiello/mesa-26.04.1`), SDK 26.6.1. Outputs (LOGS, photos) live on
Ceph under `/mnt/ceph/users/mcantiello/rednoise/` and are symlinked into each model directory.

---

## 2026-09-22 — v1 grid set up

- `grids/setup_grids.sh` builds `grids/{MW,LMC,SMC}/{w0.0,w0.2,w0.4,w0.6}/M*`: 34 masses
  (5–120 Msun) × 3 metallicities (Z = 0.014 / 0.006 / 0.002) × 4 rotation rates
  (Ω/Ω_crit = 0, 0.2, 0.4, 0.6). Template: `template_new`. Each sub-grid runs with disBatch
  through its own `submit_grid.sh`. Output goes to `/mnt/ceph/users/mcantiello/rednoise/grids`.
- Shared physics in `grids/inlist_base`: α_MLT = 1.6, exponential overshoot (f = 0.014, f0 = 0.005),
  Type2 opacities, gold tolerances, mesh_delta_coeff = 0.3. The metallicity is set in each
  `inlist_common`.

## 2026-09-23 — v1 restarted with mass loss; timestep floor lowered

- Added Dutch winds (Vink+01 hot / de Jager cool, scaling factor 1.0) to `inlist_base`, then
  restarted the whole grid from scratch. The old no-mass-loss output is archived in
  `/mnt/ceph/users/mcantiello/rednoise/grids_nomassloss_20260923` and `old_nomassloss_logs/`.
- Added the differential sub-grid `grids/MW_mltpp/w0.0`: identical to MW/w0.0 but with MLT++
  (`okay_to_reduce_gradT_excess = .true.`), to measure its bias on FeCZ properties. MLT++ is
  OFF everywhere else.
- `min_timestep_limit` lowered from 1d-1 to 1d-8: 0.1 s aborted converging models, e.g.
  MW M120 on the pre-MS.
- The rotating models' `inlist_grid` files were **edited by hand after setup**: rotation is off
  until near the ZAMS, `new_omega_div_omega_crit` moved to `&star_job`, and a "rotational mass
  loss feedback" block was added (max_mdot_redo_cnt = 5, surf_omega_div_omega_crit_limit = 0.99,
  mdot_revise_factor = 1.1, implicit_mdot_boost = 0.1). **`grids/setup_grids.sh` does not
  reproduce these files**: re-running it would write rotating inlists that fail to parse.
- Rerun job 7097002 for the first batch of failures (resumes from the newest photo).

## 2026-09-24 — v1 diagnostics and comparison with observations

- Diagnostic scripts in `grids/`: `grid_io.py` (status + incremental history cache at
  `/mnt/ceph/users/mcantiello/rednoise/cache/history_npz`), `plot_grid_status.py`, `plot_hrd.py`,
  `plot_sHRD.py`, `plot_mltpp_compare.py`.
- MLT++ bias (MW w0.0): none below ~19 Msun. From 20 Msun, v_conv,max is +0.03–0.06 dex,
  Mach +0.06–0.12, τ_conv −0.1 to −0.25 dex, and F_c,max up to +1 dex at 50–60 Msun.
  MLT++ M80–M120 are stuck on the pre-MS.
- `grids/make_rerun.sh` → `grids/submit_rerun_260924_1439.sh`: 121 models that failed under the
  old 0.1 s floor. Submitted as job **7102135**, which is pending on `QOSMaxCpuPerUserLimit`.
- Comparison with observations in `../analysis_mesa/` (see `REPORT_transfer.md`). The MESA
  section (§5.6) is drafted in `../paper_apj/`.

## 2026-09-24 — template_v2

New template `template_v2/`, kept separate from `template_new` so that running jobs and
jobs resuming from photos are unaffected. Changes to `src/run_star_extras.f90`:

- **IGW radiative damping**, from the FeCZ top and from the core boundary to the photosphere
  (quasi-adiabatic WKB, Zahn, Talon & Matias 1997): `igw_tau` gives τ_rad plus the evanescent
  tunnelling part τ_ev; `igw_nu_damp` gives the frequency where τ_rad = 1.
- **FeCZ finder** now uses `mlt_mixing_type` instead of `mixing_type`. MESA relabels
  cells as rotational mixing wherever D_rot > D_conv, which fragmented the FeCZ in rotating
  models (single-cell zones, 0.3–0.5 dex jitter). If a zone type is fragmented, the fragment
  with the largest F_conv is kept. Up to 8 regions (`max_cz`) are examined.
- **26 new history columns** (147 extra in total, 276 overall): `FeCZ_ncells`, `FeCZ_nfrag`,
  `HeII_nfrag`, `FeCZ_nHP_top`, `FeCZ_P_top`, `FeCZ_mass_above`, `FeCZ_omega_c`, `FeCZ_ell_eddy`,
  `FeCZ_tau_rad_*` (at ω_c, ω_c/3 and 3ω_c; ℓ = ℓ_eddy, 1, 5), `FeCZ_tau_ev_*`, `FeCZ_nu_damp_*`,
  `core_omega_c`, `core_ell_eddy`, `core_tau_rad_*`, `core_nu_damp_*`, `Hp_phot`.
- **Triggered profiles:** `x_ctrl(5)` = spacing in X_c on the MS, `x_ctrl(6)` = spacing in
  log Teff after the MS. Restart-safe (the state is kept in `s% xtra`, which photos save).
  They use high priority, so they survive the profile cap.
- **Optional Fe-bump mesh functions** (`use_other_mesh_functions`, `x_ctrl(1:4)`). They had no
  effect in the tests (see T1c below), so they stay off.
- **Bug fix, `conv_vel` → `mlt_vc`:** MESA overwrites `s% conv_vel` with 3·D_mix/(α H_P)
  (`star/private/mix_info.f90:362`), and D_mix includes rotational mixing. As a result, every
  velocity, Mach, turnover and ω_c column in **rotating** models was inflated. At 20 Msun with
  ω = 0.2, v_FeCZ was ×2.6 and the turnover time ×0.27 relative to ω = 0, even though the
  geometry and F_c were identical. **This affects every rotating v1 sub-grid and the rotating
  extracts in `../analysis_mesa/`.** Non-rotating results are unaffected. The pre-fix source is
  kept in `template_v2/src/run_star_extras.f90.pre_mltvc`.

### Tests (local NVMe, `/home/mcantiello/rednoise_tests`)

- **T1a/T1b/T1c** (20 Msun; ω = 0, 0.2, and 0.2 with mesh functions; stop at X_c = 0.3): the new
  FeCZ finder gives one zone of ~150–180 cells and triggered profiles work. The mesh functions
  change nothing (T1c = T1b).
- **Fix verified** with `check_mltvc.py` (T1a/T1b rerun on the rebuilt binary at 17:25): at
  X_c = 0.66 the rot/nonrot ratio is 1.02 for v_FeCZ,max, 1.03 for Mach, 0.99 for the turnover
  time and 0.86–0.88 for τ_rad. Before the fix these were 2.6, –, 0.27 and 0.005.
- **T2v0–v3** (5.6 Msun, ω = 0.6, restarted from the v1 photo at model 900; 2 h limit):
  unchanged v1 physics (v0) dies at the TAMS on `min_timestep_limit`. The rotation controls
  (v1–v3) get past the TAMS, and turning off the hard limits on lgTeff/lgL (v2, v3) makes them
  ~3× faster: 1263 and 1320 models vs 1043. None reached He burning within 2 h.

## 2026-09-24 — v2 MW pilot set up (not submitted)

Decisions: Asplund09 mixture; the rotation controls of T2v2; MW first, LMC/SMC once it is
validated.

- `grids_v2/` with `setup_grids_v2.sh`, adapted from v1. Template `template_v2`, Ceph tree
  `/mnt/ceph/users/mcantiello/rednoise/grids_v2`, job names `rn2_*`. The rotating `inlist_grid`
  mirrors the real (hand-edited) v1 files, not v1's setup script.
- `grids_v2/inlist_base` = v1 plus:
  - `initial_zfracs = 6` (AGSS09); `kap_file_prefix = 'a09'`, `kap_CO_prefix = 'a09_co'`,
    `kap_lowT_prefix = 'lowT_fa05_a09p'`. v1 used the MESA default GS98 for both, which is
    inconsistent with Z = 0.014.
  - `x_ctrl(5) = x_ctrl(6) = 0.05` (triggered profiles).
- Rotating models (`w0.2/0.4/0.6`): `surf_omega_div_omega_crit_limit = 0.98`, `tol = 0.02`,
  `max_mdot_redo_cnt = 200`, `mdot_revise_factor = 1.2`, `implicit_mdot_boost = 0.05`,
  `delta_lgTeff_hard_limit = delta_lgL_hard_limit = -1`.
- `./setup_grids_v2.sh MW` created 136 model directories. Smoke test (a copy of M20 ω = 0.2 in
  `/home/mcantiello/rednoise_tests/v2smoke`, 60 models): the inlists parse, the a09 tables load
  (their cache files are now in `$MESA_DIR/data/kap_data/cache`), and all 276 columns are written.
- To submit (not yet done; you are at the per-user CPU limit, behind rerun job 7102135):
  ```
  sbatch grids_v2/MW/w0.0/submit_grid.sh
  sbatch grids_v2/MW/w0.2/submit_grid.sh
  sbatch grids_v2/MW/w0.4/submit_grid.sh
  sbatch grids_v2/MW/w0.6/submit_grid.sh
  ```
- Later, LMC/SMC: copy `grids/{LMC,SMC}/inlist_common` into `grids_v2/` and run
  `./setup_grids_v2.sh LMC` (or `SMC`).

### Open items
- Whether the ω = 0.6 runs now reach He burning (the pilot will show).
- Redo the rotating extracts in `../analysis_mesa/` from v2. `REPORT_transfer.md` used the
  biased MW ω = 0.2/0.4 velocities.
- The abstract, conclusions, §5.2 and §5.5 of the paper have not been revised for the MESA
  results yet.

## 2026-09-24 — committed to git

Commit `d36582b` (pushed to origin/master): scripts, inlists, template sources, figures, analysis
summary tables, paper sources and handoff notes. The root `.gitignore` excludes MESA binaries and
build products, run output (LOGS/photos/rn.out, all on Ceph), Slurm/disBatch logs, the old
`models/M*/` run directories, per-model grid directories (regenerated by the setup scripts), the
large `analysis_mesa/data/mesa_{ms,post}*.csv` extracts (regenerate with `extract_mesa.py`), and
LaTeX build files. Nothing was deleted from disk.

## 2026-09-24 — repository reorganised

The Cantiello et al. 2021 material (paper/, figures/, notebooks, 15140/, magnetic_models/, notes/, and
the 2021 MESA setup that was in models/: template/, inlist_common, run_grid.sh) moved to `../cantiello2021/`.
The 2021 versions of the column lists and inlist_pgstar were restored there from commit 2251874. The
current work (models/, analysis_mesa/, paper_apj/, handoff/) did **not** move, so all live paths
(running jobs, submit scripts, symlinks, analysis scripts) are unchanged. The git tag `cantiello2021`
(= 2251874) preserves the old layout for external links. The layout is explained in `../README.md`.

## 2026-09-24 — first v2 submission failed (no models ran)

Jobs 7103427–7103430 (`rn2_MW_w0.{0,2,4,6}`), submitted 18:11, started 18:54, and all FAILED after 2 s
(exit 0:53). They were submitted from `/home/mcantiello/rednoise_tests`, which is local workstation disk
and does not exist on compute nodes, so Slurm could not use it as the working directory and no log was
written. Nothing reached `grids_v2` or Ceph. Fix: the v2 submit scripts (and `setup_grids_v2.sh`) now set
`#SBATCH --chdir=<sub-grid dir>`, so logs and tasks go to the sub-grid regardless of where `sbatch` runs.
Resubmitted: jobs 7104231 (w0.0), 7104232 (w0.2), 7104235 (w0.4), 7104236 (w0.6), pending on Priority with the
correct working directories (`grids_v2/MW/w0.*`).

## 2026-09-24 ~22:45 — v2 MW status and a crash fix

After 2 h 52 min (jobs 7104231/2/5/6): 100 of 136 done (He exhaustion), 28 running (23 still on the MS),
8 crashed. Status tooling now reads any grid: `RN_GRID=../grids_v2 python3 plot_grid_status.py` (from
`grids/`) writes `grids_v2/grid_status.{txt,png}`. The history cache for non-default grids is kept apart
(`.../cache/history_npz_grids_v2`), because cache files are keyed by sub-grid and mass only.

| sub-grid | done | running | crashed |
|---|---|---|---|
| MW ω=0 / 0.2 / 0.4 / 0.6 | 26 / 25 / 25 / 24 | 6 / 7 / 7 / 8 | 2 / 2 / 2 / 2 |

For comparison, v1 MW ω=0.6 had only 14 of 34 done after ~25 h, with 15 failures.

**Crash:** M22 and M25 in all four sub-grids, at models 1–4 (pre-MS): `Index '1466' of 's%r' above upper
bound 1465` in `get_conv_velocities` (run_star_extras.f90:1307), called from the extra history columns. Cause:
the new `get_conv_regions_mlt` scans all cells with T < 1e6 K. Early pre-MS models can be below that
everywhere and fully convective, so a region could end at the centre cell nz, and callers read r(k+1).
**Fix:** `n_limit = min(n_limit, s% nz - 1)` (backup `template_v2/src/run_star_extras.f90.pre_nzfix`).
Compile-tested in a scratch copy; M22 and M25 (ω=0) then ran cleanly past model 60 locally
(`/home/mcantiello/rednoise_tests/v2smoke/MW/w0.0/`). Other models are unaffected (the fix only changes
the pre-MS edge case).
**Rerun prepared, not submitted:** `grids_v2/MW/submit_rerun_nzfix.sh` (1 node, disBatch, 8 models from
scratch; task list `grids_v2/MW/rerun_nzfix_tasks.txt`). It requires the live `template_v2/star` to be rebuilt
first; running jobs keep their already-loaded binary.
**23:05:** live `template_v2/star` rebuilt with the fix (user approved); verified locally on v2 MW ω=0 M25 with
the live binary (20 models, no runtime error). Ready to submit: `sbatch grids_v2/MW/submit_rerun_nzfix.sh`.
