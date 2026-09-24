# Rednoise

Stochastic low-frequency ("red noise") photometric variability and surface turbulence
(macroturbulence) in massive stars, and their link to subsurface convection zones.

This repository holds two projects:

| | Project | Where |
|---|---|---|
| 1 | **Cantiello et al. 2021**, *On the origin of stochastic, low-frequency photometric variability in massive stars*, ApJ 915, 112. Published; archived. | [`cantiello2021/`](cantiello2021/) |
| 2 | **2026 follow-up**, *Surface turbulence and stochastic low-frequency variability in massive stars share a sub-surface driver*. Draft manuscript, new MESA grids, comparison with observations. Ongoing. | top level: [`paper_apj/`](paper_apj/), [`analysis_mesa/`](analysis_mesa/), [`models/`](models/), [`handoff/`](handoff/) |

Shared reading: [`literature/`](literature/) and [`references/`](references/) (paper PDFs).

## 1. Cantiello et al. 2021 — `cantiello2021/`

The FeCZ properties of MESA models (r15140, non-rotating, Z = 0.02 / LMC / SMC) are compared with the
red-noise sample of Bowman et al. (2019, 2020). ν_char and α₀ correlate with FeCZ properties; core IGWs
explain the amplitudes but not the frequencies. Notebooks, figures, the manuscript and MESA setups are in
`cantiello2021/`; the model output is on Zenodo ([10.5281/zenodo.4676427](https://doi.org/10.5281/zenodo.4676427)).
See [`cantiello2021/README.md`](cantiello2021/README.md).

**Links to the old layout.** Until September 2026 these files sat at the repository root. The git tag
[`cantiello2021`](https://github.com/matteocantiello/rednoise/tree/cantiello2021) freezes that layout, so any
link `github.com/matteocantiello/rednoise/blob/master/<path>` still works as `.../blob/cantiello2021/<path>`.

## 2. 2026 follow-up — top level

| Path | What |
|---|---|
| [`paper_apj/`](paper_apj/) | AASTeX 6.3.1 manuscript (`main.tex`), tables (incl. machine-readable tables), figures. Build: `pdflatex main && bibtex main && pdflatex main && pdflatex main`. §5.6 and the appendix describe the MESA grid. |
| [`handoff/`](handoff/) | Observational catalogues and results the manuscript is built on: `project_handoff/PROJECT_HANDOFF.md` (overview), `project_handoff/data/` (red-noise and macroturbulence catalogues, fit tables, Bowman+2020 tables, `make_figures.py`). |
| [`models/`](models/) | New MESA r26.04.1 grids. **Start with [`models/GRID_LOG.md`](models/GRID_LOG.md)**, the dated record of every grid, physics change, bug and test. |
| [`analysis_mesa/`](analysis_mesa/) | Grid vs. observations: `extract_mesa.py` (grid → tables), `compare_obs.py`, transfer models (`transfer_models.py`, `test_transfer.py`, `scenarios.py`), figure scripts; findings in `REPORT_transfer.md`. |

### `models/` layout
| Path | What |
|---|---|
| `GRID_LOG.md` | History and status of the grids (read first) |
| `grids/` | **v1 grid**: {MW, LMC, SMC} × Ω/Ω_crit ∈ {0, 0.2, 0.4, 0.6} × 34 masses (5–120 Msun), plus the MLT++ test `MW_mltpp/w0.0`. Setup and rerun scripts, `inlist_base`, per-metallicity `inlist_common`, and diagnostics (`grid_io.py`, `plot_*.py`, `figures/`). Uses `template_new`. |
| `grids_v2/` | **v2 grid** (MW pilot): Asplund09 mixture and opacities, rotation fixes, `template_v2`. `setup_grids_v2.sh`. |
| `template_new/`, `template_v2/` | MESA work directories (`src/run_star_extras.f90`). v2 adds IGW damping diagnostics and a robust FeCZ finder, and fixes rotating-model convective velocities. |
| `history_columns.list`, `profile_columns.list`, `inlist_pgstar` | Shared by both grids |
| `setup_grid.sh`, `submit_grid.sh`, `inlist_common_new` | First (Sep 2026) single-directory grid attempt, superseded by `grids/` |

Model output (LOGS, photos; tens of GB) is not in git. It lives on Ceph at `/mnt/ceph/users/mcantiello/rednoise/`
and each model directory symlinks to it. The large extracted tables in `analysis_mesa/data/` are regenerated
with `python3 analysis_mesa/extract_mesa.py`.

**Paths.** The grid directories, Slurm scripts and analysis scripts use absolute paths under
`/mnt/home/mcantiello/work/rednoise/` (Flatiron cluster). Moving `models/`, `analysis_mesa/` or `handoff/`
would break running jobs and scripts; see `GRID_LOG.md` before reorganising them.
