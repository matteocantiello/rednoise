# Cantiello et al. 2021: On the origin of stochastic, low-frequency photometric variability in massive stars

M. Cantiello, D. Lecoanet, A. S. Jermyn & L. Grassitelli, 2021, ApJ, 915, 112
([ADS](https://ui.adsabs.harvard.edu/abs/2021ApJ...915..112C)).
The inlists, processing scripts and model output are archived on Zenodo:
[10.5281/zenodo.4676427](https://doi.org/10.5281/zenodo.4676427).

High-precision photometry has revealed a ubiquitous phenomenon in early-type stars: stochastic
low-frequency photometric variability ("red noise"). It could come from internal gravity waves
launched by the convective core, or from subsurface convection. The paper shows that both the
characteristic frequency ν_char and the amplitude α₀ of the observed variability correlate well
with properties of the iron subsurface convection zone (FeCZ). The amplitude also correlates with
predicted core-excited IGW fluxes, but the observed trend of ν_char is at odds with core
turnover timescales. This suggests a unified picture in which red noise, surface turbulence,
magnetic spots and wind variability all originate in subsurface convection.

This folder is **archived**: its contents are as they were in April 2021, only moved here from
the repository root in September 2026. The original layout is preserved by the git tag
[`cantiello2021`](https://github.com/matteocantiello/rednoise/tree/cantiello2021), so links of
the form `github.com/matteocantiello/rednoise/blob/master/<file>` keep working as
`.../blob/cantiello2021/<file>`.

## Contents

| Path | What |
|---|---|
| `paper/` | Manuscript (`rednoise.tex`, AASTeX 6.3), bibliography, figures used in the paper, referee response |
| `figures/` | Figures produced by the notebooks (`figures/Old/`: earlier versions) |
| `rednoise.ipynb` | Main analysis: FeCZ and core properties vs. the red-noise sample (Bowman et al. 2019, 2020) |
| `kipp_rednoise.ipynb`, `kipp_plots.ipynb` | Kippenhahn diagrams |
| `propagation_diagrams.ipynb` | IGW propagation diagrams and wave fluxes |
| `rednoise_old.ipynb` | Earlier version of the analysis |
| `15140/Z0.02_CORE_F0/` | MESA r15140 grid setup (template, inlists, `run_grid.sh`) |
| `models/` | Earlier MESA setup: `template/` (incl. `run_star_extras.f`), `inlist_common`, the column lists and `inlist_pgstar` as used in 2021, `run_grid.sh` |
| `magnetic_models/` | MESA work directories for 3 and 30 Msun models used for the magnetic-field estimates |
| `notes/` | Handwritten notes (May 2020) |

## Notes
- The notebooks read MESA output from absolute paths on the author's laptop (`/Users/mcantiello/Dropbox/...`).
  That output is in the Zenodo archive, not in this repository. The notebooks save figures to `./figures/`,
  which is why they sit next to `figures/` here.
- The 2021 MESA templates predate the 2026 fixes in `../models/template_v2` (see `../models/GRID_LOG.md`).
  In particular, `conv_vel` includes rotational mixing in rotating models. The 2021 grids were non-rotating,
  so this did not affect the paper.
