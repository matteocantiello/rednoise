# Homogeneous MW vs SMC ν_char: one pipeline on both galaxies

*2026-09-25. Code: `tess_mw_download.py` (SPOC 2-min, ≤ 4 sectors per star, Ceph `rednoise/tess_mw/`), `vandaele_fit.py`
(tags `_mw` PDCSAP, `_mwdeg` PDCSAP with white noise added to reach the SMC median α₀/C_w, `_mwsap` SAP), `mw_smc_nuchar.py` →
`data_obs/mw_smc_{nuchar,validation}.csv`, `smc_instrumental_check.py` → `data_obs/smc_instrumental_check.csv`.*

## Galactic light curves
197 targets (Bowman 2020 + Shen 2024, duplicates removed): 171 with SPOC 2-min data, 532 sectors; 26 have none.

**Validation against the published values** (linear / log fit metric):
- Shen 2024, same sector: median Δlog ν −0.02 / +0.07, MAD 0.20 / 0.15, Spearman 0.69 / 0.77 (N = 227 / 242).
- Bowman 2020, per star (sectors ≤ 13): −0.15 / −0.06, MAD 0.09 / 0.12, Spearman 0.75 / 0.82 (N = 64).

Our pipeline reproduces the Galactic literature at the level of the 0.19 dex study-to-study systematic.

## The raw comparison
Plane log ν = c0 + c1 log ℒ + c2 log T_eff, fitted on our Galactic per-star medians (N ≈ 166). SMC stars are placed with
Bestenlehner+2025 (T_eff, log g → ℒ). **54 of the 91 match by position.** The other 37 (mostly cool supergiants) are not in the
OB table. A bug that gave them the last catalogue row's parameters was fixed on 2026-09-25, so the numbers below supersede the
first version of this report. SMC stars inside the Galactic domain: 30.

| MW variant / metric | SMC median residual [95%] | p (MWU) |
|---|---|---|
| PDCSAP / lin | −0.22 [−0.34, −0.10] | 4e-5 |
| PDCSAP / log | −0.33 [−0.47, −0.23] | <1e-6 |
| PDCSAP noise-degraded / lin | −0.22 [−0.33, −0.09] | 8e-5 |
| SAP / lin | −0.15 [−0.29, −0.01] | 0.03 |
| SAP / log | −0.29 [−0.42, −0.18] | 6e-5 |

This is a deficit of SMC ν_char at matched position, in the direction the FeCZ models predict (LMC grid: −0.16 dex; SMC grid
pending). It depends on the metric and flux type by ±0.1 dex. Adding white noise (degrade test) changes nothing.

## Instrumental check: the PSF noise floor
- **At the periodogram level the floor looks like SLF.** SMC stars that Van Daele class as showing no significant variability
  (38 stars, 104 sectors) pass our cuts as often as their SLF stars (74% vs 76%). They return a similar ν_char (log −0.50 vs
  −0.44) and only 0.25 dex lower α₀/C_w. In e⁻/s the floor amplitude scales roughly with flux (α₀ ∝ F^1.2), with 0.65 dex
  scatter and sector offsets (sector 27 about 0.5 dex quieter). Van Daele rejected 139 of 281 sectors for low-frequency excess.
- **At matched sHRD position the non-variable stars do *not* show the deficit.** The 12 placed non-variable stars in the
  Galactic domain sit at +0.11 (PDCSAP / lin), against −0.20 for the SLF stars. N is small, and those stars are brighter and
  hotter on average, so this is not conclusive. It does mean the SLF-star deficit is not simply the floor.
- Residuals do not correlate with Tmag. SAP Galactic fits move 0.05–0.1 dex towards the SMC values.
- Reading: the SLF detections are not cleanly separable from the floor, and the floor itself has ν_char ≈ 0.3 d⁻¹, below the
  Galactic plane at these positions. A floor contribution can therefore bias SMC ν_char low, but the data do not show that it
  produces all of the deficit. → noise-floor model and injection–recovery (`smc_floor_model.py`, `REPORT_floor.md`).

## Consequences
- **Not yet a clean comparison.** With our pipeline the SLF detection in the SMC PSF light curves is not separable from the
  noise floor, so neither their "similar morphology" nor our deficit can be taken at face value without modelling the floor.
  Their detection/non-detection split may rest on their visual screening more than on the fit.
- What a valid test needs: (1) an empirical noise-floor model from the "no variability" stars, background apertures or
  nearby faint sources in the same FFIs, fitted jointly (SLF + floor + white); or (2) Galactic comparison light curves made the
  same way (FFI PSF/SAP photometry of fainter Galactic OB stars, not SPOC 2-min of bright ones). With either, ν_char comparisons
  only for stars whose SLF sits clearly above the floor.
- For the paper: the SMC SLF samples do not yet constrain the metallicity dependence of ν_char. The published claim of
  metallicity-independent SLF morphology inherits this limitation. This is independent of, and adds to, the selection-bias
  argument in `Vandaele_correspondence/README.md`.
- It is also a concrete, collegial point to raise with Van Daele (it bears on their SLF/no-SLF classification) and a reason to
  ask for their per-sector fits and settings.
