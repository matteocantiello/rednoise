# Microturbulence (Markova, Cantiello & Grassitelli 2025) on the sHRD

*2026-09-25. Code: `vmic_test.py` → `data_obs/vmic_test.csv`, `figures_v2/fig_vmic.png`. Data: J/A+A/701/A297 (tables E1
non-pulsating, E2 pulsating); log ℒ uses our normalisation (checked: HD 14947 gives 4.12 as tabulated).*

Sample: 443 stars with log T_eff > 4.0 (69 pulsating); 27 with v_mic = 0 are excluded from the log fits.

## Results
1. **No onset.** v_mic rises smoothly with ℒ, from ≈1 km/s at log ℒ 1.5 to ≈15 km/s at 4.3, flattening above ≈4.1. The hinge fit
   (log T_eff covariate) prefers a single power law in every subset (ΔBIC −7 to −9).
2. **Magnitude matches the model FeCZ velocity in the middle of the range, with no normalisation.** At log ℒ ≈ 3.2–3.8 the median
   v_mic (5–15 km/s) equals the median model v_c,max at the same stars (MW v2, non-rotating). Below that the model is lower (no
   FeCZ at all below log ℒ ≈ 2.6, where v_mic is still 1–5 km/s); above, the model keeps rising to ≈60 km/s while v_mic flattens
   at ≈15 km/s. As a predictor with a floor (β = 1), model v_c,max does worse than the empirical plane (ΔBIC +42; Spearman 0.73;
   free slope 0.41).
3. **v_mic and v_macro share star-to-star fluctuations.** 143 stars in common with IACOB: Spearman +0.58, and **+0.35 after both
   are regressed on log ℒ and log T_eff**; log v_mic ∝ 1.07 log v_macro.

## Reading
- Microturbulence behaves like a weaker, smoother version of macroturbulence: same ℒ ordering, same flattening at the top,
  correlated residuals, but no threshold and a non-zero level below the FeCZ onset. That fits a picture in which the FeCZ adds to
  a background small-scale velocity field (Markova+2025 conclude v_mic is physical, not a modelling artefact), with the FeCZ
  contribution dominating only above log ℒ ≈ 3.
- The agreement in magnitude with MLT v_c,max at log ℒ 3.2–3.8, and its failure above, mirror the macroturbulence result: MLT
  velocities describe the photospheric field where the Fe zone is sub-Eddington and weakly turbulent, not near Γ_Fe ≈ 1
  (REPORT_regime.md).

## Update 2026-09-25 evening: the full sample on the sHRD, and every subsurface zone
*Code: `make_fig_vmic_hrd.py` → `figures_v2/fig_vmic_shrd.png`, `figures_v2/fig_vmic_zones.png`, `data_obs/vmic_zones.csv`.
All 1823 Markova+2025 stars (1424 non-pulsating, 366 pulsating, 33 with v_mic = 0). 1380 are cooler than 10 kK and many
lie below the model grid (log ℒ < 2.3, below 5 M☉). The zone comparison uses the ~790 stars on the MW v2 grid (ω = 0, MS +
post-MS), with each zone's v_c,max interpolated between tracks and set to 0 where the zone is absent.*

| model zone | stars with the zone | Spearman ρ | median log(v_mic / v_c,max) | ρ hot (> 20 kK) | ρ cool (< 10 kK) |
|---|---|---|---|---|---|
| **Fe bump (FeCZ)** | 544 | **+0.78** | **0.00** | +0.78 | +0.44 |
| He II | 604 | −0.31 | +1.35 | +0.42 | −0.28 |
| He I | 50 | −0.22 | +4.2 (zone negligible) | — | −0.56 |
| H I | 374 | +0.34 | −0.44 | — | +0.34 |
| largest of all zones | 771 | +0.58 | −0.35 | +0.79 | +0.42 |

- **Where a star has an FeCZ, v_mic equals the model FeCZ convective velocity on average, with no normalisation** (median log
  ratio 0.00, ρ = +0.78 over 544 stars). This is the strongest model–data match in the project. The earlier profile test
  (above) found the model too high at the very top; that remains true for the most luminous O stars, where v_c,max reaches 50–100 km/s
  and v_mic saturates at ~15–20 km/s.
- He II and He I zones are far too weak in MLT (≲ 1 km/s and ≲ 10⁻³ km/s) to set v_mic. Their anticorrelation only
  reflects where in the HRD they exist.
- For the cool supergiants (< 10 kK) the H I recombination zone is vigorous in the models (5–15 km/s), but v_mic there is
  only 1–7 km/s (median −0.44 dex) and weakly correlated (ρ = +0.34). The H I panel has two clumps. At v_mic ≈ 3–7 km/s there are
  251 stars at log T_eff ≈ 3.77 and model v_c ≈ 12 km/s, 75% of them pulsators (table E2). At v_mic ≈ 1–2.5 km/s there are 89
  cooler, non-pulsating stars (log T_eff ≈ 3.64, 1% pulsators) with model v_c ≈ 4 km/s.
- Below log ℒ ≈ 2.5 (outside the grid, and without an FeCZ at MW composition) v_mic falls to ~1 km/s, and to ≲ 0.5 km/s in
  the coolest dwarfs.
