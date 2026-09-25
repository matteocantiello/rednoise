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
