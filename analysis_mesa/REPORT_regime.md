# FeCZ regime parameters from the 3D literature, on the v2 grid

*2026-09-25. Step 1 of the plan in WAVE_MODEL.md section 7. MESA v2, MW, ω = 0 (non-rotating only: the grid
profiles store conv_vel, which includes rotational mixing). Code: `regime_params.py` →
`data_v2/regime_MWw0.0.csv` (1079 triggered profiles, 34 masses); `nu_regime_test.py`; figure
`make_fig_regime.py` → `figures_v2/fig_regime.{png,pdf}`.*

## Definitions (at the Fe opacity peak = max κ with 1e5 < T < 5e5 K)
- Γ_Fe = κ_Fe L(r)/(4πGMc): local Eddington factor.
- τ_crit = c P_rad/((P_rad + P_gas) v_c) (Jiang+2015); τ_Fe/τ_crit < 1 means radiatively lossy plumes.
- Y_Fe = L/(4πr² a T⁴)/c_iso: pseudo-Mach number (Schultz+2020).
- MLT v_c,max and Mach number in the FeCZ; thermal times ∫c_pT dm/L of the zone, above the Fe peak, and above
  the zone top (c_pT = Pδ/(ρ∇_ad), δ = (4 − 3β)/β).

## Results
**1. The v_macro saturation (log ℒ_spec = 3.70) sits where the Fe peak nears the local Eddington limit.**
Median over hot MS models, with track-bootstrap 16–84%:

| threshold | log ℒ_spec |
|---|---|
| Γ_Fe = 0.8 | 3.66 [3.65, 3.66] |
| Γ_Fe = 0.9 | 3.76 [3.75, 3.77] |
| Γ_Fe = 1.0 | 3.87 [3.86, 3.88] |
| MLT Mach = 0.2 | 3.76 [3.76, 3.78] |

The 3D models tell the same story: M13TAMS (Γ_Fe = 0.55) is weakly turbulent and matches 1D; M35ZAMS (0.82) is
quiet because radiatively lossy (τ ratio 0.02); M35MMS (0.97) is transsonic with ≈150 km/s at the photosphere.
So the saturation marks the transition to strongly turbulent, near-Eddington Fe zones, where MLT (whose
velocities keep rising) stops describing the photospheric velocity field.
*Caveat:* Γ_Fe ≈ (κ_Fe/κ_es) Γ_es and Γ_es ∝ ℒ_spec, so Γ_Fe tracks ℒ_spec nearly by construction; the
content of the result is the value (Γ_Fe ≈ 0.8–0.9, no free parameter) and its **prediction: the saturation
should move to higher ℒ at lower metallicity**, where κ_Fe is smaller (testable with LMC/SMC grids and
macroturbulence samples). Γ_Fe and the MLT Mach number cross their thresholds together within the MW grid;
metallicity would separate them.

**2. The onset (3.16 [2.94, 3.44]) coincides with several thresholds that do not discriminate:**
Γ_Fe = 0.4 at 3.10, MLT v_c = 3 km/s at 2.99, Y_Fe = 1 at 3.30 [3.21, 3.36].

**3. τ_Fe/τ_crit is not a hinge variable on the MS:** it stays at 0.06–0.18 (median) across the whole MS
range. All models are in the radiatively lossy regime, like the 3D MS models (0.02–0.25).

**4. Frequencies:** (165 of 183 primary red-noise stars covered; plane rms 0.401 dex)

| predictor | β = 1: offset, ΔBIC | free β, ΔBIC |
|---|---|---|
| ν_turn = v_c/(2π α H_P) | 10^0.22, +11.2 | 0.43, −6.9 |
| ν_th above Fe peak | 10^0.43, +194 | 0.09, +2.6 |
| ν_th above zone top | 10^−0.84, +134 | 0.17, −2.6 |
| ν_th zone | 10^1.15, +252 | 0.04, +5.5 |
| ν_turn^b1 · ν_th^b2 | — | b2 < 0 (wrong sign), −3 to −5 |

The thermal times are nearly constant across the HRD and do not track ν_char, so interpolating between
turnover and thermal frequency does not work in this form. The turnover frequency is the best single
predictor, but with β ≈ 0.4–0.5. The residual log(ν_char/ν_turn) correlates modestly with τ_Fe/τ_crit
(ρ = +0.34): ν_char sits further above turnover where the plumes are less lossy. The frequency information
is weak overall (the plane barely beats a constant).

## Next
- Γ_Fe as the saturation variable belongs in the paper (§5.6 ii) once v2 is complete; the prediction
  (saturation moves with Z) needs the LMC/SMC v2 grids.
- Rotating grids need the MLT velocity in profiles (`mlt_vc` column) or a history-based Γ_Fe/τ_crit
  (add to run_star_extras for the next grid generation).
- Frequencies remain open: the GYRE damped-mode continuum and direct 3D comparisons are the next tools.
