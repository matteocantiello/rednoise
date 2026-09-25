# Star-by-star rotation, full-coverage v2 scenarios, and the LMC regime shift

*2026-09-25. Grid v2 (template_v2, A09). Galactic extracts re-made at ~11:05 from the current runs (31–32 of
34 tracks with a complete MS per sub-grid; only 100–120 Msun still on the MS). The previous extracts, which
`REPORT_transfer_v2.md` used, are kept in `data_v2/snap_20260924_2318/`. Code: `rotation_mixed.py`,
`scenarios.py` (new S5), `regime_crossings.py`. Tables: `data_v2/rotation_mixed{,_omega}.csv`,
`data_v2/scenarios_MW_w*.{csv,txt}`, `data_v2/regime_{MW,LMC}w0.0_crossings.csv`.*

## 1. Coverage matters
With the 2026-09-24 extracts, the stars covered by all four ω maps drop to 349 v_macro / 116 α₀ / 113 ν_char
(from 460 / 166 / 163), because the massive rotating tracks were still on the MS. That alone flipped the
core-IGW verdict (α₀ ΔBIC +16 instead of +80). On the new extracts the common set is 516–522 / 177–178 /
175–178 stars. **Everything below uses the new extracts.**

## 2. Star-by-star rotation does not improve the fits
Per-star ω: v_eq = v sin i · 4/π, matched to the model surface velocity at the star's position, clipped to
[0, 0.6]. Median ω = 0.28–0.30 (16–84%: 0.06–0.60); 31–38% of stars sit at the 0.6 cap.

Joint ΔBIC against the empirical ℒ + T_eff planes (β = 1; same stars in every column):

| scenario | ω = 0 | 0.2 | 0.4 | 0.6 | per-star ω |
|---|---|---|---|---|---|
| S0 bare MLT (v_c,max; F_c/F; 1/(2πt_c)) | +14.8 | +8.4 | −0.6 | −8.1 | +6.5 |
| S1 MLT, cells-diluted flux | +7.4 | +0.9 | −8.1 | −14.9 | −1.6 |
| S3t waves F_w = 𝓜_t F_c, v saturated at c_s | +17.6 | +11.5 | +3.9 | +0.3 | +14.3 |
| **S5 v_c,max; 𝓜_t F_c/F; 1/(2πt_c)** (post hoc) | **−13.0** | **−19.7** | **−28.4** | **−34.8** | **−20.9** |
| S4 core IGW (undamped) | +87.6 | +88.7 | +93.0 | +95.5 | +91.8 |
| S2 Cantiello+09 microturbulence | +537 | +534 | +527 | +521 | +548 |

For every FeCZ scenario the per-star result lies between the ω = 0 and 0.2 grids. The trend "better with more
rotation" therefore does **not** follow the stars' own rotation: slowly rotating stars also fit better with the
ω = 0.6 models. Whatever the rotating models capture (their envelope structure, rotational mixing), it is not
the rotation of the individual star. The single-ω improvement should not be read as a rotation effect.

## 3. Per observable, and the S5 combination
Best single predictor per observable (per-star ω, ΔBIC): v_macro ← MLT v_c,max (−17.8); α₀ ← 𝓜_t F_c/F (+6.8,
i.e. as good as the 3-parameter plane with 2); ν_char ← 1/(2π t_c) (−10.0, free β = 0.67).
S5 combines them. It is chosen **post hoc**, but it has a physical reading from the 3D envelopes (Schultz+2023:
photospheric velocities ≈ FeCZ convective velocities; the brightness perturbation set by the Mach-scaled
convective flux). At ω = 0 its implied sensitivity ratios are 4.31 [3.91, 4.79] / 7.01 [6.12, 8.11] /
1.82 [1.62, 2.06] along ℒ / τ / M (observed 5.3 [4.9, 5.8] / 6.4 [4.9, 8.6] / 2.2 [1.8, 2.7]).

**Caveat:** on the hot sample (log T_eff > 4.3), where the onset and the saturation are, S5's v_macro fit is much
worse than the plane (ΔBIC +51 to +64), because v_c,max keeps rising through log ℒ = 3.7. The saturated-wave
velocity of S3t does better there (+33). This is consistent with the regime result (REPORT_regime.md): above
Γ_Fe ≈ 0.8 the MLT velocity no longer describes the photospheric field.
Core IGWs fail in every mode (α₀ ΔBIC ≈ +105, ν +16), though they fit v_macro alone well (−30).

## 4. Metallicity: the LMC grid moves the thresholds by different amounts
Median crossings over hot MS models, track bootstrap 16–84% (`regime_crossings.py`; the MW row reproduces
REPORT_regime.md to within 0.01–0.03):

| threshold | MW (Z = 0.014) | LMC (Z = 0.006) | shift |
|---|---|---|---|
| v_c = 3 km/s (onset proxy) | 3.02 [2.90, 3.05] | 3.62 [3.58, 3.65] | +0.60 |
| Γ_Fe = 0.4 | 3.11 | 3.34 | +0.23 |
| Γ_Fe = 0.8 (saturation proxy) | 3.65 | 3.89 [3.88, 3.89] | +0.24 |
| Γ_Fe = 0.9 | 3.76 | 3.99 | +0.23 |
| Γ_Fe = 1.0 | 3.88 | 4.10 | +0.22 |
| MLT Mach = 0.2 | 3.75 | 4.07 [4.04, 4.09] | +0.32 |
| Y_Fe = 1 | 3.29 | 3.12 | −0.17 |

**Prediction:** at LMC metallicity the v_macro/α₀ onset should move up by ≈0.6 dex in log ℒ and the saturation by
≈0.24 dex, so the rising segment between them narrows from ≈0.5 to ≈0.3 dex. The variables that coincide in the
Galaxy separate here: a saturation set by Γ_Fe moves by 0.24, one set by the MLT Mach number by 0.32, and Y_Fe
(no opacity dependence) moves the other way. The existing LMC red-noise stars (Bowman 2019b/2024, Ma 2024) all lie
at log ℒ ≥ 3.3 (most ≥ 3.8), and their amplitudes are not on the µmag scale, so they cannot test this yet. LMC
macroturbulence from the VLT-FLAMES surveys would, if a homogeneous set exists (not checked).

## Next
- SMC v2 regime crossings when the SMC grid (queued) has main-sequence models.
- The composition bracket (A09 at Z = 0.010 and 0.018, non-rotating, 10–30 Msun, MS only; local runs in
  `/home/mcantiello/rednoise_tests/zbracket_Z0{10,18}`, running since 11:00, 30 Msun ≈ 11 h): how far do the
  crossings and the S5 fits move for a ±0.13 dex change in Galactic Fe?
