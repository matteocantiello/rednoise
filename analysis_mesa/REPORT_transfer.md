# FeCZ → surface transfer models: what tracks v_macro, α₀ and ν_char

*2026-09-24. Grid: MW, ω/ω_c = 0 (α_MLT = 1.6, no MLT++), main sequence. Numbers from
`data/transfer_tests_MW_w0.0.csv` and `data/scenarios_MW_w0.0.csv`; figure `figures/fig_transfer.png`.*

## Method
- `transfer_models.py`: ~60 candidate predictors. Each is a closed-form function of MESA history
  quantities whose exponents come from the physical argument: velocities (MLT, Cantiello+09 microturbulence,
  wave-flux models with Goldreich–Kumar / Lecoanet–Quataert / Lighthill efficiencies, carried at c_s or at v,
  saturated at the photospheric sound speed), amplitudes (flux fractions, wave-flux fractions, surface Mach
  numbers, incoherent-cell dilution, depth-attenuated fluxes) and frequencies (turnover and crossing times,
  thermal time of the overlying layer). Controls: core-excited IGW flux, core turnover, dynamical frequency,
  ν_max scaling.
- `test_transfer.py`: the models are interpolated to each observed MS star's (log T_eff, log ℒ) position
  (IACOB v_macro N = 526; primary red-noise sample α₀ and ν_char N ≈ 175). Fit log y = log(floor + A·X^β)
  (ν: log y = a + β log X), with β fixed at 1 ("linear") or free. Each fit is compared by BIC with the
  empirical plane log y = c₀ + c₁ log ℒ + c₂ log T_eff (3 params). ΔBIC < 0 means better than the plane.
- `scenarios.py`: complete scenarios in which one physical chain predicts all three observables with β = 1.
  The predictions at the observed stars are then put through the same regressions as the data, which gives
  the implied sensitivity ratios.

**Caveat on free exponents.** On the MS, (T_eff, ℒ) fixes (M, τ), so any two model quantities with free
exponents can reproduce the observed plane. Only fixed-exponent tests and the nonlinear features (onset,
saturation) discriminate.

## Results
**Best joint scenario (S3t): wave flux F_w = 𝓜_t F_c, launched at the top of the FeCZ.**
The Mach number over the top α H_P multiplies the convective flux (Goldreich–Kumar efficiency applied to the
enthalpy flux).
- **v_macro:** set by ½ρ_s v³ = F_w, saturated at the photospheric sound speed:
  v = [(2F_w/ρ_s)^(−2/3) + c_s^(−2)]^(−1/2). Normalisation 10^−0.03 (≈ 0.93, i.e. no tuning); floor 23 km/s.
  - ΔBIC +0.8 against the plane over all MS stars (−25.6 with a quadrature floor).
  - Hot sample: free exponent β = 1.03 [0.93, 1.12], consistent with the physical value of 1.
- **α₀ ∝ F_w/F\*:** β = 1, ΔBIC +5.5 (free β = 1.41 [1.16, 1.65]); α₀ = 10^5.18 µmag × F_w/F\*.
- **ν_char = 10^0.45 × 1/(2πt_c):** ΔBIC −27.2 (free β = 0.98 [0.83, 1.14]).
- **Joint:** ΔBIC −20.9 against the three empirical planes, with 5 parameters against 9.
- **Implied sensitivity ratios (bootstrap 16–84%):**

  | | along ℒ | along τ | along M |
  |---|---|---|---|
  | S3t | 4.81 [4.46, 5.24] | 7.40 [6.62, 8.39] | 1.90 [1.70, 2.12] |
  | observed | 5.3 [4.9, 5.8] | 6.4 [4.9, 8.6] | 2.2 [1.8, 2.7] |

  The mechanism: the velocity saturates at c_s while the wave energy flux, and hence the brightness
  perturbation, keeps rising.

**Other scenarios (joint ΔBIC):**

| scenario | joint ΔBIC |
|---|---|
| S1 MLT + incoherent-cell dilution | −1.2 |
| S0 bare MLT | +4.2 (ratios 3.7 / 5.4 / 2.15) |
| S3 (𝓜_max instead of 𝓜_top) | +15.6 |
| S4 undamped core IGW | +97 |
| S2 Cantiello+09 microturbulence | +466 |

- **S4 core IGW:** fits v_macro, but its predicted ℒ-slope of α₀ is 0.46 (observed 1.55), its T_eff slope of ν
  is −0.63 (observed +1.53), and its ratios are 1.4 / 1.7 / 1.2.
- **S2 C09:** fails because ρ_s enters with the wrong mass dependence.

**Empirical facts found along the way**
- v_macro ≈ 1.5 c_s(photosphere) for all MS stars (log v_macro = a + 0.97 log c_s; rms 0.221 against 0.218
  for the plane). v_macro/c_s correlates only weakly with FeCZ properties (ρ ≤ 0.31, with 𝓜_t).
- On 95 stars with both, α₀ rises with c_s at fixed v_macro (α₀ ∝ v_macro^0.65 c_s^1.6). It does not behave
  as (v_macro/c_s)^n.
- Amplitude models built on surface density (C09, wave flux carried at c_s, surface Mach numbers) all fail
  (ΔBIC > +100 for α₀).

## Consequences for the §5.6 draft
- (iii) needs revising. Once the observed floor is included and predictions are evaluated at the observed
  stars, even bare MLT gives ratios of 3.7 / 5.4 / 2.15, and S3t gives 4.8 / 7.4 / 1.9. The earlier
  "model low by 1.6–2.5×" compared model-row slopes without the floor.
- (ii) the saturation is reproduced once v is capped at the photospheric c_s.

## Problems found
- **Rotating grids (w0.2, w0.4): FeCZ quantities are numerically unusable.** The median zone thickness is
  0.001 R (against 0.014 R non-rotating), 12–25% of hot MS models have a single-cell FeCZ, and the
  step-to-step jitter in log v_c is 0.3–0.5 dex (0.002 non-rotating). Likely cause: `run_star_extras`
  locates zones from MESA mixing regions, which are labelled by the dominant mixing type, and rotational
  diffusivities outcompete the weak FeCZ's convective D. The fix is to identify convective cells by
  gradr > grad_ad (or mlt_mixing_type) instead. This needs a rerun or a new run_star_extras; not yet done.
- **Only 1–5 MS profiles per track below 25 M☉** (max_num_profile_models = 100), so IGW radiative-damping
  integrals between the FeCZ and the photosphere cannot yet be computed.
- The α₀ scatter (rms 0.39–0.40 dex) is twice the 0.19-dex measurement systematic; the models leave half
  the variance unexplained.

---
**Caveat (added 2026-09-24):** the rotating-grid numbers above (MW w0.2/w0.4) come from v1 models whose
convective velocities include rotational mixing: MESA's `conv_vel` = 3·D_mix/(αH_P) with D_mix ⊃ D_rot.
At 20 Msun with w0.2, v_FeCZ is inflated ×2.6 and the turnover time ×0.27. The non-rotating (w0.0) results
are unaffected. Fixed in template_v2 / grids_v2; see `models/GRID_LOG.md`.
