# Outline — "Surface turbulence and stochastic low-frequency variability in massive stars: a common sub-surface driver"

**Pitch (one sentence).** Across 342 OB-to-supergiant stars with red-noise fits and 832 with
macroturbulence, photometric amplitude and turbulent velocity are organized by the same
spectroscopic-HRD coordinate, switch on at the same luminosity, grow together along the main
sequence (amplitude 50×, velocity 1.9×), and respond to rotation the same way — the signature
of a shared sub-surface convective driver, not core-excited gravity waves.

**Audience.** Massive-star / asteroseismology community (ApJ).

**Most arresting asset.** The GP posterior maps on the sHRD (Fig. 3): a continuous, data-driven
picture of *expected* v_macro, α₀, ν_char, γ at every HRD position.

## Section plan

1. Introduction — SLF variability; IGW vs FeCZ predictions; the sign/onset/evolution
   diagnostics; what a meta-catalog adds.
2. Data — 2.1 red-noise samples; 2.2 macroturbulence samples; 2.3 harmonization
   (ν_char units, four α₀ conventions, spectroscopic vs classical L); 2.4 placing stars on
   the sHRD (tier A log g, tier B MIST mass, calibration); 2.5 evolutionary parameters;
   2.6 excluded data and corrections (Shen vmacro column, Holgado kK, Dorn-Wallenstein α₀).
3. Methods — 3.1 correlations, model selection, bootstrap; 3.2 reliability weights
   (measured systematic floor, rotation penalty); 3.3 GP maps and the trustworthiness
   mask; 3.4 broken-power-law onset test; 3.5 rotation-frequency test.
4. Results — 4.1 two orthogonal HRD coordinates (Fig. 2, 3); 4.2 luminosity onset and
   saturation (Fig. 4); 4.3 evolution at fixed mass (Fig. 5); 4.4 rotation (Fig. 6);
   4.5 the cool side: ν_char floor and RSG regime; 4.6 metallicity (method-limited).
5. Discussion — 5.1 FeCZ vs core IGWs (Fig. 7); 5.2 the factor-6 amplitude-to-velocity
   sensitivity: α₀ as a flux ratio; 5.3 the sub-onset floor; 5.4 rotation as a modifier
   of sub-surface convection; 5.5 **[PLACEHOLDER] Comparison with MESA models**;
   5.6 caveats.
6. Conclusions (enumerated).
Appendix A — unit systems and cross-validation (Bowman vs Shen overlap; tier-B calibration).
Appendix B — correlation matrix; machine-readable table description.

## Figure arc (main text) — claim each figure carries

| Fig | claim | source data |
|---|---|---|
| 1 | The compiled samples cover the upper sHRD from ZAMS to supergiants, with two placement tiers | rednoise_sHRD_extended_evol.csv, macroturbulence_evol.csv, MIST |
| 2 | On the homogeneous 196-star sample, ν_char follows T_eff and α₀ follows ℒ and v_macro | master_rednoise_catalog_v2 |
| 3 | Expected v_macro, α₀, ν_char, γ at every HRD position: α₀ and v_macro share a coordinate, ν_char is orthogonal | hrd_interpolated_grids_v2 |
| 4 | Turbulence and amplitude switch on at ℒ≈3.0–3.2, above the model FeCZ appearance; v_macro saturates at 3.7 | fecz_onset_test.csv |
| 5 | At fixed mass, α₀ grows 50× across the MS, v_macro 1.9×, ν_char declines | evolution_test_full.csv |
| 6 | ν_char ≈ 17 f_rot: rotation is not the clock but a weak common multiplier on all three observables | rotation_test.csv |
| 7 | Scorecard + amplitude budget + added-variable test: FeCZ matches, core IGWs fail | fig_v2_theory inputs |

Appendix figures: A1 Bowman-vs-Shen overlap validation; A2 tier-B ℒ calibration; B1 correlation heatmap.

## Kill list (not in paper)
- Metallicity violin figure (method-confounded; a paragraph suffices).
- Outlier annotation figure (folded into a table of flagged stars).
- Multivariate coefficient lollipop (numbers go in Table 2).
