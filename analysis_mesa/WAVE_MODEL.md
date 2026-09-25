# FeCZ-excited travelling waves: from scalings to a spectral forward model

Working document for the physical model behind the transfer-model results in `REPORT_transfer.md`.
Newest status at the bottom. Started 2026-09-24.

## 1. Where the scaling results stand (MW, ω = 0, main sequence)

Best joint scenario **S3t**: waves launched at the top of the FeCZ with energy flux F_w = 𝓜_t F_c
(𝓜_t = Mach number over the top α H_P, F_c = maximum convective flux of the FeCZ).

| observable | scaling | fit |
|---|---|---|
| v_macro | ½ ρ_s v³ = F_w, capped at the photospheric c_s: v = [(2F_w/ρ_s)^(−2/3) + c_s^(−2)]^(−1/2) | normalization 0.93, free exponent 1.03 |
| α₀ | ∝ F_w / F* | δL/L ≈ 0.14 · F_w/F* (α₀ = 10^5.18 µmag · F_w/F*); free exponent 1.41 [1.16, 1.65] |
| ν_char | 10^0.45 × 1/(2π t_c), t_c = α H_P / v_c | free exponent 0.98 |

- Joint ΔBIC −21 against the empirical planes (5 parameters against 9).
- Sensitivity ratios along ℒ / τ / M: model 4.8 / 7.4 / 1.9, observed 5.3 / 6.4 / 2.2.
- Core-excited IGWs fail (wrong α₀ ℒ-slope, wrong sign of the ν T_eff-slope, ratios 1.4 / 1.7 / 1.2).

**Limits.** On the MS, (T_eff, ℒ) ↔ (M, age), so free exponents fit almost anything; only fixed
exponents and the nonlinear features discriminate. Three numbers are *fitted*, not derived: the 2.8 in
ν_char, the 0.14 in α₀, and the c_s cap. α₀ prefers exponent 1.4 and leaves half of its variance
unexplained. The basis is MW non-rotating only; the v1 rotating grids are biased (see
`../models/GRID_LOG.md`).

## 2. Physical chain and equations

**Excitation** (Goldreich & Kumar 1990; Lecoanet & Quataert 2013):
- F_w ≈ 𝓜 F_c for a sharp convective/radiative interface (weaker 𝓜 dependence for a smooth one).
- Spectrum peaks at ω_c = v_c / (α H_P), falls roughly as dF/d ln ω ∝ (ω/ω_c)^(−13/2) above it (sharp
  interface). Horizontal scale ~ eddy size: ℓ_eddy ≈ r / (α H_P) ≈ 100–160 at 20 Msun.

**Propagation** through the radiative layer between the FeCZ top and the photosphere (~5.5 H_P, `FeCZ_nHP_top`):
- WKB dispersion: k_r² = k_h² (N²/ω² − 1), k_h = √(ℓ(ℓ+1)) / r.
- Wave luminosity: 4π r² ρ v_h² v_g,r e^(τ) = const, v_g,r ≈ ω² / (N k_h) for N ≫ ω. The amplitude grows
  ~ ρ^(−1/2) outward (×~15 over the layer), so the waves become nonlinear (k_r ξ_r → 1) near the surface.
- Radiative damping (Zahn, Talon & Matias 1997; coded in `template_v2` as `igw_tau`):
  τ(ω, ℓ) = [ℓ(ℓ+1)]^(3/2) ∫ K N³ / (ω⁴ r³) · √(N²/(N²−ω²)) dr,  K = 16σT³ / (3κρ²c_p).
  Evanescent zones (N < ω) add τ_ev = 2 ∫ k_h √(1 − N²/ω²) dr (negligible in the tests: ≲ 0.01).
  τ ∝ ℓ³ ω⁻⁴: the overlying layer is a high-pass filter in ω and a low-pass filter in ℓ.

**Surface response:**
- Velocity: nonlinear saturation near the surface is presumably what the fitted "½ρv³ = F_w, capped at
  c_s" stands in for.
- Brightness: at the photosphere the thermal time ≪ wave period, so the response is non-adiabatic
  (Dupret et al. 2003; Townsend 2002, both in `../literature/`). Reading of α₀ ∝ F_w/F*: the wave energy
  flux is converted to emergent-flux fluctuations; 0.14 = the fraction surviving disk integration
  (favours low ℓ).
- Travelling, not standing: with τ ≫ 1 per round trip at the relevant frequencies, no standing modes
  form between the FeCZ and the surface, so the spectrum is continuous (red noise). Prediction: coherent
  peaks could appear only at the highest frequencies, where τ < 1.

## 3. Hint from the damping columns (test T1a, 20 Msun, ω = 0, `/home/mcantiello/rednoise_tests`)

| X_c | ν_c = 1/(2π t_c) | 2.8 ν_c (fitted ν_char) | ν_damp ℓ=1 | ν_damp ℓ=5 | ν_damp ℓ_eddy | τ(ω_c, ℓ=1) |
|---|---|---|---|---|---|---|
| 0.68 | 0.48 d⁻¹ | 1.36 | 1.75 | 4.9 | 43 | 172 |
| 0.50 | 0.51 | 1.42 | 1.38 | 3.8 | 31 | 54 |
| 0.30 | 0.51 | 1.42 | 1.03 | 2.9 | 19 | 17 |

(ν_damp = frequency where τ = 1.)

- Waves at the turnover frequency are fully damped, even at ℓ = 1. The surface sees only the source
  spectrum's tail above ν_damp. A steep source × a sharp cutoff peaks near ν_damp, which could explain the
  fitted factor ~3 between ν_char and ν_c. On this track, ν_damp(ℓ=1) ≈ 2.8 ν_c.
- ν_damp falls with age while ν_c is flat: a testable prediction for the observed age trend of ν_char.
- **Caveats:** one track; the quasi-adiabatic τ is least reliable just below the photosphere, where much of
  the damping may accumulate.

## 4. Plan

1. **Grid test of the damping cutoff** (when v2 MW has MS models): regress observed ν_char on
   ν_damp(ℓ = 1…3), ν_c and both, in the `test_transfer.py` BIC framework; check the age trend against
   `../handoff/project_handoff/data/evolution_test_full.csv`.
2. **Spectral forward model** (Python, on profiles; v2 writes one every ΔX_c = 0.05):
   P(ν) = Σ_ℓ S(ω, ℓ) · e^(−τ(ω,ℓ)) · V_ℓ² · R(ω). Fit P(ν) with the same red-noise profile as the
   data → α₀, ν_char and the steepness γ with no fitted factors. γ is a new test (source slope × cutoff).
3. **Brightness response R(ω)**: non-adiabatic surface response (Dupret et al. 2003 style, or GYRE) to
   replace the fitted 0.14.
4. **Saturation**: evaluate k_r ξ_r along the path; saturate where it reaches ~1; check v_macro ≈ 1.5 c_s.
5. **Calibration** against 3D radiation-hydro envelope models (Jiang et al.; Schultz, Bildsten & Jiang
   2022) for the excitation efficiency and spectral shape. Anders et al. (2023) found core IGWs too weak
   photometrically.
6. **Rotation and Z** once the v2 grids exist; the SMC stars below the model FeCZ onset (12 of 23) are a
   direct test.

## 5. Status log

- **2026-09-24:** plan written. Starting step 2: code in `wave_spectrum/` (profile reader, FeCZ locator,
  damping integral validated against `run_star_extras`, source spectra, visibilities, red-noise fit),
  developed on the local T1 profiles.

### 2026-09-24 — step 2 scaffold built; first forward spectra (20 Msun, T1a, ω = 0)

**Code** (`wave_spectrum/`, run from `analysis_mesa/` as `python3 -m wave_spectrum.<script>`):
- `mesa_io.py`: profile/history reader; FeCZ locator mirroring run_star_extras (MLT-convective runs above
  T = 1e6 K, type set by the T of the bottom cell in (1e5, 5e5) K, max L_conv/L fragment); `Track.pairs()`
  pairs each MS profile with its history row. The source (ω_c, 𝓜_t, F_c) comes from HISTORY, because
  profile `conv_vel` includes rotational mixing and `mlt_vc` is not in the v2 profile columns (the v2 grid
  was already running, so the column list was left unchanged).
- `damping.py`: `tau` (same discretisation as `igw_tau`), `nu_damp`, `tau_profile` (cumulative),
  `q_nonad` (local non-adiabaticity q = K k_r²/ω), and the option `stop_at_q` (damp only below the layer
  where q reaches a threshold).
- `source.py`: parametric excitation spectrum: L_w = 𝓜_t F_c/F · L, frequency slope a (13/2) above ω_c,
  eddy-size cutoff ℓ_max = ℓ_eddy (ω/ω_c)^(3/2), low-ℓ slope p (unknown; key for the visible fraction).
- `visibility.py`: disk-integration weights W_ℓ (Eddington limb darkening; W_1 = 0.50, W_2 = 0.11,
  W_3 = 0.004).
- `spectrum.py`: PSD = Σ_ℓ W_ℓ R dL_w,ℓ/dν e^(−τ); red-noise fit α₀/(1 + (ν/ν_char)^γ).
- `validate_damping.py`: **the Python τ matches the MESA history columns exactly** (all printed digits,
  all T1a profiles).
- `demo_T1.py`: the scans below plus the figure `figures/wave_spectrum_T1.png`.

**Findings**
1. **The quasi-adiabatic damping integral is invalid where it matters.** 75% of τ accrues in the outermost
   pressure scale height below the photosphere (log T 4.5–4.6), where q = K k_r²/ω ≫ 1: median 12 at
   1 d⁻¹ and 450 at 0.3 d⁻¹, for ℓ = 1. Over 99% of τ at ν ≤ 1 d⁻¹ comes from cells with q > 1. There the
   wave's temperature fluctuations are radiated within a period: the lost wave energy becomes
   emergent-flux fluctuations (i.e. the photometric signal) rather than heat. Damping and the brightness
   response are the same process in that layer.
2. Damping to the photosphere (as coded in run_star_extras): the spectrum is band-pass, peaked at
   0.8–1.4 d⁻¹ (≈ ν_damp, ℓ=1), with zero power at 0.05 d⁻¹ and transmission 10⁻¹¹–10⁻⁸. It cannot be
   red noise; the ν_damp ≈ 2.8 ν_c coincidence of section 3 should be treated with caution.
3. Damping only below the transition layer (`stop_at_q` = 1–10): 4–40% of the wave flux reaches the
   layer, and the spectrum is essentially the assumed source spectrum. **Caveat, from the figure:** the
   apparent red-noise fit (ν_char ≈ ν_c, rms 0.02 dex) is mostly the sharp corner of the broken-power-law
   source at ν_c, and below ν_c the source rises as ν^(−1/2) instead of flattening. At `stop_at_q` =
   100–1000 the hard threshold creates dips and secondary bumps (the set of cells kept jumps with ω; at low
   ν nothing is damped). **The ν_char and γ numbers from the threshold scan are artifacts, not results.**

**Consequences / next**
- Step 3 moves ahead of step 2. The propagation through the outer ~1 H_P has to be solved
  non-adiabatically: linear non-adiabatic wave equations for a forced travelling wave (outgoing boundary
  condition), in the style of Dupret et al. 2003, or GYRE's non-adiabatic machinery. That gives, in one
  calculation, the damping, the conversion into emergent flux, and hence R(ω) and α₀ without the 0.14
  factor.
- Replace the broken-power-law source with a smooth spectrum (e.g. a Lorentzian/GK-like rollover at ω_c,
  with the low-frequency behaviour from Lecoanet & Quataert 2013) so that the fitted ν_char is not set by
  a kink.
- Tension to keep in mind: linear waves give amplitude ∝ (wave flux)^(1/2), but the scaling fit wanted
  α₀ ∝ F_w^1 (free exponent 1.41). Transmission that rises with F_w, or non-linearity near the surface,
  would be needed.

### 2026-09-24 — diagnostic figures for one profile

`python3 -m wave_spectrum.plot_profile [--run DIR] [--xc 0.5]` → `figures/wave_envelope_<tag>.{png,pdf}` and
`figures/propagation_<tag>.{png,pdf}` (first made for T1a, 20 Msun, X_c = 0.51: tag `M20_Xc0.51`).
- **Envelope figure** (x = log T, surface on the right; FeCZ/HeII shaded; dashed line = FeCZ top, where
  waves are launched; dotted line = τ = 2/3), eight panels: (a) κ (Fe and HeII bumps); (b) ∇_rad vs ∇_ad;
  (c) v_conv (MLT) and c_s (𝓜_t = 0.022, v_max = 4.2 km/s); (d) H_P/R; (e) L_conv/L, launched wave flux
  (L_w/L = 2.6e-5), and wave flux vs depth under quasi-adiabatic damping, for all ℓ and for ℓ ≤ 3;
  (f) envelope propagation diagram: N/2π, S_1/2π, ν_c, and the frequency below which waves are
  non-adiabatic (q = 1) for ℓ = 1 and 5; (g) cumulative τ_rad (ℓ = 1) from the FeCZ top, at 0.3/1/3 d⁻¹;
  (h) q at the same frequencies.
- **Propagation diagram** (x = r/R, whole star): N, S_1, S_2, the ℓ = 1 g-mode cavity, the convective
  core, the FeCZ, and the turnover frequencies (FeCZ 0.51 d⁻¹; core 0.0039 d⁻¹, i.e. t_c ≈ 41 d).
- What the figure shows at a glance: eddy-scale waves (ℓ ~ 100+) die within the first cells above the
  FeCZ; only ℓ ≲ 3 carry flux outward; above log T ≈ 5.0 the ν ≲ 1 d⁻¹ waves are non-adiabatic (q > 1),
  which is where the quasi-adiabatic τ piles up (panels e, g, h).
- `damping.py` gained `tau_cumulative` (vectorised τ vs depth) and `nu_nonad` (ω³ = K k_h² N² curve).

## 6. Next: non-adiabatic propagation with GYRE (decided 2026-09-24)

**Why GYRE:** the outer ~1 H_P is non-adiabatic (q ≫ 1) at the observed frequencies, and it is also where
the wave energy turns into emergent-flux fluctuations. A non-adiabatic oscillation code handles both at once
and gives δL/L at the photosphere directly (replacing the fitted 0.14 and the ad hoc damping).

**Availability:** GYRE 8.1 ships with MESA r26 (`$MESA_DIR/gyre/gyre-8.1.tar.gz`), but only as libraries;
the standalone `gyre` executable is not built, and there is no Lmod module. Build it from the tarball in a
separate directory (e.g. `~/software/gyre-8.1`, with the MESA SDK loaded), not inside `$MESA_DIR/gyre`.

**Input data (missing):** GYRE needs full structure files; MESA must write them with
`write_pulse_data_with_profile = .true.`, `pulse_data_format = 'GYRE'`, and — because the transition layer
sits just below the photosphere — `add_atmosphere_to_pulse_data = .true.`. Neither the v2 grid nor the T1
tests write them, and the v2 grid (running since 2026-09-24 ~20:40) should not be changed mid-run.
Plan: restart a few v2 models locally from their photos (or rerun T1a) with pulse-data output at the
triggered-profile points.

**Physics route (to decide):** GYRE solves free oscillations (and tidally forced ones), not a travelling wave
launched at an inner boundary. Two options:
1. **Stochastically excited, strongly damped non-adiabatic g modes of the envelope cavity.** Compute
   non-adiabatic modes (complex ω) over the relevant range; take the damping rates η and the surface
   δL/L per unit mode energy; set mode energies from the FeCZ excitation (power ∝ source spectrum at ω,
   overlap with the FeCZ top); sum Lorentzians. Where linewidths exceed mode spacings, the sum becomes the
   travelling-wave continuum. This is the standard stochastic-excitation formalism and fits GYRE's
   capabilities directly.
2. **Transfer function at real frequency:** use GYRE's non-adiabatic equations (`nad`) as a boundary-value
   problem at fixed real ω, with the wave injected at the FeCZ top and an outgoing/radiative surface
   condition. Closer to the travelling-wave picture, but needs custom driving (library use from Fortran
   rather than the stock executable).
Start with option 1 on one profile (20 Msun, X_c ≈ 0.5), and check that the damping rates reproduce the
quasi-adiabatic τ where q ≪ 1.
