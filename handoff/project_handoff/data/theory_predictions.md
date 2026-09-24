# Theoretical Predictions for Red-Noise (SLF) Scaling Relations

**Purpose.** Compile the predicted scalings of the three red-noise parameters —
α₀ (amplitude at zero frequency), ν_char (characteristic frequency), and γ
(high-frequency logarithmic slope) — with macroscopic stellar parameters under
the two competing excitation mechanisms, and identify the **diagnostic**
observables where the two scenarios diverge.

The red-noise / stochastic low-frequency (SLF) variability is universally
parameterised by a (super-)Lorentzian:

    α(ν) = α₀ / [ 1 + (ν/ν_char)^γ ] + C_w

with ν_char = 1/(2πτ) the inverse of a characteristic timescale τ.

---

## Mechanism A — Internal Gravity Waves (IGWs) excited by the convective core

Sources: Rogers et al. (2013); Aerts & Rogers (2015); Edelmann et al. (2019);
Ratnasingam et al. (2020); **Anders et al. (2023, Nature Astronomy)** —
the most quantitative modern prediction (3D Dedalus wave-generation +
wave-propagation simulations at realistic stellar luminosities, transfer
function via GYRE).

Physical picture: turbulent core convection stochastically excites a spectrum
of gravity waves that propagate through the radiative envelope and produce
surface brightness fluctuations. The wave luminosity spectrum is set by the
core convective turnover frequency (~0.1 d⁻¹ for a 15 M⊙ ZAMS star).

**Predicted scalings (Anders et al. 2023):**

| Parameter | Dependence | Notes |
|---|---|---|
| α₀ | **increases** with L and M | but predicted amplitude ≈ 0.06 µmag (15 M⊙), 3–4 orders of magnitude **below** the observed ≳10 µmag |
| ν_char | **decreases** with increasing L and M | signal shifts to lower ν and higher amplitude as mass rises |
| γ | high-frequency tail set by convective excitation efficiency + cancellation of high-ℓ modes | resonant g-mode peaks superposed |
| rotation | non-rigid (differential) rotation shifts power to **lower** frequency (Rogers D11 case); rigid rotation → weak effect | Anders: a P_rot=10 d run raises α₀ only to 0.21 µmag; **no** predicted correlation of surface amplitude with v sin i |
| Z | wave *generation* nearly Z-independent (set by core convection, interior) | envelope propagation only weakly Z-dependent |

**Key qualitative signatures of IGWs:**
1. ν_char **anti-correlates** with L (and mass): more luminous ⇒ lower ν_char.
2. Amplitude is essentially **independent of rotation** and of surface
   properties (Teff), because the driver is deep in the interior.
3. Because generation is interior, the mechanism is largely **metallicity-independent**
   and should persist unchanged from Milky Way → LMC → SMC.
4. Anders et al.'s central result: the predicted **amplitude is far too small**
   (0.06–0.2 µmag vs observed ≳10 µmag), so pure core-excited IGWs
   **cannot** be the dominant source of the ubiquitous red noise — unless
   amplitude-boosting physics missing from the simulations (strong rotation,
   evolved structure) intervenes.

---

## Mechanism B — Sub-surface (Fe/He opacity-bump) convection zones

Sources: Cantiello et al. (2009); **Cantiello et al. (2021, ApJ)** — the most
quantitative modern prediction (1D MESA models, MLT properties of the iron-bump
convection zone, FeCZ); Schultz et al. (2022); Jermyn et al. (2022).

Physical picture: the opacity peak from iron ionisation at T ≈ 150 kK drives a
thin, inefficient near-surface convection zone (FeCZ; a HeCZ at lower L). Its
turbulent motions perturb the surface directly (and/or excite locally-driven
waves), producing both SLF photometric variability and spectroscopic
macroturbulence.

**Predicted scalings (Cantiello et al. 2021):**

- ν_char = 1/(2πτ_c), with τ_c = α_MLT·H_P / v_c the **convective turnover
  time** in the FeCZ. Typical τ_c ≈ 0.1–2 d ⇒ ν_char ≈ 0.1–2 d⁻¹, matching the
  observed range.
- α₀ ∝ FeCZ convective flux ratio F_c/F_* (∝ ρ_c v_c³), i.e. the vigour of the
  sub-surface convection.
- All FeCZ properties depend almost exclusively on **position in the
  spectroscopic HR diagram**, ℒ ≡ T_eff⁴/g (Langer & Kudritzki 2014).

| Parameter | Dependence | Notes |
|---|---|---|
| α₀ | **increases** with L and Teff (toward high ℒ) | traces F_c/F_* of the FeCZ |
| ν_char | **increases** toward higher L and higher Teff | via shorter τ_c; Bowman+2020 values ≈ factor 3 higher than model normalisation (MLT: ν_char ∝ α_MLT⁻¹) |
| γ | not strongly predicted; related to convective spectrum | — |
| rotation | second-order; rotation can modify FeCZ vigour | macroturbulence correlates with FeCZ velocity |
| Z | **strong** dependence: FeCZ vigour ∝ opacity ∝ metallicity. FeCZ is **absent below log(L/L⊙) ≈ 2.5** at Galactic Z, and **weakens/disappears at SMC Z** | the metallicity lever is the sharpest sub-surface diagnostic |

**Key qualitative signatures of sub-surface convection:**
1. ν_char **increases** with L and Teff — **opposite sign** to the IGW prediction.
2. α₀ and ν_char organise cleanly on the **spectroscopic HRD** (ℒ = T_eff⁴/g).
3. Amplitude correlates with **macroturbulence** (same driver).
4. **Metallicity dependence:** signal should **weaken at low Z** (SMC) and
   vanish for stars below the FeCZ stability threshold.

---

## The diagnostic observables (where the scenarios diverge)

| Observable test | IGW (core) | Sub-surface (FeCZ) | Discriminating power |
|---|---|---|---|
| **sign of dν_char/dL** | **negative** (ν_char falls with L) | **positive** (ν_char rises with L, Teff) | **HIGH** — opposite signs |
| ν_char–Teff | weak/none | **increases** with Teff | HIGH |
| absolute α₀ | ≈0.06–0.2 µmag (too small) | ≳ observed (right order of magnitude) | HIGH (amplitude budget) |
| α₀ vs rotation (v sin i) | **no** correlation | weak/second-order | MEDIUM |
| **metallicity (MW→LMC→SMC)** | signal **unchanged** | signal **weakens** at low Z; absent below FeCZ threshold | HIGH — the cleanest lever |
| correlation with macroturbulence | not required | **expected** (shared driver) | MEDIUM |
| γ vs evolutionary stage | resonant-peak structure, mass-dependent | mild | LOW-MEDIUM |

**Primary observables to test in the data:**
1. **The sign of the ν_char–ℒ (and ν_char–Teff) relation.** This is the single
   most powerful discriminator: IGWs predict ν_char decreasing with L, FeCZ
   predicts it increasing.
2. **The organisation of α₀ and ν_char on the spectroscopic HRD** (ℒ=T_eff⁴/g),
   which the FeCZ picture predicts should be smooth.
3. **The metallicity dependence** (Galactic vs LMC vs SMC ν_char/amplitude),
   which cleanly separates an interior (Z-independent) from an envelope
   (Z-dependent) driver.
4. **Any α₀–v sin i correlation**, which IGW theory says should be absent.

## What has already been measured

- **Bowman et al. (2019a, CoRoT; 2019b/2020, TESS+spectroscopy):** established
  the ubiquity of SLF variability across the OB main sequence; found ν_char and
  α₀ correlate with position on the spectroscopic HRD (ℒ) and with
  main-sequence age. Sample sizes: ~30 CoRoT OBAF (Paper I), 70 Galactic OB
  (Paper II, N=75 with full fits).
- **Cantiello et al. (2021):** showed the FeCZ ν_char and α₀ contours reproduce
  the observed HRD trends (within a factor ~3 in normalisation).
- **Anders et al. (2023):** showed core-IGW amplitudes are far too small ⇒
  argues against core IGWs as the dominant driver.
- **Bowman & Dorn-Wallenstein (2022):** GP-regression (celerite2) method for
  robust ν_char; 30 Galactic stars.
- **Bowman et al. (2024, SMC/LMC):** extended to low Z (23 SMC + 21 LMC); found
  SLF variability **persists** at low Z with the **same** mass/age trends,
  challenging a pure sub-surface picture for young stars; proposed a
  **transition**: core-excited IGWs dominate for younger/low-Z stars lacking a
  FeCZ, sub-surface convection dominates for more evolved stars.
- **Ma et al. (2024, LMC blue supergiants):** 20 BSGs show a universal modified-
  Lorentzian signal (ν_peak≈0.2 d⁻¹, γ≈2–5); favour the FeCZ but find no strong
  correlation of fit parameters with spectroscopic properties in their narrow
  sample.
