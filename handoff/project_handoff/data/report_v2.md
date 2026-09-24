# Red noise and surface turbulence in massive stars — expanded meta-analysis

*486 OB stars with red-noise fits and 927 with homogeneous macroturbulence, assembled from
15 published tables. Update and extension of the first-pass (169-star) analysis.*

---

## Abstract

I expanded the aggregated red-noise catalogue for massive stars from 169 to **486 stars**
and built a companion catalogue of **927 stars with homogeneous IACOB macroturbulent
velocities**, then re-tested the internal-gravity-wave (IGW) versus sub-surface-convection
(FeCZ) scalings. The key addition is a **primary homogeneous Galactic sample of 196 stars**
(Bowman et al. 2020 ∪ Shen et al. 2024) with red noise in identical units (µmag, d⁻¹) and
spectroscopic luminosities ℒ = T_eff⁴/g, of which 144 also have macroturbulence, 126 have
masses and fractional main-sequence ages. Two independent studies that fit 24 of the same
stars agree to Δlog α₀ = +0.03 and Δlog ν_char = −0.06 dex (scatter 0.19 dex), validating
the merge. On the tripled sample the central diagnostic survives: **ν_char increases toward
hotter stars**, log ν_char = −0.23 log ℒ **+1.53 log T_eff** (all bootstrap-significant),
the sign opposite to the IGW prediction. The amplitude is governed by luminosity and
evolutionary state: log α₀ = **+1.56 log ℒ** − 3.16 log T_eff (R²=0.53). In a *separate*
regression on the 126 stars that have evolutionary parameters,
log α₀ = +1.91 τ + 1.83 log M (R²=0.45) — α₀ rises steeply with both fractional
main-sequence age and mass; α₀ also falls with gravity (Spearman ρ=−0.64). The new macroturbulence data confirm the FeCZ
shared-driver signature — α₀ ∝ v_macro^1.37 (ρ=+0.52, p=3×10⁻¹¹, N=144; the log–log slope
uses the 143 stars with v_macro>0) — and, importantly,
v_macro retains **independent** predictive power after ℒ and T_eff are removed
(partial ρ=+0.18, p=0.03; coefficient +0.28 [+0.15,+0.46]), though most of the raw
correlation is luminosity-mediated. Core IGWs remain 3–4 orders of magnitude too faint
(predicted 0.06–0.21 µmag vs observed median 350 µmag). A **clean division of labour**
emerges: α₀ is set by luminosity/mass/evolutionary state, ν_char by T_eff and rotation —
the two-parameter structure the FeCZ picture predicts. One genuinely unexplained result
is a significant residual rotation dependence of ν_char (partial ρ=+0.20, p=0.01) that
neither mechanism predicts.

---

## 1. New data acquired

Discovery was done by keyword search of the full VizieR table metadata via its TAP/ADQL
endpoint (searching descriptions for *macroturbulen\**, *IACOB*, *line-broadening*,
*stochastic low-frequency*, and OB-variability terms) rather than by guessing catalogue
identifiers. Fifteen tables were downloaded; ten carried new usable data.

### 1.1 Red-noise / SLF additions

| Source | VizieR | N | Units | New parameters |
|---|---|---|---|---|
| **Shen et al. 2024** (ApJS 275, 2) | `J/ApJS/275/2` | **150** | α₀ µmag, ν_char d⁻¹ | **spectroscopic ℒ, log g, M, R, τ, v sin i**; 298 per-sector fits (TESS S1–65) aggregated to 150 stars |
| **Bowman et al. 2019b** (NatAs 3, 760) | `J/other/NatAs/3.760/tabled13` | 114 | α₀ nominally ppm, ν_char d⁻¹ | K2 OB stars; SpType, Gaia BP−RP |
| Bowman et al. 2019b (LMC) | `.../tabled24` | 53 | as above | LMC TESS OB stars |
| Burssens et al. 2020 | `J/A+A/639/A81` | 98 | — (frequency analysis) | parent spectroscopy for Bowman 2020 |

Shen et al. is the single most valuable addition: it is **method- and unit-matched to
Bowman et al. 2020** and supplies the stellar parameters (mass, gravity, age fraction)
that the earlier compilation lacked. I verified that its `logLs` column reproduces
ℒ = T_eff⁴/g under my own normalisation to within 0.005 dex, so it merges directly onto
the existing luminosity axis.

### 1.2 Macroturbulence additions

| Source | VizieR | N with v_macro | Notes |
|---|---|---|---|
| Simón-Díaz et al. 2017 (IACOB III) | `J/A+A/597/A22` | 431 | dedicated macroturbulence study; also logT_eff, log L |
| de Burgos et al. 2024 | `J/A+A/687/A228` | 476 | B supergiants; T_eff, log g |
| Holgado et al. 2018 | `J/A+A/613/A65` | 128 | O stars; T_eff, log g |
| Simón-Díaz & Herrero 2014 (IACOB-broad) | `J/A+A/562/A135` | 199 | v sin i (FT) + v_macro (GOF) |
| Burssens et al. 2020 | `J/A+A/639/A81` | 72 | adopts IACOB values |

These four IACOB-family sources are **mutually consistent** (median pairwise offsets
≤1 km s⁻¹, MAD 3–8 km s⁻¹, ρ=0.70–0.88; Burssens vs Simón-Díaz 2017: ρ=0.97, offset 0.0),
so they were merged under a priority scheme into a single catalogue of **927 unique stars**
(`macroturbulence_catalog.csv`), median v_macro = 54 km s⁻¹.

### 1.3 A data-quality exclusion

Shen et al.'s own `vmacro` column was **excluded** from the macroturbulence catalogue.
Against the IACOB scale it is offset by −18 to −25 km s⁻¹ and barely correlated
(ρ=0.08–0.50); 32% of its values exceed 150 km s⁻¹ (maximum 440 km s⁻¹) versus 0.38% in
the IACOB family, and 32 stars have a `vmacro` value but no v sin i. This is the signature
of a degenerate total-line-broadening fit with rotation absorbed into the macroturbulent
component, not a radial–tangential macroturbulence. Shen's red noise and stellar
parameters — both independently validated — are retained; only this one column is dropped.

### 1.4 Independent validation of the merge

24 stars are fitted independently by Bowman et al. 2020 (Lorentzian amplitude-spectrum
fit) and Shen et al. 2024 (same functional form, TESS S1–65):

| Quantity | Median offset (Shen − Bowman) | Scatter (MAD) | Rank correlation |
|---|---|---|---|
| log α₀ | +0.031 dex | 0.186 dex | ρ=0.93 (p=8×10⁻¹¹) |
| log ν_char | −0.058 dex | 0.192 dex | ρ=0.69 (p=2×10⁻⁴) |
| γ | −0.076 | 0.391 | — |
| log ℒ | 0.000 | 0.000 | (same spectroscopy) |

The offsets are negligible and the correlations strong, so the two samples are merged.
**0.19 dex is adopted as the empirical systematic** on individual red-noise parameters.

### 1.5 Catalogue totals

| | before | after |
|---|---|---|
| stars with red-noise fits | 169 | **486** |
| stars with spectroscopic ℒ | 97 | **264** |
| stars with macroturbulence | 59 | **186** |
| primary homogeneous Galactic sample | 70 | **196** |
| separate macroturbulence catalogue | — | **927** |

---

## 2. Results on the expanded sample

### 2.1 Univariate (primary Galactic sample, N up to 196)

| Relation | ρ | p | N |
|---|---|---|---|
| **ν_char vs T_eff** | **+0.39** | 1.7×10⁻⁸ | 196 |
| ν_char vs ℒ | −0.02 | 0.79 | 196 |
| ν_char vs log g | +0.18 | 0.044 | 126 |
| ν_char vs v sin i | +0.24 | 3.0×10⁻³ | 152 |
| ν_char vs v_macro | +0.15 | 0.076 | 144 |
| **α₀ vs ℒ** | **+0.68** | 4.8×10⁻²⁸ | 196 |
| α₀ vs T_eff | −0.09 | 0.23 | 196 |
| **α₀ vs log g** | **−0.64** | 6.1×10⁻¹⁶ | 126 |
| **α₀ vs v_macro** | **+0.52** | 2.7×10⁻¹¹ | 144 |
| **α₀ vs τ** (MS age) | **+0.51** | 1.6×10⁻⁹ | 126 |
| **α₀ vs M** | **+0.51** | 8.1×10⁻¹⁰ | 126 |
| α₀ vs v sin i | +0.32 | 4.7×10⁻⁵ | 152 |
| γ vs anything | ≤0.17 | ≥0.10 | 196 |

### 2.2 Multivariate power laws (bootstrap median [16–84%])

ℒ, T_eff and log g are **algebraically degenerate**: log ℒ = 4 log T_eff − log g − const is
the *definition* of ℒ, so the three cannot be independent predictors and only the
{ℒ, T_eff} basis is used (including log g as well would give a rank-deficient fit).
The *published* values honour the identity closely but not perfectly: over the 126 stars
carrying all three, the residual has median 0.0005 dex and MAD 0.003 dex, and 116 of 126
lie within 0.01 dex. Ten stars exceed that, two of them by >0.1 dex (largest 0.90 dex) —
θ¹ Ori, HD 108, HD 191612 and HD 15137 among them, all known magnetic/peculiar or binary
objects whose published T_eff, log g and ℒ come from different sub-analyses and are
internally inconsistent. These residuals reflect catalogue inhomogeneity, not a failure of
the identity.

- **log ν_char = −0.23 [−0.32,−0.14] log ℒ + 1.53 [+1.20,+1.88] log T_eff**  (N=196, R²=0.09)
- with rotation: −0.26 log ℒ + 1.20 log T_eff **+0.34 [+0.25,+0.42] log v sin i** (N=152, R²=0.13)
- **log α₀ = +1.56 [+1.43,+1.68] log ℒ − 3.16 [−3.83,−2.55] log T_eff**  (N=196, R²=0.53)
- with turbulence: +1.43 log ℒ − 3.04 log T_eff **+0.28 [+0.15,+0.46] log v_macro** (N=144, R²=0.58)
- evolutionary (separate model, τ and mass available for the Shen subset only):
  **log α₀ = +1.91 [+1.62,+2.23] τ + 1.83 [+1.57,+2.10] log M**  (N=126, R²=0.45)
- ν_char shows **no** dependence on τ (−0.15 [−0.49,+0.20]) or mass (−0.18 [−0.49,+0.12]);
  the same model for ν_char has R²=0.005

These τ/mass coefficients are regression slopes from the `tau+log_mass` model and are **not**
the rank correlations of §2.1 (where α₀ vs τ and α₀ vs M give ρ=+0.51 each); all four rows are
in `multivariate_fits_v2.csv` under `model = tau+log_mass`.

All coefficients except the ν_char–τ and ν_char–mass terms are bootstrap-significant.

### 2.3 The one honest weakening

The ν_char–T_eff relation **persists in sign and slope** (+1.53 vs +1.39 on the original 70
stars) but its explanatory power collapses: R² = 0.09 on N=196 versus 0.46 on the curated
N=70. The Bowman 2020 sample was more homogeneous (one fitting pipeline, one spectroscopic
analysis, narrower selection). The larger sample is noisier, so **the diagnostic sign is
robust while the tightness of the original relation was partly a small-sample effect.**
This is the most important caveat added by the new data.

### 2.4 Macroturbulence: how much is independent?

The raw α₀–v_macro correlation is strong (ρ=+0.52, N=144) and the power law is
α₀ ∝ v_macro^1.37 (fitted on the 143 stars with v_macro>0; one star has v_macro=0 and is
excluded from the log fit only). But α₀ and v_macro both rise with luminosity, so most of the raw
correlation is luminosity-mediated. Removing ℒ and T_eff first:

- partial ρ(α₀, v_macro | ℒ, T_eff) = **+0.18, p=0.027** (N=144)
- partial coefficient +0.28 [+0.15, +0.46]

So the shared-driver signature is **real but modest** once luminosity is controlled — a
weaker statement than the raw correlation implies, and stated here deliberately.

---

## 3. What this does to the mechanism question

**The FeCZ picture gains on three fronts and the IGW picture on none.**

1. **The sign test survives tripling the sample.** ν_char rises with T_eff
   (+1.53 in log–log), as expected if ν_char = 1/(2πτ_c) traces a convective turnover time
   that shortens toward hotter, denser iron-bump convection zones. IGW theory predicts the
   opposite sign.
2. **A clean division of labour.** α₀ is set by luminosity, mass and evolutionary state
   (τ, log g) and is *independent* of T_eff at fixed ℒ; ν_char is set by T_eff and rotation
   and is *independent* of τ, mass and ℒ. This is precisely the two-parameter structure
   Cantiello et al. (2021) predict, where α₀ ∝ F_c/F_* and ν_char ∝ 1/τ_c are separate
   functions of position in the spectroscopic HRD.
3. **α₀ grows as stars evolve** (regression slope +1.91 per unit fractional MS age,
   N=126; rank correlation ρ=+0.51). As the
   envelope expands the FeCZ becomes more vigorous and sits closer to the surface. Core
   IGW amplitudes are not predicted to track envelope structure this way.
4. **Surface turbulence and photometric red noise share a driver** — now on 144 stars
   rather than 59 — with v_macro retaining independent predictive power after luminosity
   and temperature are removed. An interior mechanism has no reason to produce amplitude
   correlated with *surface* turbulent velocity at fixed L and T_eff.
5. **The amplitude budget is unchanged and still decisive**: 0.06–0.21 µmag predicted for
   core IGWs versus an observed median of 350 µmag (10–90%: 64–3607 µmag).

**Two things remain unexplained by either mechanism.** ν_char depends significantly on
rotation after ℒ and T_eff are removed (partial ρ=+0.20, p=0.013, N=152), and α₀ correlates
with v sin i (ρ=+0.32). Neither the FeCZ models (non-rotating 1D) nor the core-IGW
simulations predict this. Rotation is the most promising remaining discriminant, and it is
the one axis where both current theories are silent.

The metallicity conclusion from the first-pass analysis is unchanged: SLF variability
persists to SMC metallicity, which still argues for a metallicity-independent floor —
plausibly core-excited IGWs — in young, metal-poor stars that lack a vigorous FeCZ. The
new data are all Galactic (plus 53 LMC stars without spectroscopic parameters), so they do
not sharpen the metallicity lever.

---

## 3b. Property maps on the spectroscopic HR diagram

Following Cantiello et al. (2021), each property was mapped on the sHRD (ℒ = T_eff⁴/g vs
T_eff) with MIST v1.2 solar-metallicity tracks overlaid — first as coloured points, then as
a **reliability-weighted Gaussian-process posterior** rendered as filled contours with no
data points (Figures 4–5).

### 3b.1 HRD gradients (point maps)

A plane fit log(prop) = a·log ℒ + b·log T_eff + c gives:

| property | ∂/∂log ℒ | ∂/∂log T_eff | N | R² |
|---|---|---|---|---|
| v_macro | **+0.29** | +0.10 | 832 | 0.44 |
| α₀ | **+1.54** | −3.13 | 196 | 0.53 |
| ν_char | −0.23 | **+1.53** | 196 | 0.09 |
| γ (linear) | +0.15 | −2.10 | 196 | 0.10 |

v_macro and α₀ organise along the *same* HRD direction (upward, with ℒ); ν_char organises
*orthogonally* (leftward, with T_eff). α₀ climbs ~5× faster in log than v_macro per dex of
ℒ, so the photometric amplitude is far more sensitive to HRD position than the turbulent
velocity — consistent with α₀ tracking a flux *ratio* rather than a velocity, and explaining
why the direct α₀ ∝ v_macro^1.37 relation is shallower than the ratio of HRD gradients.

Building the v_macro map exposed a unit error in the macroturbulence catalogue: Holgado et
al. (2018) tabulate T_eff in kK, which had been logged as K for 62 stars. Corrected
(`macroturbulence_catalog.csv` v2); the red-noise catalogue and all correlation/fit results
were unaffected (no affected star entered any of them).

### 3b.2 Reliability-weighted GP maps

Each property was modelled as a GP (Matérn-5/2, anisotropic length scales) over
(log T_eff, log ℒ) with per-star measurement variance σᵢ² entering the likelihood. Two
choices matter:

- **The published formal errors are unusable as weights.** Bowman et al. 2020's `e_alpha0`
  is in µmag (ReadMe-confirmed), giving relative errors of 10⁻⁵–10⁻² — formal MCMC errors.
  The 24 stars fitted independently by Bowman and Shen disagree by 0.19 dex. The measured
  systematic (0.186 dex α₀, 0.192 ν_char, 0.391 γ, ~20% v_macro) is used as the floor.
- **Rotation penalty.** σ_extra = 0.35·v sin i/(v sin i + 100 km s⁻¹) dex, added in
  quadrature, targets the v sin i–v_macro degeneracy: 0.26 dex at 300 km s⁻¹. It
  down-weights 73% of the v_macro sample appreciably but only 5% of the red-noise sample.
  Turning it off shifts the maps by a median of 0.02–0.11 dex (locally up to 1.6 dex in
  sparse regions) — it refines rather than drives the result.

A GP with long correlation lengths extrapolates *confidently*: the v_macro fit reported
posterior sd ≈ 0.035 dex over 91% of the diagram including empty corners. Posterior variance
measures model confidence, not data support, so the trustworthy region requires **both**
posterior sd < 0.6× the data scatter **and** ≥3 stars within Δlog T_eff = 0.06,
Δlog ℒ = 0.25. This cuts trustworthy area to 69% (v_macro), 34% (α₀), 31% (ν_char), 21% (γ).

Fitted correlation lengths are themselves informative: v_macro is genuinely smooth
(≈2.2–2.7 in standardised units); ν_char's sit near the lower bound (≈0.3), so its map is
patchy — a faithful picture of a noisy field (R²=0.09), not a hidden smooth trend; γ's
T_eff length scale hit the upper bound, i.e. the fit found no T_eff dependence at all.

### 3b.3 Extending the maps to the cool, evolved side

The OB samples end at log T_eff ≈ 4.1. Two routes placed 146 additional stars on the sHRD:

- **Tier A — spectroscopic log g (47 stars).** Parent analyses of the Magellanic red-noise
  samples located by VizieR metadata search: Bestenlehner et al. 2025 (XShootU SMC+LMC,
  `J/A+A/695/A198`), Vink et al. 2023 (ULLYSES, `J/A+A/675/A154`), Urbaneja et al. 2017
  (LMC supergiants, `J/AJ/154/102`; publishes flux-weighted gravity log g_F = log g −
  4 log(T_eff/10⁴ K), inverted exactly), McEvoy et al. 2015 (VFTS BSGs, `J/A+A/575/A70`).
  Name matching (with alias splitting) plus a 2″ coordinate pass gave 34/64 Bowman 2024 +
  Ma 2024 stars and 13/53 Bowman 2019b LMC stars. New T_eff agree with the red-noise
  papers to median 0.000 dex (MAD 0.005): the same parameter sources Bowman used.
- **Tier B — MIST-track mass (99 stars).** log ℒ = log L − log M exactly, so stars with
  only classical L were placed via a mass from the nearest MIST track point (masses 5–120,
  phases 0–6). Calibrated on 33 stars with both routes: median offset −0.09 dex (applied),
  scatter 0.19 dex (propagated into the weights via each property's HRD gradient). This
  added the remaining Magellanic stars and Dorn-Wallenstein et al. 2020's 76 cool
  supergiants (`J/ApJ/902/24`; 48 RSG + 28 YSG; same functional form, ν_char = 1/2πτ).

Additional systematics entered the weights for GP-regression vs Lorentzian ν_char (+0.20
dex for Bowman 2024) and cross-method samples (+0.15). Two α₀ sets were **excluded** from
the amplitude map: Bowman 2024 (a PSD-maximum proxy) and Dorn-Wallenstein 2020 (unit field
"10⁻³" with white-noise values 1.5–52, too large for a ppt amplitude — quantity ambiguous).

Result (342 stars, log T_eff 3.53–4.76):

| region | N | ν_char (d⁻¹) median [16–84%] | γ median |
|---|---|---|---|
| hot OB (log T_eff ≥ 4.3) | 235 | 1.47 [0.45, 3.80] | 1.91 |
| evolved OB / BSG (4.0–4.3, ℒ > 3.6) | 22 | **0.27** [0.20, 0.46] | 2.58 |
| YSG (3.75–4.0) | 20 | **0.30** [0.17, 0.50] | 1.63 |
| RSG (< 3.75) | 51 | 0.59 [0.33, 2.33] | 1.72 |

The ν_char–T_eff trend persists over the full range (ρ=+0.49, p=5×10⁻²², N=342) but
**flattens to a floor of ≈0.3 d⁻¹ for blue and yellow supergiants**: the drop from ≈1.5 to
≈0.3 d⁻¹ occurs across log T_eff ≈ 4.3→4.1 and then stays flat. If ν_char traces a
convective turnover time, the floor marks where the FeCZ gives way to the deeper He II/H
convection zones of cooler envelopes. **RSGs break the monotonic trend** (median higher
than YSGs, 16–84% range 0.33–2.33 d⁻¹, some τ as short as 0.03 d): a different regime
(granulation on a fully convective envelope) that should be read descriptively, not within
the FeCZ framework. γ is steepest in the evolved OB/BSG region (≈2.6 vs ≈1.6–1.9
elsewhere). The α₀ map's boxy feature at log T_eff ≈ 4.1–4.2, ℒ ≈ 3.6–3.9 is a sparse
cluster of Bowman 2019b LMC and Ma BSG stars — the coverage mask working as intended, not
structure.

---

## 4. Caveats

- **Sample heterogeneity dominates the ν_char scatter.** R² for ν_char fell from 0.46
  (N=70) to 0.09 (N=196). The sign is robust; the normalisation and tightness are not.
- **α₀ units across samples.** Five conventions now (µmag amplitude, nominal ppm, ppm²/µHz
  PSD, mmag²/d⁻¹ PSDmax proxy). All amplitude fits are restricted to the µmag-matched
  primary sample. For the Bowman 2019b tables the VizieR byte-by-byte description leaves
  the α₀ unit field blank (`---`), so the unit was not confirmed from the catalogue
  metadata and I did not consult the journal article's own text; I carry it as *nominally*
  ppm, which would be a ≤0.04 dex offset from µmag — immaterial in log space, but these
  stars are excluded from the amplitude fits regardless.
- **Bowman 2019b K2 stars mostly lack stellar parameters**: only 17 of 114 cross-matched to
  spectroscopy within 2″ (K2 fields are ecliptic, IACOB is northern-bright), so they
  contribute to distribution-level statistics but not the scaling fits.
- **Shen's `vmacro` excluded** (§1.3); its red noise and stellar parameters are retained.
- **Mass and τ are model-dependent** (from Shen's evolutionary-track fits), so the
  α₀–τ and α₀–M slopes inherit those tracks' systematics.
- **v_macro is not a clean turbulence measurement.** It is a line-broadening fit parameter
  degenerate with rotation; the IACOB sources handle this consistently, which is why only
  they were merged.
- Selection effects (bright-star bias), T_eff-scale offsets between pipelines, and
  coherent-pulsator contamination all persist from the first-pass analysis.

---

## 5. Next steps that would matter

1. **Uniform re-analysis** of all 486 light curves with one method (GP regression) — the
   single largest remaining systematic, and the reason the ν_char relation loosened.
2. **Spectroscopic log g for the LMC/SMC samples** to put every star on the ℒ axis and
   make the metallicity test genuinely matched.
3. **Rotating FeCZ and IGW models**, since rotation is now a significant, unexplained term
   in both ν_char and α₀ and is where the two theories are equally silent.
4. Cross-match the 927-star macroturbulence catalogue against TESS full-frame photometry
   to build a red-noise sample selected *on* surface turbulence rather than on brightness.

---

## Appendix: artifacts

- `master_rednoise_catalog_v2.csv` / `.parquet` — 486 stars, 8 samples, harmonised.
- `macroturbulence_catalog.csv` — 927 stars, homogeneous IACOB v_macro + v sin i.
- `correlations_summary_v2.csv`, `multivariate_fits_v2.csv`, `correlation_matrix_v2.csv`.
- Figures: `fig_v2_scalings.png`, `fig_v2_heatmap.png`, `fig_v2_theory.png`,
  `fig_hrd_properties.png` (point maps), `fig_hrd_interpolated_v2.png` (GP maps, extended).
- HRD maps: `hrd_gradients.csv`, `rednoise_sHRD_extended.csv` (342 stars with sHRD position,
  tier, systematics), `hrd_interpolated_grids_v2.csv.gz` / `.npz` (fitted fields on both grids),
  `hrd_predictions_tracks.csv` (expected values along MIST tracks, OB-only maps).
- First-pass analysis retained: `report.md`/`.pdf`, `theory_predictions.md`, and the
  original seven figures.

*New sources: Shen et al. (2024, ApJS 275, 2); Bowman et al. (2019b, Nat. Astron. 3, 760);
Burssens et al. (2020, A&A 639, A81); Simón-Díaz et al. (2017, A&A 597, A22);
Simón-Díaz & Herrero (2014, A&A 562, A135); Holgado et al. (2018, A&A 613, A65);
de Burgos et al. (2024, A&A 687, A228). All data public via VizieR/CDS.*
