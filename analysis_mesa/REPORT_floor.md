# A noise-floor model for the BLOeM SMC PSF light curves: how far can ν_char be pushed?

*2026-09-25. Code: `smc_floor_model.py` (steps calib / fits / inject), `floor_analysis.py` (also `expected_R()`,
`expected_detection()`). Outputs in `data_obs/floor_*.csv`. Figure: `Vandaele_correspondence/fig3_floor_model.png`.*

## Method
- **Likelihood:** Whittle likelihood on the independent frequencies (step 1/T, 0.1–40 d⁻¹) of the periodogram power. This
  replaces the linear-vs-log least-squares ambiguity (a 0.17 dex effect on ν_char, see REPORT_vandaele.md).
- **Model power:** S = [α₀/(1+(ν/ν_c)^γ)]² + [k·F_s]² + C_w².
- **Floor F_s:** an empirical semi-Lorentzian, calibrated on the 104 sectors of the 38 stars Van Daele class as non-variable.
  - Fractional amplitude = sector offset + 0.09·(Tmag − 13), with scatter **0.47 dex**. The floor is nearly proportional to flux.
  - Shape (median ν_f, γ_f): sector 27 has 0.82 d⁻¹, 2.3; sectors 28, 67 and 68 have 0.23–0.35 d⁻¹, 1.2–1.4.
  - k has a lognormal prior with σ = 0.47 dex. The calibration leaves the fitted star out.
- **Detection:** SLF above the floor means ΔBIC(floor → SLF + floor) > 10.
- **Injection–recovery:** random-phase SLF (Rayleigh amplitudes following the semi-Lorentzian, γ = 2) is injected into all 104
  non-variable sectors, on a grid of ν_in = 0.3, 0.7, 1.5, 3 d⁻¹ and R = α₀,in / floor amplitude = 0.1–16 (3328 injections).
- **Galactic comparison:** the same Whittle fits (no floor needed) on our 532 SPOC 2-min sectors.

## Results
1. **Almost all SMC "SLF" is indistinguishable from the floor.** Against white noise alone, SLF is detected in 100% of sectors
   of *every* class. Against the floor, it is detected in 24% of SLF-class sectors, 27% of binary/rotational, 50% of pulsators
   and **17% of the non-variable class** (the false-positive rate). The naive median log ν_char is −0.70 for SLF-class and −0.62
   for non-variable stars.
2. **Recovery map (injections):**
   - SLF at 1.5–3 d⁻¹ is detected in ≥ 69% of cases once R ≥ 1, and in ≥ 97% for R ≥ 4. The floor-aware ν_char is then unbiased
     (|bias| ≤ 0.05 dex, MAD ≤ 0.07).
   - SLF at 0.3–0.7 d⁻¹ sits on top of the floor. Detection needs R ≳ 4–8, and the naive fit is biased by up to −0.6 dex.
   - For weak SLF (R ≤ 0.25) the naive fit is pulled towards the floor: −0.3 to −0.8 dex at ν_in = 1.5–3.
3. **Where Galactic-like SLF would sit.** The Galactic fractional α₀ plane (our Whittle SPOC fits, N = 159, rms 0.42 dex)
   evaluated at the 85 SMC star-sectors in the Galactic domain gives median log α₀/F = −2.60. The SMC floor is −2.41. So
   **R_exp ≈ 0.9 (median), 16–84% [0.18, 3.8]**, and only 48% of sectors have R_exp ≥ 1. The SMC PSF light curves sit right at
   the sensitivity needed.
4. **Forward model** (fold each sector's expected (ν, R) through the injection map):

| SMC SLF relative to Galactic | expected detection above floor | expected naive median log ν |
|---|---|---|
| same amplitude, same ν_char | 0.33 | −0.35 |
| amplitude ×0.5 | 0.24 | −0.38 |
| amplitude ×0.1 | 0.14 | −0.49 |
| ν_char ×0.63 (−0.2 dex) | 0.28 | −0.47 |
| ν_char ×0.40 (−0.4 dex) | 0.27 | −0.64 |
| ν_char ×0.25 (−0.6 dex) | 0.27 | −0.84 |
| **observed (85 sectors)** | **0.25** | **−0.69** |
| non-variable stars (floor alone, naive) | 0.17 | −0.62 |

## Reading
- The SMC light curves in the Galactic sHRD domain look like the floor. **Galactic-like SLF (Galactic amplitude and Galactic
  ν_char at the same ℒ, T_eff) is disfavoured**: it would give a naive median log ν ≈ −0.35 against the −0.69 observed, and a
  somewhat higher detection rate. Consistent alternatives are SLF weaker than about 0.1× Galactic at Galactic ν_char, or
  ν_char lower by ≳ 0.4 dex at comparable amplitude, or a mix. **Both are the direction the FeCZ picture predicts at low Z** (a
  weaker, slower Fe zone; LMC grid Δlog ν_c = −0.16, SMC grid pending). Neither supports "similar morphology to Galactic stars".
- **Is the floor instrumental or astrophysical?** This matters, but the conclusion above survives either way:
  - Instrumental (strong sector dependence: sector 27 has a 0.5 dex lower amplitude and a 3× higher ν_f; Van Daele rejected
    139 of 281 sectors for low-frequency excess): then SMC SLF is weak or slow, as above.
  - Astrophysical (the "non-variable" stars have SLF too, and F ∝ flux): then the floor *is* the SMC SLF, with log ν ≈ −0.5
    to −0.6 at log ℒ 3.5–4.2. That is again 0.4–0.5 dex below the Galactic plane.
  - Either way the SMC stars do not show Galactic-like ν_char.
- **Caveats:**
  - The floor is calibrated on 38 stars, with 0.47 dex scatter and sector-dependent shape.
  - Only 54 of 91 SMC stars have spectroscopic (T_eff, log g), so 30 stars / 85 sectors are in the domain.
  - The Galactic comparison uses SPOC PDCSAP (crowding-corrected) photometry of bright stars; the SMC uses PSF photometry.
  - The forward model assumes the Galactic plane in α₀/F and ν transfers star by star, and uses γ = 2 injections only.
  - A decisive version needs the floor measured independently (background/neighbour sources in the same FFIs, or PSF light
    curves of Galactic OB stars made the same way) and the SMC v2 grid for the prediction.
