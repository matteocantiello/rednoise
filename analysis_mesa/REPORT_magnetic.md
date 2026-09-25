# Magnetic stars: is the red noise suppressed?

*2026-09-25. Code: `magnetic_test.py` → `data_obs/magnetic_test.csv`, `data_obs/magnetic_stars.csv`. Data: Shen+2023
(J/ApJ/955/123; same SLF fit and µmag units as Shen 2024); B_d from Shultz+2022 (J/MNRAS/513/1429, A1), otherwise B_p from
Petit+2013 Table 1 (MNRAS 429, 398; parsed from the arXiv source into `handoff/data_survey_raw/petit2013/`, '>' = lower limit).*

## Setup
- Baseline: the primary sample (Bowman 2020 + Shen 2024, a0_ok), **after removing 10 known magnetic stars** that were in it
  (HD 108, HD 191612, HD 37742, HD 57682, HD 46328, HD 44743, HD 52089, …): 186 stars. Plane log y = c0 + c1 log ℒ + c2 log T_eff.
- Scale check: 9 stars fitted by both Shen 2023 and the primary sample agree in α₀ (median +0.02 dex, MAD 0.13); ν_char −0.15 dex.
- Magnetic sample: 30 Shen 2023 stars inside the baseline domain (log ℒ 2.6–4.3, log T_eff 4.14–4.64). All sit where the MW v2
  models have an FeCZ. 26 have a field value (16 measured, 10 lower limits).
- Model B_shutoff: the equipartition field of the FeCZ (history FeCZ_B_shutoff_conv, MW v2, ω = 0, MS), interpolated at each
  star's position: 1.6 kG at log ℒ 2.6, rising to 7.7 kG at 3.65.

## Results
| | log α₀ residual | log ν_char residual |
|---|---|---|
| all magnetic (N = 30), median | **+0.22** (Mann-Whitney p = 0.03 vs baseline) | −0.17 (p = 0.13) |
| rotationally modulated (N = 6) | +0.64 | −0.08 |
| other classes (N = 24) | +0.14 | −0.17 |
| Spearman vs log B_d, measured only (N = 16) | −0.29 (n.s.) | +0.07 |
| Spearman vs log B_d, incl. lower limits (N = 26) | +0.15 | −0.30 (p = 0.14) |
| B_d > B_shutoff, measured (N = 3: HD 37479, HD 184927, HD 64740) | −0.16 (values +1.21, −0.26, −0.61) | −0.18 |
| baseline rms | 0.45 | 0.40 |

1. **No suppression.** Magnetic OB stars have red-noise amplitudes equal to or *higher* than non-magnetic stars at the same
   (ℒ, T_eff). The excess comes from the rotationally modulated and multi-periodic/binary stars (e.g. HD 47129, HD 149277), which
   points to unremoved coherent or rotational power leaking into the low-frequency fit, not to extra convective driving.
2. **No dependence on field strength.** α₀ residuals show no trend with B_d, including across the ~50× range from HD 37742
   (0.06 kG) to NGC 1624-2 (> 20 kG).
3. **The critical test has little power.** Only 3 stars with a measured B_d exceed the model B_shutoff; 6 do when lower limits are
   counted. Two of the three measured cases lie below the plane (−0.26, −0.61) and one lies well above (σ Ori E, +1.21, a
   rotational variable). With a 0.45 dex baseline scatter, N = 3 cannot detect even a factor-of-3 suppression.

## Reading
- If the photometric SLF came from FeCZ convection that a surface field above B_shutoff shuts off (Cantiello & Braithwaite 2011), the strongest-field stars should drop below the plane. They do not as a group,
  though the sample able to test it is 3–6 stars.
- The result is equally awkward for a simple core-IGW picture in which strong near-surface fields damp the waves
  (a prediction discussed in the literature; reference to check before citing). What the data say is that **α₀ does not know about the surface field**. Either the photometric
  SLF is set below the magnetically dominated layers, or ~kG dipoles do not quench the FeCZ velocity field in the way the
  equipartition criterion assumes. Macroturbulence may be the sharper test (reported low v_macro in some strongly magnetic O stars;
  literature to verify). A v_macro/B_d compilation for magnetic B stars is the natural next step (MOBSTER spectroscopy), and it would
  separate the spectroscopic and photometric signals.
- Paper: one paragraph in §5 (tests) with the table above. It is not a discriminant at present sample size, but it rules out a
  strong (> 0.5 dex) blanket suppression of α₀ in magnetic stars.
