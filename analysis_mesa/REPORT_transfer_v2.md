# Transfer-model tests on the v2 MW grid: does the FeCZ wave model hold?

*2026-09-24, ~23:40. Grid v2 (template_v2, Asplund09 mixture and opacities, rotation fixes; see
`../models/GRID_LOG.md`), MW, ω/ω_c = 0, 0.2, 0.4, 0.6. Same code and observational samples as
`REPORT_transfer.md` (v1). Outputs: `data_v2/`, `figures_v2/`. Scripts run with
`RN_GRID=.../models/grids_v2 RN_DATA=data_v2 RN_FIG=figures_v2`; side-by-side tables from
`python3 compare_v1_v2.py [DIR_A] [DIR_B]`.*

## Coverage caveat and the matched comparison
When this was run, v2 was incomplete: 5–21, 23, 24 and 30 Msun were done; M22 and M25 were missing (crash, rerun
pending); 40–120 Msun were mid-MS (X_c 0.16–0.69). v2 therefore covers fewer observed stars (v_macro
460 against 523; ν_char 163 against 173). To separate coverage from physics, `data_v1matched/` holds v1 tracks
restricted to v2's coverage (same masses, each cut at the X_c v2 has reached). **All v1-vs-v2 statements below
use that matched set.**

## Headline: S3t (F_w = 𝓜_t F_c, v saturated at c_s) no longer beats the empirical planes

MW, ω = 0, same stars:

| | v1 physics | v2 physics |
|---|---|---|
| v_macro ΔBIC (β = 1) | +4.5 | +19.2 |
| α₀ ΔBIC | +4.4 | +6.7 |
| ν_char ΔBIC (ν = 10^a/(2π t_c)) | **−22.4** | **−1.4** |
| **joint ΔBIC** | **−13.6** | **+24.5** |
| ratios L / τ / M (observed 5.3 / 6.4 [4.9, 8.6] / 2.2) | 4.8 / 7.0 / 1.8 | 5.7 / 9.2 / 2.0 |

(The unmatched v1 numbers in `REPORT_transfer.md`: joint −20.9.)

Free exponents (what changed):
- **ν_char vs 1/(2π t_c): β = 0.93 [0.79, 1.10] → 0.52 [0.40, 0.66]**; the normalization is unchanged
  (10^0.45 → 10^0.43). The v2 model turnover frequency varies about twice as much across the HRD as the
  observed ν_char. Over the OB region, d log ν_c/d log T_eff = 2.85 (model) against 1.47 (GP map of ν_char);
  the median offset log ν_char − log ν_c is 0.40 dex (v1: 0.49).
- v_macro vs the saturated wave velocity: β = 1.67 → 1.92.
- α₀ vs 𝓜_t F_c/F: β = 1.58 → 1.67 (essentially unchanged).
- Bare MLT v_c,max now fits v_macro better: β = 0.73 → 0.81; ΔBIC −2.6 → −15.7.

## Cause: the Asplund09 mixture, not the code changes
20 Msun, ω = 0 (T1a = template_v2 code with the v1 GS98 mixture):

| X_c = 0.5 | v_max [km/s] | t_c [s] | ν_c [d⁻¹] | F_c/F |
|---|---|---|---|---|
| v1 grid (v1 code, GS98) | 4.27 | 26,460 | 0.520 | 1.33e-3 |
| T1a (v2 code, GS98) | 4.37 | 27,174 | 0.506 | 1.34e-3 |
| v2 grid (v2 code, A09) | 6.01 | 20,996 | 0.655 | 3.32e-3 |

The code changes (FeCZ finder, mlt_vc, new columns) move non-rotating quantities by ≤ 2–3%. The mixture change
(at fixed Z = 0.014, A09 is ~solar Fe while GS98 is sub-solar) raises F_c/F ×2.5, v ×1.4 and ν_c ×1.3 at
20 Msun. The effect varies across the HRD: the FeCZ onset moves down (first FeCZ at log ℒ_spec 2.58 against
2.73; v_c = 3 km/s at 3.08 against 3.33; observed onset 3.16).

## What still holds
- Binned medians of v_macro, α₀ and ν_char are tracked by the FeCZ-based predictions (`figures_v2/fig_transfer.png`).
- The onset is still reproduced by a v_c ≈ 1–3 km/s threshold.
- α₀ ∝ 𝓜_t F_c/F is stable between v1 and v2 (same exponent, normalization and ΔBIC to within ~2).
- Core-excited IGWs still fail (ratios ~1.6 / 2.0 / 1.3; α₀ ΔBIC +80 at ω = 0).
- The ratios along ℒ and M remain close to observed; the τ ratio (9.2) is now just above the observed range.

## Rotation (now usable: the FeCZ is resolved at all ω)
Resolved FeCZ in rotating v2 models: median thickness 0.016–0.023 R, <1% thin (≤0.002 R), step jitter
0.0036 dex, all the same as ω = 0. In v1 at ω = 0.2/0.4: 0.0004/0.0006 R, 68/67% thin, 0.10/0.08 dex.

Joint ΔBIC in v2, by scenario:

| ω | S0 bare MLT | S3t waves | S4 core IGW |
|---|---|---|---|
| 0.0 | +6.5 | +24.5 | +58.5 |
| 0.2 | −3.6 | +19.1 | +31.7 |
| 0.4 | −13.7 | +12.1 | +4.3 |
| 0.6 | **−31.1** | +19.5 | −15.8 |

With v2 physics the plain MLT scalings (v_c,max; F_c/F; 1/(2π t_c)) do best, and better with more rotation.
Rotation is a population property, however (the observed stars have a range of v sin i), so a single-ω grid is
not the right comparison; a population-weighted mix of ω should be tested.

## Verdict
The FeCZ connection survives: onset, binned trends, and α₀ scaling hold, and core IGWs still fail. The specific
S3t transfer model does not survive the move to a consistent (A09) composition. Its success in v1 rested mainly
on the ν_char–turnover relation having a slope of 1, and with A09 that slope drops to 0.5. **The results are
sensitive to the Fe content through the Fe opacity bump**, so composition is a first-order systematic, not a
detail. The frequency physics is also exactly where the GYRE work (WAVE_MODEL.md) says the simple turnover
picture is incomplete.

## Next
1. Finish the v2 grid (M22/M25 rerun; massive tracks), then rerun to confirm with full coverage.
2. Composition bracket: a small v2 grid with GS98 vs A09 at fixed Z, or A09 at Z = 0.010 and 0.018, to see how
   strongly the conclusions depend on Fe. Galactic B stars are close to the Nieva & Przybilla (2012) cosmic
   abundance standard (≈ solar A09), which favours v2.
3. Population-weighted rotation (mix of ω by the observed v sin i distribution).
4. Frequencies: use GYRE (damped-mode continuum) rather than 1/(2π t_c).
5. Revise §5.6 of the paper accordingly (it currently quotes the v1 S3t result).
