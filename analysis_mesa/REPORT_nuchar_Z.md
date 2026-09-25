# GP ↔ Lorentzian ν_char, and the metallicity test it unblocks

*2026-09-25. Code: `digitize_bdw22.py` → `data_obs/bdw22_nuchar.csv`; `nuchar_Z_test.py` → `data_obs/nuchar_Z_test.csv`
(`MODEL_SHIFT=1` for the model prediction). Data: Bowman & Dorn-Wallenstein 2022 (A&A 668, A134), arXiv source in
`handoff/data_survey_raw/bowman_dw2022/`.*

## 1. Cross-calibration (30 Galactic O stars fitted both ways)
The paper has no per-star table. Each appendix panel marks ν_char from both methods as vertical lines on a log axis. We read
them at the pixel level, calibrating the axis per image from the tick marks (227.5 px per decade, 1 px = 0.004 dex).
- **Validation:** the Bowman 2020 line reproduces the tabulated B20 ν_char (J/A+A/640/A36) for all 30 stars, with median −0.001
  and max 0.003 dex.
- **GP vs Lorentzian:** median log(GP/Lor) = **−0.02 dex**, scatter 0.19 dex, Spearman 0.73; log GP = −0.01 + 1.18 log Lor.
  The differences depend on the star (the Q/ν_damp subgroups of BDW22). The low-Q, ν ≈ 1 d⁻¹ stars have GP ≈ 0.8× Lorentzian;
  several 2–3 d⁻¹ stars have GP ≈ 2× Lorentzian. **No large global offset.**

## 2. Magellanic vs Galactic ν_char at matched sHRD position
Galactic Lorentzian plane (Bowman 2020 + Shen 2024, N = 196), log ν = c0 + c1 log ℒ + c2 log T_eff. Bowman+2024 GP stars inside
the Galactic domain: 35 (LMC 15, SMC 20). The 95% interval resamples both the calibration stars and the Magellanic stars.

| GP → Lorentzian mapping | LMC+SMC median residual [95%] | LMC | SMC |
|---|---|---|---|
| none | −0.12 [−0.22, −0.06] | −0.15 | −0.10 |
| median offset (+0.02) | −0.10 [−0.29, −0.03] | −0.12 | −0.08 |
| inverse linear fit | −0.12 [−0.28, −0.04] | −0.12 | −0.12 |
| running median (10 nearest in ν) | −0.03 [−0.14, +0.05] | −0.02 | −0.04 |

- The paper's raw contrast (1.83 vs 1.07 d⁻¹ in 4.4 < log T_eff < 4.65, −0.23 dex) is roughly **half explained by the
  Magellanic stars sitting at different (ℒ, T_eff)**. What remains is −0.03 to −0.12 dex, depending on how GP values are mapped.
- **Model prediction (v2, ω = 0):** FeCZ turnover frequency log ν_c(LMC) − log ν_c(MW) = **−0.16 dex** at the Magellanic star
  positions (16–84%: −0.50 to −0.03). There is no SMC grid yet.
- Reading: the observed shift is in the predicted direction and of comparable size for most mappings, but zero is not excluded
  under the mapping that best follows the ν-dependent method difference. The SMC is not below the LMC (−0.08 vs −0.12 under the
  median offset), whereas the FeCZ picture expects the SMC to sit lower. The test has now moved from "blocked by method" to
  "limited by N and by the star-dependent method scatter (0.19 dex)". It becomes decisive only with a homogeneous refit (the Van
  Daele SMC light curves; next).

## Paper
- §4.6 and limitation (1): replace "the two methods are known to return systematically different ν_char … no cross-galaxy
  comparison" with the BDW22 calibration (median −0.02 dex, 0.19 dex star-to-star) and the matched-position result above.
- Abstract, last sentence: the metallicity lever is no longer blocked by a method offset; it is limited by sample size.
