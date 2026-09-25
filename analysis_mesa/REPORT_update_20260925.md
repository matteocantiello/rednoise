# Analysis update on the near-complete v2 grids, with all datasets (2026-09-25 evening)

*Grid state: 304/408 done at 17:00 (GRID_LOG); every model ≤ 25 Msun is complete except SMC w0.4 M24/M25, and 30–120 Msun
are still running. Extracts re-made for all 11 available sub-grids (`data_v2/mesa_{ms,post}_<Z>w<w>.csv`; previous ones in
`data_v2/snap_20260925_1100/`). Regime tables for MW, LMC and SMC (ω = 0). Log: `data_v2/extract_20260925_1715.log`,
`data_v2/chain_20260925_1720.log`. New figure script: `make_fig_grids_obs.py` → `figures_v2/fig_shrd_grids_obs`,
`fig_profiles_Z`, `fig_fecz_window`, and `data_v2/grids_obs_summary.csv`.*

## 1. What did not change
The models added since 11:00 are high-mass tracks past the main sequence, where there are no data, so every MW result is stable:
- **Scenarios** (joint ΔBIC, negative = better than the empirical planes): S5 −9.0 / −16.8 / −20.9 / −34.4 for ω = 0 / 0.2 /
  0.4 / 0.6 (previously −8.9 / −16.8 / −20.8 / −37.1). S0 +18.7 … −10.6. Core IGW (S4) +88 to +91 everywhere. S5's
  hot-star v_macro is still poor (+54 to +64).
- **MW regime crossings** are identical to 0.01 dex: v_c = 3 km/s at log ℒ 3.02; Γ_Fe = 0.8 at 3.65 and 1.0 at 3.88.
- **LMC v_macro shift** Δ = +0.57 (method offset fixed at 0) or +0.40 (free); **v_mic** has no onset, and the model v_c,max
  → v_mic test gives ΔBIC +42. Both unchanged.

## 2. New: the SMC grid
| threshold (hot MS, ω = 0) | MW | LMC | SMC |
|---|---|---|---|
| onset proxy v_c = 3 km/s | 3.02 | 3.62 | **4.08** |
| Γ_Fe = 0.4 | 3.11 | 3.34 | 3.56 |
| Γ_Fe = 0.8 (saturation proxy) | 3.65 | 3.89 | **4.08** |
| Γ_Fe = 1.0 | 3.88 | 4.10 | 4.28 |
| MLT Mach 0.2 | 3.75 | 4.07 | 4.26 |
| FeCZ present in 50% of models (from the maps) | ≈ 2.6 | ≈ 3.05 | ≈ 3.55 |

Going MW → LMC → SMC, the FeCZ appears 0.45 and 0.95 dex higher in log ℒ, and the vigorous-convection proxies shift by
0.6 and 1.06 dex. Saturation (Γ_Fe = 0.8) shifts by 0.24 and 0.43 dex.

## 3. Models against data, per galaxy (`grids_obs_summary.csv`, hot stars log T_eff > 4.3)
- **Scale factors, set on the MW only:** ν_char ≈ 3.39 × ν_c (the turnover frequency), and α₀ ≈ 2.4 × 10⁵ µmag × 𝓜_t F_c/F.
  These factors are then applied unchanged to LMC and SMC.
- **ν_char, observed − model:** MW −0.02 dex (by construction); **LMC +0.01** (26 stars, Bowman 2019b/2024 GP and Lorentzian,
  matching the MW-calibrated FeCZ prediction); **SMC +0.74** (17 Bowman 2024 stars with a model FeCZ; 5 more have none).
- **v_macro, observed − model v_c,max:** MW +0.35, LMC +0.25. Their difference, 0.10 dex, compares with the model's LMC−MW
  drop of ~0.3–0.5 dex. This is the LMC v_macro shift result again.
- **Profiles (`fig_profiles_Z`):**
  - v_mic tracks the MW model v_c,max at log ℒ 3.3–3.8 and flattens above that.
  - The model α₀ predictor drops by 1–2 dex from MW to LMC at fixed ℒ. There are no µmag-comparable Magellanic amplitudes to
    test this.
  - The model LMC ν_c is lower than MW at log ℒ 3.5–3.9; the observed LMC values follow the LMC model.

## 4. The Magellanic samples against the SMC FeCZ window (`fig_fecz_window`; correspondence items 1 and 2)
| sample (placed with its own-Z grid) | N on grid | where its own-Z models have no FeCZ |
|---|---|---|
| BLOeM spectroscopic sample (Bestenlehner+2025) | 757 | **62%** |
| Van Daele PSF light curves, all | 46 | 6.5% |
| – their SLF stars | 24 | **0%** |
| – their no-variability stars | 15 | 20% |
| Bowman+2024 SMC | 23 | 22% (5 stars) |
| LMC: Serebriakova v_macro | 115 | 1.7% |
| LMC: Bowman 2019b/2024 | 29 | 0% |

- **Selection bias, quantified.** 62% of the spectroscopic BLOeM stars sit where the SMC models have no FeCZ, but *none* of
  the 24 SLF detections in the PSF light curves do. The magnitude-limited light curves do not probe the window, so the Van
  Daele SMC sample cannot test FeCZ-driven SLF.
- **The Bowman & Van Daele 2024 window stars** (the actual basis of the "insensitive to Z" claim): 5 of 23 lie in the SMC
  no-FeCZ window, all at log ℒ 3.2–3.5, just below the SMC boundary at ≈ 3.55. They are AzV 468, 2dFS 3780, NGC 346 ELS 25,
  NGC 346 ELS 46 and [M2002] SMC 81469.
  - Their GP amplitude proxy is 0.9 dex lower than for the 9 stars just above the boundary (median log −4.33 vs −3.40), but
    this is not significant (MWU p = 0.52, N = 5 vs 9).
  - Their ν_char is higher (median 2.0 vs 1.0 d⁻¹).
  - Weak, high-frequency red noise is what a marginal detection looks like, which fits the user's recollection (e.g. TIC
    181887485), but N is too small to claim it.
- **A real tension for the FeCZ picture: ν_char near the SMC threshold.** The 9 Bowman 2024 stars just above the SMC
  boundary (log ℒ 3.5–3.8) have GP ν_char ≈ 1 d⁻¹. The SMC models there have weak FeCZs (v_c 0.1–0.6 km/s) and turnover
  frequencies ×3.4 ≈ 0.05–0.15 d⁻¹, which is 1 dex lower. This drives the SMC +0.74 dex offset. Well above the boundary
  (log ℒ > 3.8) the SMC ν_char (0.3–1 d⁻¹) and the SMC model agree within ≈ 0.3 dex. Possible readings:
  - in stars with a marginal FeCZ, the red noise is not set by the FeCZ turnover (another source, or the core);
  - the 1D MLT turnover is not the relevant frequency near threshold;
  - the GP ν_char of weak signals is biased high (compare the SMC PSF floor results, where weak SLF biased fitted ν *low*, so
    the direction is not obvious).
  This needs Bowman 2024's per-star detection significance, and a check with LMC stars near the LMC boundary (none sit there
  now).

## 5. Consequences
- For the paper and the Van Daele exchange, §4 gives the quantitative form of the selection-bias argument. It also gives an
  honest caveat: near the SMC threshold the observed frequencies are not what the FeCZ turnover predicts.
- The LMC is consistent with the FeCZ predictions in both ν_char (+0.01 dex) and the direction of the v_macro shift. The SMC
  test is limited to 23 GP stars, 14 of them within 0.3 dex of the FeCZ boundary.
- The MW-based conclusions (scenario ranking, onset, v_mic) stand as reported.
