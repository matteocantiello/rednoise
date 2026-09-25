# Pedersen & Bildsten 2025: calibration, and the red-noise amplitude onset revisited

*2026-09-25. Code: `ingest_pb25.py` → `data_obs/pb25_calibration.csv`, `data_obs/pb25_stars.csv`; `onset_pb25.py` →
`data_obs/onset_pb25.csv`. Data: Pedersen & Bildsten 2025, MNRAS 539, 2742 (Zenodo 15261328), tables A1–B2.*

## Calibration onto our scale (43 Bowman 2020 stars refitted by PB25 with valid errors)
| relation | slope | intercept | scatter | Spearman |
|---|---|---|---|---|
| log α₀ [µmag] (Bowman) vs log PSD α₀ [ppm²/µHz] (PB25) | **0.505** | −0.296 | 0.32 dex | 0.93 |
| log α₀ [µmag] vs log rms [ppm] (model-independent) | 1.40 | −2.05 | 0.33 dex | 0.91 |
| log ν_char (Bowman) − log ν_char (PB25, µHz → d⁻¹) | — | +0.151 | 0.32 dex | 0.88 |
| γ (Bowman) − γ (PB25) | — | −0.80 | — | — |

The amplitude slope of 0.5 is what an amplitude ∝ √(power density) conversion gives. The 0.32 dex scatter is larger than the
0.19 dex Bowman–Shen study-to-study systematic.

## The amplitude onset
Hinge fit with a log T_eff covariate (as in the paper; `lmc_vmac_test.hinge_boot`):

| sample | N | stars below hinge | hinge | ΔBIC |
|---|---|---|---|---|
| paper set: a0_ok, log T_eff ≥ 4.0 (MW µmag + 26 LMC ppm) | 222 | 19 | 3.16 [2.96, 3.34] | **+10.9** (paper: +16.3) |
| same, **without the 26 LMC ppm stars** | 196 | 19 | 3.16 | **+0.9** |
| primary + 49 CygOB (calibrated) | 245 | 15 | 2.82 | +38.4 |
| same, log ℒ ≥ 2.0 (drops one 2.5 Msun, 10.7 kK star with α₀ ≈ 1.8e4 µmag) | 244 | — | — | **−2.7** |
| same, also without HD 37711 and HD 27563 (SPB-like spectra) | 242 | — | — | −3.3 |
| PB25 internal (B20 refit + CygOB), log rms, log ℒ ≥ 2.0 | 91 | — | — | −1.6 |
| PB25 internal, log rms, log ℒ ≥ 2.0, no SPB-like | 89 | — | — | −6.0 |

- **The red-noise amplitude onset in the paper is not robust.** It needs the 26 LMC stars whose ppm amplitudes (Bowman+2019b, Ma+2024)
  sit at high ℒ (log ℒ 3.3–4.4) and are 0.34 dex above Galactic stars at the same ℒ. Those stars steepen the bright end and create
  the hinge. Galactic µmag amplitudes alone, the homogeneous PB25 set (with its model-independent rms), and the combination all
  prefer a **single power law in ℒ over log ℒ ≈ 2–4.3** (slope 1.4–1.8 with the T_eff covariate).
- A strong break appears only because of one low-mass star at log ℒ 1.31 that is not a massive-star SLF object; the two
  SPB-like B20 stars raise the low-ℒ end further. Neither is a real floor.
- **The macroturbulence onset stands** (log ℒ 2.96 [2.90, 3.10], ΔBIC +81 on 890 Galactic stars; REPORT_lmc_vmac.md).

## Consequences for the paper
- Abstract and §4.2: "Both v_macro (hot stars) and α₀ switch on at log ℒ ≈ 3.0–3.2" should become: v_macro switches on; the red-noise
  amplitude of Galactic stars rises as a single power law without a detectable break in homogeneous data (the apparent α₀ hinge
  comes from mixing the LMC ppm sample). The "α₀ ~ 30 µmag floor" statement (§5.3) is not supported.
- The discussion point "both signals switch on" (§5.1 *Onset*) holds for v_macro only. For the FeCZ interpretation this is not
  fatal: the model F_c/F and 𝓜F_c/F rise smoothly from the first FeCZ appearance (log ℒ ≈ 2.6 in v2), so a smooth power law in α₀
  with no sharp threshold is also what the models give; the v_macro onset is the sharper, threshold-like feature.
- Table tab_onset: add the rows above (at least "without LMC ppm stars" and "with CygOB, log ℒ ≥ 2").

## Also available now (not yet used)
- 49 CygOB stars with calibrated α₀ (µmag) and ν_char (d⁻¹), T_eff and ℒ, for the scaling fits and the transfer-model tests
  (`data_obs/pb25_stars.csv`, sample `PB25_CygOB`; the calibration scatter is 0.32 dex, so they need that weight).
- AV 232 (SMC) is in PB25 table A2 (classical L only).
