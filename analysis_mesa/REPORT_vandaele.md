# Van Daele+2026 (BLOeM SMC): can we reproduce their SLF fits?

*2026-09-25. Code: `vandaele_fit.py` → `data_obs/vandaele_fits.csv` (0.1–40 d⁻¹, e⁻/s; the earlier full-range mmag run is in
`vandaele_fits_fullrange_mmag.csv`), `digitize_vandaele.py` → `data_obs/vandaele_digitised{,_raw}.csv`, `vandaele_compare.py`
→ `data_obs/vandaele_compare.csv`. Light curves: Zenodo 20540863 on Ceph (`/mnt/ceph/users/mcantiello/rednoise/vandaele2026`),
281 PSF star-sectors, 21 empty → 260 usable, 91 stars. Placement: Bestenlehner+2025 A1 by position (84 with T_eff and L).*

## What we did
- Our fit follows their Sect. 6: one sector at a time, amplitude spectrum, least-squares semi-Lorentzian + white noise. The
  settings their figures reveal are adopted: flux units, 0.1–40 d⁻¹. Their manual Period04 prewhitening is replaced by an automatic
  S/N ≥ 5 (against the SLF model) iterative prewhitening. Cuts: ν_char ≥ 0.15, and ΔBIC against a white-noise line. Their
  visual rejection of 17 sectors cannot be reproduced. The fit metric (residuals in A or in log A) is not stated, so both were run.
- Their per-sector values were digitised from the Teff/L regression figure (3500 px original, markers found by sector colour,
  axes calibrated from ticks). A value is attributed to one of our star-sectors when it sits at that star's log T_eff and
  log L, with the same ordinate in both panels. N = 54 (ν_char) and 71 (α₀/C_w) star-sectors.

## Results
1. **Star by star: not reproduced.** Our log ν_char vs theirs has median +0.07 to +0.13 dex, MAD 0.19–0.26, 28–37% within 0.1 dex,
   Spearman ≤ 0.33 in every variant. For α₀/C_w the median offset is ≈ 0 with MAD 0.28 and Spearman 0.40–0.46.
2. **Why: single-sector ν_char is not a robust quantity.** On *identical* data, switching the fit metric alone changes log ν_char
   by 0.17 dex (MAD; only 32% agree within 0.1 dex). The frequency range moves it by 0.08. Prewhitening barely matters for SLF-
   dominated stars. Our sector-to-sector scatter for the same star is ≈ 0.24 dex, and their figure shows the same. The mismatch
   with Van Daele is the same size as our own choice-to-choice scatter, so star-level agreement cannot be expected without
   their code.
3. **Ensemble trends: partly reproduced, and dependent on the variant.**

| | N sectors | ν~T_eff slope | ν~L slope | α₀/C_w~L slope | α₀/C_w~T_eff |
|---|---|---|---|---|---|
| Van Daele (105 sectors, incl. cool stars) | 105 | **0.39** (p<0.01) | **0.08** (p=0.05) | **0.26** (p<0.01) | 0.05 (p=0.78) |
| ours, linear, prewhitened (42 OB stars) | 96 | **0.45** (p=0.007) | 0.09 (p=0.18) | 0.07 (p=0.75) | −1.28 (p=0.02) |
| ours, linear, raw | 98 | 0.65 | 0.10 | 0.05 | −1.21 |
| ours, log, prewhitened | 109 | 0.04 (p=0.74) | 0.01 | **0.34** (p=0.001) | −0.38 |
| ours, log, raw | 114 | 0.40 (p=0.006) | 0.03 | **0.28** (p=0.005) | −0.49 |

*(Corrected 2026-09-25: an earlier version of this table included 37 stars wrongly placed at one catalogue row; only the 54
stars matching Bestenlehner+2025 by position are used now. Van Daele also place cool supergiants, from other sources, which we
lack.)*

   No single variant reproduces all their trends. The linear fit recovers the ν_char–T_eff slope but not α₀/C_w–L; the log fit
   recovers α₀/C_w–L, and with prewhitening it loses the ν_char–T_eff trend.
   Medians: log ν_char −0.44 (ours) vs −0.54 (theirs); log α₀/C_w 1.34 vs 1.49.

## Consequences
- The Van Daele SMC values can be used only at the ensemble level, with ~0.2 dex per-sector noise and a method dependence of
  ≈ ±0.1 dex in the median. The metallicity comparison needs the **same pipeline on MW light curves**. Our own fitter can now do
  that: run `vandaele_fit.py` on TESS FFI/SPOC light curves of the Bowman 2020 / Shen 2024 Galactic stars, with the same range,
  metric and cuts. That is the homogeneous refit the paper calls for, and its infrastructure is ready.
- For the red-noise literature more broadly: published single-sector ν_char values carry an undocumented method systematic
  comparable to the 0.19 dex study-to-study scatter we measured earlier. The paper's 0.19 dex floor is appropriate, and it should
  not be read as the measurement error of any one study.
