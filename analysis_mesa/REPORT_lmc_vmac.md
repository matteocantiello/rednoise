# LMC macroturbulence: the first metallicity test of the FeCZ prediction

*2026-09-25. Code: `lmc_vmac_test.py` (hinge fits), `lmc_vmac_shift.py` (method offset, luminosity shift, figure).
Tables: `data_v2/lmc_vmac_test.csv`, `data_v2/lmc_vmac_shift.csv`. Figures: `figures_v2/fig_lmc_vmac.png`,
`figures_v2/fig_lmc_vmac_shift.png`. Data sources: `../handoff/DATA_SURVEY_2026-09-25.md`.*

## Data
- **MW:** the IACOB compilation (832 stars) plus the **58 Holgado+2022 O stars** (J/A+A/665/A150) not in it, which were measured with
  the same IACOB-BROAD method: 890 stars.
- **LMC:** Serebriakova+2024 (J/A+A/692/A245), galaxy assigned from coordinates: **141 LMC stars** (122 from Serebriakova+2023 and 19
  low-luminosity LMC members of the Gebruers+2022 set); 136 have v_macro, 135 of them within 4.0 < log T_eff < 4.53. (Correction to the
  survey: the S23 subsample also contains 2 Galactic stars, and G22 contains 19 LMC stars.)
- **No star is in both the Serebriakova and IACOB samples,** so the two v_macro scales cannot be tied star by star.

## Results
**1. The LMC sample cannot locate its onset.** It has 12 stars at log ℒ 1–2.5 and the rest at 3.6–4.65, with none in 2.5–3.6. The
hinge fit that finds the Galactic onset at log ℒ = 2.96 [2.90, 3.10] (ΔBIC +81, N = 890) finds no break in the LMC (ΔBIC −3.9).

**2. At fixed ℒ, LMC v_macro is lower by 0.26 dex** (log ℒ 3.6–4.0 and 4.0–4.4, same T_eff range), and it is still rising over
3.7–4.3, where the Galactic relation has flattened.

**3. Fitted as the Galactic relation shifted in luminosity,** log v_LMC(ℒ) = log v_MW(ℒ − Δ) + log k:

| fit | Δ log ℒ | log k |
|---|---|---|
| k free | +0.40 (bootstrap median +0.49; 75% in 0.2–0.7, 25% > 0.7, 0% < 0.2) | −0.08 [−0.11, +0.22] |
| k = 0 (no method offset) | **+0.57 [0.53, 0.62]** | 0 |

The v2 grids predict Δ = +0.60 for the onset proxy (v_c = 3 km/s) and +0.24 for the saturation proxy (Γ_Fe = 0.8)
(REPORT_rotation_Z.md §4). A metallicity-independent driver (e.g. core-excited waves) predicts Δ ≈ 0, which no bootstrap sample
gives.

**4. Method offset.** Serebriakova's Galactic stars against IACOB at matched log ℒ 2.0–2.8 and T_eff (below the onset in both):
−0.15 [−0.27, +0.03] dex (11 against 133 stars). This is consistent with zero but only weakly constrained; an offset of −0.26 dex,
about 1σ from this estimate, would by itself explain the whole LMC–MW difference.

## Reading
The LMC macroturbulence relation sits ≈0.4–0.6 dex higher in ℒ than the Galactic one, as the FeCZ grids predict for the
lower Fe opacity: close to the onset-proxy shift, and more than the saturation proxy alone predicts. It is the first observational
metallicity signature in this project that is not blocked by a fitting-method systematic, **with one caveat of the same kind**: the
v_macro measurement methods are not cross-calibrated, and the low-luminosity comparison constrains their offset only to ±0.15 dex.

## Caveats
1. The method offset (point 4). Settling it needs stars measured by both methods, or IACOB-BROAD applied to the LMC spectra (beyond
   the "published tables only" rule).
2. The LMC onset itself is not sampled (gap at log ℒ 2.5–3.6); the shift is measured on the rising and saturated part.
3. The samples differ in selection: Serebriakova's LMC stars are mostly B giants/supergiants, many beyond the MS; the MW set is
   IACOB's O/B mix. Both are restricted to 4.0 < log T_eff < 4.53, but no evolutionary-state or v sin i matching is done.
4. With k free, shift and scale are partly degenerate where the Galactic relation flattens (the upper tail Δ > 0.7).

## Next
- Look for LMC/Galactic stars measured by both methods (Serebriakova+2023 may compare with literature values), or a second LMC
  v_macro set (VFTS, XShootU) measured with IACOB-BROAD.
- Repeat with the SMC (NGC 346 O stars with IACOB-BROAD v_macro; SMC v2 grid when complete): predicted shift larger still.
- Predict the LMC relation directly from the transfer scenarios at each LMC star's position (LMC v2 grid), not only the shift.
