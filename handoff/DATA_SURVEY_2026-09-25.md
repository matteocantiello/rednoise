# Data survey: published datasets not yet in the compilation

*2026-09-25 (corrected the same day: Serebriakova counts; the earlier arXiv link for it was wrong, the reference is A&A 692, A245). Searches: VizieR TAP keyword sweep of table descriptions (27 terms, 370 catalogs; `data_survey_raw/vz_sweep.txt`),
an author sweep (43 names; `vz_authors.txt`), the CDS ReadMes of every candidate (`data_survey_raw/readme/`), and a web/arXiv search
for 2022–2026 papers. Downloads so far are in `data_survey_raw/` (VizieR tables, Pedersen & Bildsten Zenodo tables). Nothing has been
merged into the analysis catalogs yet.*

The current compilation (handoff §2): red noise, MW 196 OB (µmag Lorentzian) + 76 YSG/RSG; LMC 47 and SMC 23 (GP / ppm, amplitudes
not comparable). Macroturbulence: 832 Galactic IACOB stars only. No microturbulence, no magnetic-star sample, no LMC/SMC
macroturbulence.

## Priority 1: directly unlock tests the paper cannot do now

| dataset | what it adds | size | access | caveats |
|---|---|---|---|---|
| **Serebriakova, Tkachenko & Aerts 2024**, A&A 692, A245 (`J/A+A/692/A245`) | **LMC macroturbulence** (LMC values from Serebriakova+2023, A&A 676, A85), measured consistently with a Galactic set; v sin i (two methods), T_eff, log g, spectroscopic luminosity | 141 LMC by coordinates (122 of the S23 subsample + 19 LMC members of G22; 136 with v_macro; log ℒ 1–4.65 with a gap at 2.5–3.6) + 69 Galactic (mostly low mass, log ℒ 0.95–2.8) | VizieR, downloaded (`vizier/serebriakova24.tsv`) | `Ls` is already log ℒ/ℒ⊙ (ReadMe unit "[K4/g]" is misleading); method differs from IACOB-BROAD, so cross-calibrate on Galactic stars in common before merging |
| **Van Daele+ 2026**, MNRAS (arXiv 2605.15757), BLOeM TESS PSF photometry | **Lorentzian SLF fits for 51 SMC stars** (ν_char in d⁻¹), with BLOeM spectroscopic parameters (Bestenlehner+2025) | 51 SMC, log L 4.5–5.8 | light curves on Zenodo 20540863; per-star parameters in the paper appendix (table A) | amplitude given only as α₀/C_w (relative to white noise), so ν_char and γ are comparable to the MW Lorentzian set, the amplitude is not |
| **Pedersen & Bildsten 2025**, MNRAS 539, 2742 (arXiv 2504.15861) | 50 new SLF stars (49 Cygnus OB + AV 232 SMC) **reaching log L = 1.7 (below the onset)**, plus a **refit of the Bowman 2020 sample (54 stars)** with the same method; model-independent rms (ppm) and ν₅₀ | 50 + 54 (refit) | Zenodo 15261328, downloaded (`pedersen_bildsten2025/`, MRT tables A1–B2) | α₀ is a power density (ppm²/µHz) and ν_char in µHz; the B20 overlap calibrates them onto the µmag/d⁻¹ primary sample |
| **Markova, Cantiello & Grassitelli 2025**, A&A 701, A297 (`J/A+A/701/A297`) | **microturbulence** with T_eff, log g, log ℒ: a third surface-velocity observable on the same axis (Cantiello+2009 relate v_mic to FeCZ velocities) | 1455 non-pulsating (374 hot, T_eff > 10 kK, log ℒ 1.2–4.3) + 368 pulsating (69 hot) | VizieR, downloaded | the compiled database also has v_macro and v sin i, but VizieR carries only T_eff, log g, v_mic, log ℒ, SpType |

## Priority 2: extend existing tests or add a clean differential test

| dataset | what it adds | size | access | caveats |
|---|---|---|---|---|
| **Holgado+ 2022**, A&A 665, A150 (`J/A+A/665/A150`) | IACOB-BROAD v_macro + v sin i for Galactic O stars; **62 O stars with v_macro not in our catalog** (same method as the IACOB compilation) | 285 (234 with v_macro) | VizieR, downloaded | merge with the existing priority order |
| **Shen+ 2023**, ApJ 955, 123 (`J/ApJ/955/123`) | SLF fits (µmag, same method as Shen 2024) for **118 magnetic hot stars**; they report ν_char anti-correlated with the dipole field B_p | 118 (189 sector fits) | VizieR | masses mostly 1.5–20 M☉, i.e. largely below the FeCZ onset; B_p from the magnetic literature (Shultz+2018 `J/MNRAS/475/5144`, Petit+2019 `J/MNRAS/489/5669`, MOBSTER VI `J/MNRAS/513/1429`). Test: at fixed ℒ, T_eff, magnetic stars with B_p above the FeCZ shut-off field (our `FeCZ_B_shutoff_conv`) should have lower α₀/ν_char |
| **Bowman & Dorn-Wallenstein 2022**, A&A 668, A134 (arXiv 2211.08347) | the **same 30 Galactic stars fitted with both Lorentzian and GP** (celerite2): the cross-calibration needed to compare MW Lorentzian ν_char with the Bowman+2024 LMC/SMC GP values | 30 | not in VizieR; per-star tables to be extracted from the paper (A&A pages block scraping; use arXiv) | without it the metallicity comparison of ν_char stays blocked (paper §4.6) |
| **Kourniotis+ 2025**, A&A 697, A152 (arXiv 2503.20860) | 41 Galactic blue supergiants with FEROS parameters and SLF (modified Lorentzian) fits; peak ~0.2 d⁻¹ (cool-side test) | 41 | tables in the paper appendix; Zenodo 15106993 has only figures | modified Lorentzian (like Ma+2024) |
| **Crawford+ 2026** (arXiv 2605.23209) | SLF (GP) in **16 extreme helium stars**; metal-poor EHe stars show none. FeCZ-like zones at low mass: an independent metallicity signature | 16 | Zenodo (tables) | GP SHO kernel; a different stellar class, use as a side test |

## Priority 3: context, placement, or future
- **BLOeM parameters** (Bestenlehner+2025, `J/MNRAS/540/3523`, 778 SMC OB stars; Lennon+2026 v sin i, `J/A+A/707/A204`): spectroscopic
  T_eff/log g/L for SMC stars, so tier-A placement of Bowman 2024 SMC stars and of the Van Daele sample.
- **Castro+ 2018** (SMC RIOTS4, `J/ApJ/868/57`, 329 OB) and **Castro+ 2021** (LMC NGC 2070 MUSE, `J/A+A/648/A65`, 333 OB): Magellanic sHRDs with
  log ℒ and v sin i (population context, placement).
- **Thomson-Paressant+ 2026** (MNRAS, arXiv 2608.03383): HERMES + TESS for 873 OB stars, 268 classified SLF-dominated, TLUSTY T_eff/log g/v sin i;
  SLF fits deferred to "Van Daele et al. in prep". Watch for it: it would be the largest homogeneous Galactic SLF sample.
- **Gebruers+ 2022** (`J/A+A/665/A36`): FEROS parameters incl. v_mic for 98 Galactic B stars with TESS light curves (low mass).
- **Holgado+ 2025** (`J/A+A/703/A175`): updated O-star calibrations (no broadening).
- **Pedersen+ 2019** (`J/ApJ/872/L9`): TESS variability classes for 154 OB stars (qualitative).
- **NGC 346 SMC O stars** with IACOB-BROAD v_macro (A&A 2025, "The O Vz stars in NGC 346"; Rickard+2022): small SMC v_macro sets, worth
  checking for tables.
- **LMC yellow supergiants, turbulent broadening** (arXiv 2508.14971): cool-side macro/micro-turbulence at LMC metallicity.
- **Gaia DR3 photometric dispersions** (Maíz Apellániz+2023, 145 M sources): a crude amplitude proxy for thousands of OB stars; noise- and
  cadence-limited, probably not competitive.
- Not useful: `J/MNRAS/466/1290` (SX Phe stars), VFTS O-star v sin i (Ramírez-Agudelo 2013; no v_macro in the CDS table).

## Suggested ingestion order
1. Serebriakova 2024 (LMC v_macro) + Holgado 2022 (MW O stars): cross-calibrate against IACOB on common Galactic stars, then test the LMC
   onset/saturation shift predicted by the v2 LMC grid (onset proxy log ℒ 3.62, Γ_Fe saturation 3.89; REPORT_rotation_Z.md §4).
2. Pedersen & Bildsten 2025: calibrate on the B20 overlap, add the 50 stars (below-onset coverage).
3. Van Daele 2026 + BLOeM parameters: SMC Lorentzian ν_char against the MW Lorentzian set.
4. Markova 2025 v_mic on the sHRD against the model FeCZ velocities.
5. Shen 2023 magnetic stars + B_p compilations: the magnetic suppression test.
6. Extract the Bowman & Dorn-Wallenstein 2022 two-method table (unlocks GP ↔ Lorentzian conversion for the Bowman 2024 LMC/SMC sets).


## Status
- 2026-09-25: Serebriakova 2024 and Holgado 2022 ingested for the LMC macroturbulence test: `../analysis_mesa/REPORT_lmc_vmac.md`.
- 2026-09-25: Pedersen & Bildsten 2025 ingested and calibrated (`../analysis_mesa/REPORT_pb25.md`): the paper's α₀ onset is not robust once the LMC ppm stars are removed.
- 2026-09-25: Van Daele 2026 SLF parameters are **not published per star** (figures only; Zenodo has light curves). Using them requires refitting the light curves, which the project's published-tables-only rule excludes; left for the user to decide.
- 2026-09-25: Markova 2025 microturbulence analysed (`../analysis_mesa/REPORT_vmic.md`).
- 2026-09-25 (later): the user lifted the published-tables-only rule for Van Daele. Their light curves are refitted
  (`../analysis_mesa/REPORT_vandaele.md`), the Galactic SPOC light curves too (`REPORT_mw_smc.md`), and a noise-floor model is
  built (`REPORT_floor.md`). Bowman & Dorn-Wallenstein 2022 digitised from their figures (`REPORT_nuchar_Z.md`). Magnetic
  test with Shen 2023 + Shultz 2022 + Petit 2013 (`REPORT_magnetic.md`). Bestenlehner+2025 A1 and Shen 2024 t1/t2 are in `vizier/`.
