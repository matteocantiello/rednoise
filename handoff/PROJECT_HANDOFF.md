# Project handoff — Red noise and macroturbulence in massive stars

*Written 2026-09-22 for transfer to a second Claude Science instance that will run the MESA model comparison. Everything below was checked against the saved artifacts, not written from memory. Where a number appears, the CSV it comes from is named.*

---

## 0. How to use this document

1. Read §1–§3 to understand the question, the conventions, and what has been established.
2. The data products you need are listed in §5 and are packaged with this file in `project_handoff.tar.gz` (see §5.1). Load them; do not rebuild them.
3. §8 is the concrete MESA task: what to compute, what to compare it to, and where it goes in the manuscript.
4. §9–§10 contain the reproduction recipes and the pitfalls that cost time in this project. Read the pitfalls before touching any catalog.
5. Appendix A is a skill (`vizier-catalog-harmonize`) that exists only in the original project's skill store; re-create it on the new machine from the verbatim text if catalog work is needed.

---

## 1. Research question and standing decisions

**Question.** Does the stochastic low-frequency (red-noise) photometric variability of OB stars, and the macroturbulent line broadening seen in their spectra, scale with stellar parameters in the way predicted by (a) internal gravity waves (IGWs) excited by the convective core, or (b) turbulence in sub-surface convection zones driven by the iron opacity bump (FeCZ)? Can the two be distinguished?

**Decisions made by the user (binding):**
| decision | value |
|---|---|
| data source | aggregate **published tables only** — no new light-curve reduction or red-noise fitting |
| deliverable | report + figures + data table; then an **ApJ (AASTeX) manuscript** |
| luminosity axis | spectroscopic luminosity ℒ ≡ T_eff⁴/g with the user's normalisation `ell_sun = 5777**4/(274*100)`; `ell = (10**logTeff)**4/(10**logg)`; `logL_spec = log10(ell/ell_sun)` |
| HRD style | Cantiello et al. (2021) — properties mapped on the spectroscopic HRD with MIST tracks |
| interpolated maps | 2-D filled contours of expected values, **no data points**, reliability-weighted (rapid rotators down-weighted) |
| next step | user is collecting **MESA models**; steps 1–3 of the "next steps" list were run (onset, rotation, evolution); MESA comparison is the open item |

**Conventions adopted by the agent (documented, should be kept):**
- All ν_char in d⁻¹ (CoRoT µHz → d⁻¹ via 1 d⁻¹ = 11.574 µHz).
- α₀ has four incompatible unit systems across samples; amplitude comparisons are restricted to µmag-compatible stars (flag `a0_ok`).
- Primary homogeneous Galactic sample = Bowman et al. 2020 (70) ∪ Shen et al. 2024 (126 not already in Bowman) = **196 stars**. All scaling *coefficients* come from this sample.
- ℒ, T_eff, log g are algebraically degenerate (log ℒ = 4 log T_eff − log g − const); only the {ℒ, T_eff} basis is used in regressions.
- Published formal red-noise errors are unusable (10⁻⁵–10⁻² relative); the study-to-study scatter measured on 24 stars fitted by both Bowman 2020 and Shen 2024 (0.19 dex in log α₀ and log ν_char, 0.39 in γ) is the per-star systematic used everywhere.
- Z proxies: MW 0.014, LMC 0.006, SMC 0.002.

---

## 2. Data assembled

### 2.1 Red-noise samples (486 stars total; 342 placed on the sHRD)
| sample id in catalogs | source | N on sHRD | fit method | α₀ unit | used for α₀? |
|---|---|---|---|---|---|
| `Bowman2020_Galactic` | Bowman+2020, A&A 640, A36 (VizieR J/A+A/640/A36) | 70 | Lorentzian | µmag | yes |
| `Shen2024_Galactic` | Shen+2024, ApJS 275, 2 (J/ApJS/275/2 table1+table2) | 126 (+24 dup with Bowman) | Lorentzian | µmag | yes |
| `Bowman2019b_LMC` | Bowman+2019b Nat.Astron 3, 760 (J/other/NatAs/3.760/tabled24) | 13 | Lorentzian | ppm (unit field blank in ReadMe) | map/onset only, +0.05 dex sys |
| `Bowman2024_LMC`, `_SMC` | Bowman+2024, A&A 692, A49 (tables transcribed from PDF) | 21, 23 | GP (celerite2) | PSD-max proxy | **no** |
| `Ma2024_LMC_BSG` | Ma+2024, ApJ 966, 196 | 13 | modified Lorentzian | ppm | map/onset only, +0.15 dex sys |
| `DornWallenstein2020_YSG`, `_RSG` | Dorn-Wallenstein+2020, ApJ 902, 24 (J/ApJ/902/24/sample) | 28, 48 | damped random walk (τ) → ν_char = 1/(2πτ) | unconfirmed | **no** |
| `Bowman2019_CoRoT` (in master catalog only) | Bowman+2019a A&A 621, A135 | — | Lorentzian | PSD ppm²/µHz | no; not on sHRD (A/F pulsators) |
| `Bowman2019b_K2` (in master catalog only) | Bowman+2019b tabled13 | 17 matched | Lorentzian | ppm | no |

### 2.2 Macroturbulence (927 stars; 832 on sHRD) — `macroturbulence_catalog.csv`, `macroturbulence_evol.csv`
Priority order when a star is in several: Simón-Díaz+2017 (J/A+A/597/A22, 382 used) → de Burgos+2024 (J/A+A/687/A228, 376) → Holgado+2018 (J/A+A/613/A65, 62) → Burssens+2020 (J/A+A/639/A81, 12) → Simón-Díaz & Herrero 2014 (J/A+A/562/A135). IACOB-family sources agree pairwise to median Δ ≲ 1 km/s, ρ = 0.70–0.88.
**Excluded:** Shen+2024's `vmacro` column — offset −18 to −25 km/s vs IACOB, ρ ≤ 0.5, 32% of values > 150 km/s, 32 stars have vmacro but no v sin i → it is a total-broadening parameter with rotation absorbed. Shen's red-noise and atmospheric parameters *are* used.

### 2.3 Spectroscopic log g sources for Magellanic stars (tier A placement, 47 stars)
Bestenlehner+2025 (J/A+A/695/A198 tablea1/a2), Vink+2023 ULLYSES (J/A+A/675/A154/tableb2), Urbaneja+2017 LMC supergiants (J/AJ/154/102 BA/OB1/OB2/OB3; publishes flux-weighted log g_F → log g = log g_F + 4 log(T_eff/10⁴ K)), McEvoy+2015 VFTS BSGs (J/A+A/575/A70/stars). Name match after normalisation + 2″ coordinate pass. Matched T_eff agree with the red-noise papers' T_eff to median 0.000 dex (MAD 0.005).

### 2.4 Tier B placement (99 stars: Bowman2024 unmatched, Ma unmatched, all Dorn-Wallenstein)
ℒ/ℒ⊙ = (L/L⊙)/(M/M⊙) exactly. Mass from nearest MIST point in (log T_eff, log L), 8-nearest inverse-distance-squared weighting, scaling (0.02, 0.10). Calibrated on 33 stars with both routes: **offset −0.089 dex (applied), scatter 0.186 dex (propagated)**. Figure A2 of the paper.

### 2.5 MIST tracks
MIST v1.2, [Fe/H]=0, v/v_crit=0.4, EEP tracks. Downloaded from `https://mist.science/data/tarballs_v1.2/MIST_v1.2_feh_p0.00_afe_p0.0_vvcrit0.4_EEPS.txz` (112 MB; `waps.cfa.harvard.edu` redirects there). **Filename encoding is mass × 100** (`01500M.track.eep` = 15 M⊙). Columns read from the `#` header line containing `star_age`: `log_L`, `log_Teff`, `log_g`, `star_mass`, `phase`. MS = phase 0 (and −1 pre-MS excluded). Fractional MS age τ ∈ [0,1] on the MS, >1 beyond. No 25 M⊙ track — 24 used.

### 2.6 Evolutionary parameters assigned to every star (`ev_*` columns)
τ, log M, log R, phase, initial mass, and a distance-to-track `ev_dist` (exclude if ≥ 0.5; 9 of 342). Validated vs Shen+2024's own published τ and M on 126 stars: Δlog M median 0.00 (MAD 0.02), Δτ median +0.08 (MAD 0.04; different track physics). Projected rotation frequency `f_rot = vsini/(2πR)` in d⁻¹ (lower limit).

---

## 3. Results established (numbers from the saved tables)

### 3.1 Scaling coefficients — `multivariate_fits_v2.csv` (primary sample, bootstrap 16–84%)
| target | model | coefficient |
|---|---|---|
| log α₀ | ℒ + T_eff (N=196, R²=0.53) | ∂/∂log ℒ = **+1.56** [+1.43, +1.68]; ∂/∂log T_eff = −3.16 [−3.83, −2.55] |
| log α₀ | + log v_macro (N=144) | ∂/∂log v_macro = **+0.28** [+0.15, +0.47] (independent of ℒ, T_eff) |
| log ν_char | ℒ + T_eff (N=196, R²=0.09) | ∂/∂log T_eff = **+1.53** [+1.20, +1.88]; ∂/∂log ℒ = −0.23 [−0.32, −0.14] |
| log ν_char | + log v sin i (N=152) | ∂/∂log v sin i = +0.34 [+0.25, +0.42] |
| log α₀ | τ + log M, Shen's own τ,M (N=126, R²=0.45) | τ: +1.91 [+1.62, +2.23]; log M: +1.83 [+1.57, +2.10] |

### 3.2 HRD-plane gradients — `hrd_gradients.csv`
| property | ∂/∂log ℒ | ∂/∂log T_eff | N | R² |
|---|---|---|---|---|
| v_macro | +0.29 | +0.10 | 832 | 0.44 |
| α₀ | +1.54 | −3.13 | 196 | 0.53 |
| ν_char | −0.23 | +1.53 | 196 | 0.09 |
| γ (linear) | +0.15 | −2.10 | 196 | 0.10 |
**Interpretation:** α₀ and v_macro share the ℒ coordinate; ν_char is set by T_eff and flat in ℒ; γ ≈ 1.9 (16–84%: 1.43–2.45) nearly universal.

### 3.3 Luminosity onset — `fecz_onset_test.csv` (broken power law in log ℒ with log T_eff covariate; ΔBIC = single − broken; break claimed if > 6)
| quantity | N | hinge log ℒ [16–84%] | slope below → above | ΔBIC | verdict |
|---|---|---|---|---|---|
| v_macro, log T_eff > 4.3 | 579 | 3.80 [3.70, 3.84] (one hinge); **two-hinge 3.16 & 3.70 preferred by a further ΔBIC=17.5**, slopes +0.37/+0.81/+0.29 | +0.66 → +0.27 | +14.9 | break |
| v_macro, all IACOB | 832 | 2.96 [2.90, 3.10] | −0.18 → +0.44 | +95 | break |
| v_macro, 4.0 < log T_eff ≤ 4.3 | 253 | 2.78 [2.68, 3.08] | −0.08 → +0.23 | +4.8 | **not significant** |
| α₀ (OB) | 222 | **3.16 [2.94, 3.44]** | −0.02 → +1.66 | +16.3 | break (19 stars below; bimodal) |
| ν_char (OB) | 271 | — | — | −0.6 | no break |
Floor below onset: hot-star v_macro ≈ 16 km/s (median at ℒ < 2.9); α₀ ~ 30 µmag (poorly constrained). 1-D models (Cantiello+2009, 2021) first form an FeCZ at log ℒ ≈ 2.5 at solar Z, i.e. **0.6 dex below the observed onset**. On MIST, log ℒ = 3.16 at log T_eff > 4.3 ≈ 12 M⊙, log L/L⊙ ≈ 4.0.

### 3.4 Evolution at fixed mass — `evolution_test_full.csv` (MS stars, τ < 1, ev_dist < 0.5; model log(·) = a τ + b log M)
| target | N | R² | ∂/∂τ | ∂/∂log M |
|---|---|---|---|---|
| log α₀ | 190 | 0.51 | **+1.73** [+1.57, +1.91] (×50 ZAMS→TAMS) | +1.32 [+1.11, +1.52] |
| log ν_char | 221 | 0.03 | −0.33 [−0.48, −0.19] | +0.12 [+0.03, +0.22] |
| log v_macro | 562 | 0.45 | +0.27 [+0.22, +0.32] (×1.9) | +0.60 [+0.57, +0.63] |
Post-TAMS: growth stops — 28 post-MS OB stars median α₀ = 877 µmag vs 1211 µmag for 82 MS stars with τ > 0.8 (MWU p = 0.22).
**Amplitude-to-velocity sensitivity ratio:** 5.3 along ℒ, 6.4 along τ, **but only 2.2 along mass**. Direct power law α₀ ∝ v_macro^1.37 (N=143).

### 3.5 Rotation — `rotation_test.csv`
- ν_char / f_rot: median **16.7** (16–84%: 6–40); no pile-up at f_rot or 2 f_rot → rotational modulation is not the turnover. 11 stars within 0.15 dex of f_rot or 2 f_rot flagged `rotmod_candidate`.
- Residual slopes vs log v sin i after ℒ, T_eff: ν_char +0.34 (partial ρ +0.19, p 0.018, N 164); α₀ +0.19 (ρ +0.17, p 0.031, N 164); v_macro +0.14 (ρ +0.29, p 3×10⁻¹⁷, N 832 — opposite to the fitting-degeneracy direction, hence secure).

### 3.6 Cool side — regional medians (from `rednoise_sHRD_extended_evol.csv`)
| region | N | ν_char d⁻¹ [16–84%] | γ |
|---|---|---|---|
| hot OB (log T_eff ≥ 4.3) | 235 | 1.47 [0.45, 3.80] | 1.91 |
| evolved OB/BSG (4.0–4.3, ℒ > 3.6) | 22 | 0.27 [0.20, 0.46] | 2.58 |
| YSG (3.75–4.0) | 20 | 0.30 [0.17, 0.50] | 1.63 |
| RSG (< 3.75) | 51 | 0.59 [0.33, 2.33] | 1.72 |
BSG/YSG share a ~0.3 d⁻¹ floor; RSGs are a distinct regime (read descriptively).

### 3.7 Metallicity — **method-confounded, unresolved**
LMC vs SMC (both GP-fit): medians 1.04 vs 0.96 d⁻¹, indistinguishable. MW (Lorentzian) higher even at matched T_eff (1.85 vs 1.04 in 4.4 < log T_eff < 4.65, MWU p = 0.011) but Lorentzian vs GP ν_char differ systematically (Bowman & Dorn-Wallenstein 2022). No cross-galaxy frequency comparison is made. A homogeneous refit of all light curves would unlock this — out of scope under the "published tables only" decision.

**Update 2026-09-25:** superseded. The GP↔Lorentzian offset is only −0.02 dex (BDW22), and a homogeneous refit (SPOC Galactic + Van Daele SMC PSF light curves) with a noise-floor model now exists. See `../STATUS.md` and `../analysis_mesa/REPORT_{nuchar_Z,mw_smc,floor}.md`.

### 3.8 Theory verdict
FeCZ predicts: amplitude ∝ position toward Eddington limit (ℒ), ν_char set by turnover time (T_eff), an onset, growth with evolution, α₀–v_macro link. Core IGW predicts: amplitude grows slowly with mass, spectrum softens/shifts lower with age, no α₀–rotation link, amplitude 0.06–0.21 µmag for 15 M⊙ (Anders+2023) vs observed 10–90% of 64–3600 µmag (median 350) → 3–4 dex short. **Sub-surface convection is the dominant driver of both phenomena in evolved OB stars.** A metallicity-independent IGW floor below the onset remains possible and untested. Full prediction compilation: `theory_predictions.md`.

---

## 4. Manuscript status

`manuscript.pdf` (17 pp, AASTeX 6.3.1, two-column) and `paper_bundle.tar.gz` (main.tex, refs.bib, aastex631.cls, aasjournal.bst, tables/, figures/, README.md, outline.md). Compiles clean: 0 errors, 0 overfull boxes, 0 unresolved refs. 79 numerical claims cross-checked by script against the CSVs above.

**Title:** Surface turbulence and stochastic low-frequency variability in massive stars share a sub-surface driver.

**Figures:** 1 samples on sHRD; 2 scalings (primary sample); 3 GP maps 2×4; 4 onset; 5 evolution at fixed mass; 6 rotation; 7 theory scorecard + amplitude budget + added-variable test; A1 Bowman-vs-Shen overlap; A2 tier-B calibration; B1 correlation matrix. All regenerated from artifacts by `make_figures.py`.

**Tables:** tab_samples, tab_vmac, tab_fits, tab_onset, tab_regions (deluxetable); `table_rednoise_shrd.mrt` (342 rows), `table_macroturbulence.mrt` (832 rows) in AAS MRT format.

**Placeholders (red `\todo{}`):** author list/affiliations; **§5.6 "Comparison with envelope models"** (this is the MESA slot); acknowledgments; `\facilities` list to check.

**Bibliography caution:** 37 refs, all DOIs verified via CrossRef. Nine DOIs recalled from memory pointed to the *wrong papers* and were fixed by title search. Never add a reference without resolving its DOI. `pauli2025`, `hawcroft2024`, `maeder2008` were dropped (not resolvable / wrong hit); Bestenlehner+2025 A&A 695, A198 covers the XShootU parameters.

---

## 5. Data products to transfer

### 5.1 The transfer package
`project_handoff.tar.gz` contains this document plus the files below. On the new machine, untar into the workspace and `save_artifacts` the data files so they persist across workspace cleanups.

| file | rows | what it is | key columns |
|---|---|---|---|
| `rednoise_sHRD_extended_evol.csv` | 342 | **the red-noise analysis sample** | `star, sample, lT (log Teff), lL (log ℒ), tier, alpha0, nuchar (d⁻¹), gamma, vsini, a0_ok, sys_nu, sys_a0, ev_tau, ev_logM, ev_logR, ev_phase, ev_Mini, ev_dist, ev_M, ev_R, f_rot, log_alpha0, log_nuchar, rotmod_candidate` |
| `macroturbulence_evol.csv` | 832 | **the macroturbulence analysis sample** | `key, vmac, vsini, vmac_src, SpType, logTeff_sp, logL_sp, logg_sp, ev_tau, ev_logM, ev_logR, ev_phase, ev_Mini, ev_dist, ev_R, f_rot` |
| `master_rednoise_catalog_v2.csv` | 486 | full red-noise compilation incl. CoRoT/K2, `dup_B2020_Shen` flag, `alpha0_kind` | |
| `macroturbulence_catalog.csv` | 927 | unified IACOB compilation (v2, Holgado kK fix applied) | |
| `hrd_interpolated_grids_v2.npz` / `.csv.gz` | 150×150 ×2 grids | GP posterior fields: `{vmac,alpha0,nuchar,gamma}_{mean,std,ok}`, `logTeff_OB, logL_OB` (4.03–4.72 × 1.90–4.50), `logTeff_ext, logL_ext` (3.50–4.78 × 1.90–4.70). `mean` is log10 for vmac/alpha0/nuchar, linear for gamma. `ok` = trustworthy mask | |
| `hrd_predictions_tracks.csv` | 15 | expected v_macro/α₀/ν_char/γ ± sd at near-ZAMS/mid-MS/TAMS on 9, 15, 20, 32, 60 M⊙ MIST tracks (OB-grid maps only; NaN where untrustworthy) | |
| `multivariate_fits_v2.csv`, `hrd_gradients.csv`, `fecz_onset_test.csv`, `evolution_test_full.csv`, `rotation_test.csv`, `correlations_summary_v2.csv`, `correlation_matrix_v2.csv` | — | the result tables quoted in §3 | |
| `theory_predictions.md` | — | IGW vs FeCZ prediction compilation with diagnostic table | |
| `report_v2.md` | — | the long-form analysis report (v4) with method details | |
| `paper_bundle.tar.gz` | — | manuscript source, figures, tables | |
| `make_figures.py` | — | regenerates all paper figures (references artifact IDs — see §9.2) | |
| `bowman2020_A36_tablea1.tsv`, `bowman2020_A36_tablea2.tsv` | 70 | raw VizieR downloads of the flagship sample | |

### 5.2 Artifact IDs (valid only inside the original project `proj_ff77e5035617`)
Full inventory with `latest_version_id` is in `artifact_inventory.csv` inside the package. The ones `make_figures.py` hard-codes:
`rednoise_sHRD_extended_evol.csv` 114af09b-6e33-4161-a82e-8948b966202a · `macroturbulence_evol.csv` dd29c64d-df65-4121-9718-725a60a3a01e · `master_rednoise_catalog_v2.parquet` 000f7d1a-4ac6-4eec-81a1-3629b5611e2f · `hrd_interpolated_grids_v2.npz` b8244da2-cc76-42f7-8a0d-5fce81fd5677 · `fecz_onset_test.csv` ea99664e-a9aa-462f-bc51-f4d1a5da0aaa · `evolution_test_full.csv` bb4fef82-2184-4fc8-a891-23a8f0a854b0 · `rotation_test.csv` 9f5de9cc-6249-4398-a75d-4f4283c9cf6a · `multivariate_fits_v2.csv` 8bb974e0-2801-4ae8-92cb-b0f326559d01.

---

## 6. Methods in enough detail to reproduce

**Spearman + BIC + bootstrap.** Rank correlations per sample. Linear models in log space; exhaustive subset selection by BIC; 3000–5000 bootstrap resamples for 16–84% intervals. VIF < 2.5 on {ℒ, T_eff, v sin i}.

**Reliability weights (GP maps).** σ_i² = σ_sys² + σ_rot,i² + σ_B,i². σ_sys: 0.186 dex (log α₀), 0.192 dex (log ν_char), 0.391 (γ), 0.08 dex (log v_macro, ~20% IACOB pairwise). σ_rot = 0.35 · vsini/(vsini + 100 km/s) dex; unknown v sin i → 0.175. σ_B (tier B only) = 0.19 dex × |∂property/∂log ℒ| (GRAD = {alpha0: 1.54, nuchar: 0.23, gamma: 0.15}). Plus per-sample `sys_a0` (0.05 Bowman2019b, 0.15 Ma) and `sys_nu` (0.15 for Dorn-Wallenstein τ→ν conversion). Sensitivity: dropping the rotation term changes maps by median < 0.02 dex.

**GP.** sklearn `GaussianProcessRegressor`, kernel `C(1.0) * Matern(length_scale=[1,1], nu=2.5)`, length-scale bounds (0.3, 10) in standardized units, `alpha = (σ_i / y_std)²`, **standardize y explicitly and divide σ by y_std** (using `normalize_y=True` while passing raw σ² makes the noise ~13× too small and collapses the length scales — this bug was hit). No WhiteKernel (it inflates the reported sd on the latent mean). `n_restarts_optimizer=4`. Trustworthy = (posterior sd < 0.6 × y_std) AND (≥ 3 stars within ellipse semi-axes 0.06 in log T_eff, 0.25 in log ℒ).

**Onset test.** log(·) = c₀ + c₁(x − h) + c₂ max(0, x − h) + c₃ log T_eff, x = log ℒ; h scanned 2.2–4.3 step 0.02 with ≥ 10 stars each side; 300 bootstrap resamples for h interval; ΔBIC = BIC_single − BIC_broken; claim if > 6 and hinge interior.

**Rotation test.** f_rot = vsini · 86400 / (6.957e5 · 2π · R[R⊙]) d⁻¹. Partial residuals: regress y and log vsini on [log ℒ, log T_eff], correlate residuals.

**Figure conventions.** `figure-style` skill applied; viridis for T_eff / ℒ colouring, plasma for mass, cividis for ν_char, magma for γ; double-column 7.25 in; PDF + PNG; text-overlap self-check after `savefig`.

---

## 7. Pitfalls encountered (read before touching any catalog)

1. **Guessed VizieR IDs resolve to unrelated tables.** Discover by TAP keyword search of `TAP_SCHEMA.tables.description` (`tapvizier.cds.unistra.fr`). `astroquery.Vizier.get_catalogs` votable path **hangs** in this sandbox — use `curl` to `https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=CAT&-out.max=unlimited&-out.all=2&-out.add=_RAJ2000,_DEJ2000`.
2. **Always read the CDS ReadMe for units** (`https://cdsarc.cds.unistra.fr/ftp/<CAT>/ReadMe`). Bugs found: Holgado+2018 T_eff in **kK**; Urbaneja+2017 T_eff in **10⁴ K** (a "<100 → kK" heuristic mis-scaled it); Bowman+2019b α₀ unit field blank.
3. **Formal errors in red-noise tables are meaningless** as weights; measure the systematic from stars in common between studies.
4. **Shen+2024 `vmacro` is not macroturbulence** (see §2.2).
5. **VizieR TSV row counting**: `grep -vc '^#'` overcounts by the blank/header/unit/dash lines (70 stars, not 75). Parse by finding the dash separator row.
6. **Name normalisation** needed for: `Cl* NGC 346 MPG 12`, `AzV 083`→`AV83`, `HDE`→`HD`, `Sk -67 171`, `SK--65 2` (double dash), LaTeX-laden aliases, comma/`=` compound alias fields.
7. **MIST filenames are mass×100**; no 25 M⊙ track; extract tarball with `tarfile`, files may land mode 000 — `chmod -R u+rwX`.
8. **astropy config dir** is not writable in the sandbox: set `XDG_CONFIG_HOME`, `XDG_CACHE_HOME`, `ASTROPY_CONFIGDIR`, `ASTROPY_CACHEDIR` to workspace `.cache/` subdirs before importing astropy.
9. **Matplotlib**: `errorbar` does not accept `ec=`; mathtext has `\leq` not `\le`; never `plt.savefig` (use `fig.savefig`); re-check text overlaps *after* `savefig` (layout finalises then); keep `fig/axes` names unique when several figures are alive (a panel was overwritten once).
10. **Bibliography**: DOIs recalled from memory were wrong for 9 of 36 papers. Resolve every DOI via CrossRef (`https://api.crossref.org/works/<DOI>`), title-search (`query.bibliographic=`) when a DOI 404s, and check the returned title. CrossRef returns CJK author annotations for some names — strip them or bibtex emits U+FFFD.
11. **LaTeX toolchain in the sandbox**: conda `texlive-core` (macOS) is binaries-only; `tlmgr` broken; `tectonic` cannot verify the sandbox proxy TLS certificate (uses macOS Security framework, ignores `SSL_CERT_FILE`); local sockets cannot be bound (no relay). **Working recipe**: `curl` the bundle tar from `https://data1b.fullyjustified.net/tlextras-2022.0r0.tar` (2.9 GB; `relay.fullyjustified.net` redirects there), extract to a directory, `chmod -R u+rwX`, then `tectonic -X compile --bundle <dir> main.tex` with `TECTONIC_CACHE_DIR`/`XDG_CACHE_HOME` pointed into the workspace. The bundle directory contains `aastex631.cls` and `aasjournal.bst`. Delete `main.bbl`/`main.aux` after editing `refs.bib`.
12. **Network domains that had to be requested**: vizier.cds.unistra.fr, tapvizier.cds.unistra.fr, cdsarc.cds.unistra.fr, waps.cfa.harvard.edu, mist.science, journals.aas.org, relay.fullyjustified.net, data1b.fullyjustified.net. A&A (aanda.org) serves a DataDome captcha — unscrapable; use arXiv or Unpaywall via `fetch_article_fulltext(doi=…)`. The `literature` MCP connector was **not** attached in the last session; CrossRef REST worked without credentials.
13. **Kernel resets / workspace sweeps** happen between sessions. Everything needed is an artifact; restore with `host.artifacts(filename=…, exact=True)` → `host.artifact_path(vid)`. Do not rebuild.

---

## 8. The MESA task (open item — what the second instance should do)

### 8.1 Purpose
Fill manuscript §5.6 ("Comparison with envelope models", currently a `\todo{}`) by asking whether a grid of 1-D MESA models with an FeCZ reproduces the three empirical numbers. The result decides whether the sub-surface interpretation is "favoured" (current wording) or "quantitatively confirmed".

### 8.2 Model quantities to extract, per model and per timestep (MS + early post-MS, phase ≤ 2)
From each profile/history: T_eff, log g, log L, M, age → **log ℒ = 4 log T_eff − log g − log ℓ⊙** with `ell_sun = 5777**4/(274*100)` (i.e. log ℒ/ℒ⊙ = log(T_eff⁴/g) − log(5777⁴/27400)); τ = fractional MS age (define TAMS as X_c < 10⁻³ or by MESA's `center_h1`). For the FeCZ (the convective region whose base lies near log T ≈ 5.3, the Fe opacity bump):
- presence / absence;
- **v_c,max** (maximum MLT convective velocity in the zone) and the mixing-length used (α_MLT);
- **F_c/F_tot,max** (maximum convective flux fraction);
- turnover time **t_c = α_MLT · H_P / v_c** at the location of v_c,max → **ν_c = 1/(2π t_c)** in d⁻¹ (Cantiello+2021 definition);
- zone depth below the photosphere (in R and in mass), density at v_c,max;
- optionally the turbulent pressure fraction P_turb/P and the He II convection zone properties for the cool side.

### 8.3 The three comparisons (in order of diagnostic value)
1. **Onset.** Plot v_c,max and F_c/F_tot against log ℒ for models with log T_eff > 4.3. The observations show v_macro flat (≈16 km/s) below log ℒ = 3.16 [2.94, 3.44] and α₀ flat (≈30 µmag) below the same point; 1-D models first form the FeCZ at ≈ 2.5. **Question:** at what log ℒ does v_c,max cross ~1–3 km/s, and does F_c/F_tot become non-negligible (say > 1%)? If that happens at ≈ 3.0–3.2, the onset is a *threshold* effect and the 0.6-dex gap is explained. Compare on the same axes as Fig. 4 (`fecz_onset_test.csv`; the binned medians can be recomputed from the CSVs).
2. **Saturation.** v_macro flattens from ∝ ℒ^0.81 to ∝ ℒ^0.29 above log ℒ = 3.70 [3.70, 3.84] while α₀ keeps rising as ℒ^1.66. **Question:** does v_c,max saturate near there (sound speed? envelope inflation?) while F_c/F_tot continues to rise? Report the model slopes d log v_c/d log ℒ and d log(F_c/F)/d log ℒ in the two regimes.
3. **Sensitivity ratios.** Observed: ∂log α₀/∂log ℒ ÷ ∂log v_macro/∂log ℒ = 5.3; along τ = 6.4; along mass = 2.2. Naïve expectation with α₀ ∝ F_c/F ∝ ρ v_c³/T_eff⁴ and v_macro ∝ v_c gives 3. **Question:** what does the model give for d log(F_c/F)/d log v_c along each axis, and does ρ at the FeCZ increase along ℒ and τ enough to raise it to 5–6 (and why not along mass)?

Secondary checks: (4) **ν_c vs ν_char** on the sHRD — Cantiello+2021 found model ν_c ≈ 3× lower than observed; test whether the offset is constant across the diagram (compare to the GP `nuchar_mean` field in `hrd_interpolated_grids_v2.npz`, trustworthy mask `nuchar_ok`) — a constant offset supports the turnover-time identification, a varying one does not; the model should also give ∂log ν_c/∂log T_eff ≈ +1.5 and ≈ 0 in ℒ. (5) **Evolution at fixed mass**: model d log v_c/dτ and d log(F_c/F)/dτ vs observed +0.27 and +1.73 (α₀ ×50 across the MS). (6) **Cool side**: where does the FeCZ hand over to the He II / H zones, and do their turnover frequencies give ~0.3 d⁻¹ at 4.0 < log T_eff < 4.3?

### 8.4 Deliverables for the paper
- A figure (proposed Fig. 8, double column, 3 panels): (a) model v_c,max and F_c/F vs log ℒ with the observed binned v_macro/α₀ overplotted and the onset/saturation hinges marked; (b) model ν_c contours on the sHRD over the observed ν_char GP field; (c) model vs observed sensitivity ratios along ℒ, τ, M.
- A table: for each of the three targets, observed value ± interval, model value, and the α_MLT (and any other free parameter) required.
- Text for §5.6 replacing the `\todo{}` (keep the three-numbered structure already there). Update the abstract's last sentence ("identify … as quantitative targets") to state the outcome.
- Add `mesa_*.csv` outputs and the model grid description (masses, Z, rotation, overshoot, α_MLT, MESA version, opacity tables) to §2 or an appendix; cite MESA instrument papers (Paxton et al. 2011, 2013, 2015, 2018, 2019; Jermyn et al. 2023) — resolve DOIs via CrossRef before adding to `refs.bib`.

### 8.5 Regenerating and recompiling
`make_figures.py` hard-codes artifact IDs from the original project; on the new machine, replace each `host.artifact_path("…")` with the path to the corresponding file from the transfer package (the mapping is in §5.2). Compile with the recipe in §7.11 or any standard TeX Live (`latexmk -pdf main`).

---

## 9. Reproduction recipes

### 9.1 Restore state on a fresh kernel (original project)
```python
import pandas as pd, numpy as np
ap = host.artifact_path
EXT = pd.read_csv(ap("114af09b-6e33-4161-a82e-8948b966202a"))   # 342 red-noise stars
VM  = pd.read_csv(ap("dd29c64d-df65-4121-9718-725a60a3a01e"))   # 832 macroturbulence stars
GR  = np.load(ap("b8244da2-cc76-42f7-8a0d-5fce81fd5677"))       # GP grids
X   = pd.read_parquet(ap("000f7d1a-4ac6-4eec-81a1-3629b5611e2f"))# 486-star master
```
On the new machine: same, with local file paths from the package.

### 9.2 Primary sample definition
```python
P = X[(X["sample"]=="Bowman2020_Galactic") | ((X["sample"]=="Shen2024_Galactic") & ~X.dup_B2020_Shen)]   # N = 196
```

### 9.3 Spectroscopic luminosity
```python
ell_sun = 5777.0**4/(274*100)
spec_ell = lambda logT, logg: np.log10((10**np.asarray(logT))**4/(10**np.asarray(logg))/ell_sun)
# check: spec_ell(np.log10(35800), 3.86) == 3.746 (Shen publishes 3.75)
```

### 9.4 MIST-based τ, M, R (as used)
Build a point cloud from all EEP tracks 5–120 M⊙, phases 0–6, columns (log_Teff, log_L or ℒ, star_mass, phase, EEP index → τ); `cKDTree` on scaled coordinates (0.02 in log T_eff, 0.10 in log L); k=8 inverse-square-distance weights; `ev_dist` = nearest scaled distance. τ on the MS from the EEP number between ZAMS (EEP 202) and TAMS (EEP 454) — see `report_v2.md` §3b for details.

---

## 10. People/agent roles and review process used
- The user set the science direction and three binding decisions (§1). The agent did all data work.
- A background reviewer ("Auditor") repeatedly checked numeric claims against tool outputs; findings that were correct were fixed (sample-size miscount, unbacked coefficients, mis-labelled subsample row, over-general "ratio ≈6" claim, inconsistent amplitude-inclusion statement); one was rebutted with evidence. **Expect the same on the new machine: keep every number traceable to a CSV.**
- Workflow habit that paid off: write long analysis code to a `.py` file and `exec(open(...).read())`; check figure text overlaps after saving; script-verify manuscript numbers against tables.

---

## Appendix A — Skill `vizier-catalog-harmonize` (verbatim; recreate via `skill("customize")` → `host.skills.edit(name, "SKILL.md", …)`, `host.skills.edit(name, "kernel.py", …)`, `host.skills.publish(name)`)

### A.1 SKILL.md
```
---
name: vizier-catalog-harmonize
description: "Discover, download, unit-check and cross-match published stellar catalogs from VizieR/CDS by keyword search of table metadata (TAP/ADQL), then harmonize them onto a common axis (spectroscopic luminosity L=Teff^4/g, flux-weighted gravity, name/coordinate matching). Use when adding a literature table (red noise, macroturbulence, Teff/logg, any star-by-star parameters) to an analysis without guessing catalog identifiers."
---

# VizieR catalog discovery → download → harmonize

Kernel helpers (loaded automatically): `vz_tap_search`, `vz_fetch`, `vz_read_tsv`,
`vz_readme_units`, `star_key`, `split_aliases`, `alias_index`, `spec_lum`,
`logg_from_loggF`, `xmatch_coords`, `astropy_config_fix`.

## 0. Network domains (request once per project)
`vizier.cds.unistra.fr` (asu-tsv download), `tapvizier.cds.unistra.fr` (TAP metadata
search), `cdsarc.cds.unistra.fr` (ReadMe unit definitions). For MIST tracks:
`waps.cfa.harvard.edu` redirects to `mist.science` — request both.

## 1. Discover by metadata search — never guess catalog IDs
```python
hits = vz_tap_search(["macroturbulen", "IACOB", "red noise", "stochastic low-frequency"])
# -> list of (table_name, description). Also try author surnames ("Bowman D.M.").
```
Guessed `J/A+A/vol/page` IDs resolve to unrelated tables; `astroquery.Vizier.find_catalogs`
does fuzzy matching and returns hundreds of junk hits. Search `TAP_SCHEMA.tables.description`
with `LIKE '%term%'` instead. To list every table of one catalog:
`vz_tap_search([], table_like="J/ApJS/275/2")`.

## 2. Download with the plain asu-tsv endpoint
```python
df = vz_fetch("J/A+A/640/A36/tablea2", "bowman2020_t2", add_coords=True)
```
`astroquery.Vizier.get_catalogs` (votable) hangs in the sandbox; curl to
`viz-bin/asu-tsv` works. `add_coords=True` appends `_RAJ2000,_DEJ2000` — many tables
lack positions by default. Row count: count rows *after* the dashed separator line; a
naive `grep -vc '^#'` also counts header/unit/blank lines.

## 3. Verify units from the ReadMe before merging anything
```python
vz_readme_units("J/A+A/613/A65", ["Teff", "logg", "alpha0", "vmac"])
```
Cases that bit: Teff in **kK** (Holgado+2018, Bestenlehner+2025, Pauli+2025), in **K**
(Vink+2023), in **10^4 K** (Urbaneja+2017) — a `<100 → kK` heuristic mis-scales the
10^4 K case. Amplitude columns come as µmag, ppm, ppt (`10-3`), PSD (ppm²/µHz) or a
PSD-maximum proxy — flag each with an `alpha0_kind` column and only compare within kind.
A blank unit field (`---`) means *unconfirmed*, say so. Published `e_` columns may be
formal MCMC errors (10⁻⁵ relative): when two studies fit the same stars, measure the
real systematic from the overlap (0.19 dex for red-noise parameters) and use that as the
weight floor instead.

## 4. Name matching
```python
idx = alias_index(df, "Alias")           # handles "AV14,SK9", LaTeX, "SK--65" double dash
row = idx.get(star_key("Cl* NGC 346 MPG 12"))
```
`star_key` upper-cases, strips whitespace/`Cl*`/LaTeX, maps AzV→AV, HDE→HD, collapses
`--`, and drops leading zeros in numeric IDs. Always report matched/total per sample.

## 5. Coordinate matching (2″ default)
```python
astropy_config_fix()                     # sandbox blocks ~/.astropy — redirect first
m = xmatch_coords(df_a, "_RAJ2000", "_DEJ2000", df_b, "_RAJ2000", "_DEJ2000", radius_arcsec=2.0)
```
Drops NaN positions before matching (astropy raises otherwise). Returns index into
`df_b` and separation, NaN where unmatched.

## 6. Common luminosity axis
Spectroscopic luminosity ℒ = Teff⁴/g, normalised so `spec_lum(log Teff, log g)` returns
log(ℒ/ℒ⊙) with ℒ⊙ = 5777⁴/(274·100). Facts that let you place stars from partial data:
- `log ℒ = log L − log M` **exactly** (classical L and ℒ differ only by mass). With only
  classical L, estimate M from MIST tracks (nearest point in (log Teff, log L) among
  phases 0–6) → tier-B placement; calibrate against stars that have real log g
  (measured: −0.09 dex offset, 0.19 dex scatter) and propagate that scatter.
- Flux-weighted gravity `log g_F = log g − 4 log(Teff/10⁴ K)` (Urbaneja+2017, Kudritzki):
  `logg_from_loggF(loggF, Teff_K)` inverts it exactly.
- log ℒ, log Teff, log g are **algebraically degenerate** — never put all three in one
  regression. Check the published triple honours the identity (residual MAD ~0.003 dex;
  outliers flag inconsistent literature compilations).

## 7. Merge hygiene
Record per row: `sample`, `tier` (A = spectroscopic log g, B = track mass), `*_kind`
unit flags, and per-source systematics. Validate any merge on stars fitted by both
sources (offset, MAD, rank correlation) before combining. Keep a column for the source
of every derived quantity.
```

### A.2 kernel.py
```python
import os, re, subprocess
import numpy as np
import pandas as pd

VIZIER_TSV = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv"
VIZIER_TAP = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"
CDS_README = "https://cdsarc.cds.unistra.fr/ftp/{cat}/ReadMe"
ELL_SUN = 40649846254.21317  # (5777**4)/(274*100): spectroscopic-luminosity normalisation (cgs Teff^4/g)

def vz_tap_search(terms, table_like=None, timeout=120):
    """Keyword search of VizieR TAP_SCHEMA.tables.description. Returns [(table_name, description)]."""
    conds = [f"description LIKE '%{t}%'" for t in terms]
    if table_like:
        conds.append(f"table_name LIKE '%{table_like}%'")
    where = " OR ".join(conds) if not table_like else " AND ".join(
        ["(" + " OR ".join(conds[:-1]) + ")"] * bool(terms) + [conds[-1]])
    q = f"SELECT table_name, description FROM TAP_SCHEMA.tables WHERE {where}"
    r = subprocess.run(["curl", "-sS", "--max-time", str(timeout), "-G", VIZIER_TAP,
                        "--data-urlencode", "request=doQuery", "--data-urlencode", "lang=ADQL",
                        "--data-urlencode", "format=tsv", "--data-urlencode", f"query={q}"],
                       capture_output=True, text=True)
    out = []
    for l in r.stdout.strip().split("\n")[1:]:
        p = l.split("\t")
        if len(p) >= 2 and p[0].strip():
            out.append((p[0].strip("\"'"), p[1].strip('"')))
    return out

def vz_fetch(catalog, name, outdir="data/vizier", add_coords=False, timeout=180):
    """Download a VizieR table via asu-tsv to <outdir>/<name>.tsv and parse it."""
    os.makedirs(outdir, exist_ok=True)
    url = f"{VIZIER_TSV}?-source={catalog}&-out.max=unlimited&-out.all=2"
    if add_coords:
        url += "&-out.add=_RAJ2000,_DEJ2000"
    path = os.path.join(outdir, f"{name}.tsv")
    subprocess.run(["curl", "-sS", "--max-time", str(timeout), url, "-o", path],
                   capture_output=True, text=True)
    return vz_read_tsv(path)

def vz_read_tsv(path):
    """Parse a VizieR asu-tsv file (header / units / dashes / data). Numeric columns auto-cast."""
    lines = [l for l in open(path).read().split("\n") if not l.startswith("#")]
    dash = [i for i, l in enumerate(lines) if l and set(l.replace("\t", "")) <= set("- ") and "-" in l]
    if not dash:
        raise ValueError(f"no separator row in {path} (empty result or error page?)")
    d = dash[0]
    cols = [c.strip() for c in lines[d - 2].split("\t")]
    recs = [(l.split("\t") + [""] * len(cols))[:len(cols)] for l in lines[d + 1:] if l.strip()]
    df = pd.DataFrame(recs, columns=cols)
    for c in df.columns:
        s = pd.to_numeric(df[c].str.strip(), errors="coerce")
        df[c] = s if s.notna().sum() > 0.5 * max(len(df), 1) else df[c].str.strip()
    df.attrs["units"] = dict(zip(cols, [u.strip() for u in lines[d - 1].split("\t")]))
    return df

def vz_readme_units(catalog, labels, timeout=60):
    """Return ReadMe byte-by-byte lines mentioning any of `labels` (units are in column 3)."""
    r = subprocess.run(["curl", "-sS", "--max-time", str(timeout), CDS_README.format(cat=catalog)],
                       capture_output=True, text=True)
    return [l for l in r.stdout.split("\n") if any(k in l for k in labels)]

def star_key(name):
    """Normalise a star designation for matching across catalogs."""
    if not isinstance(name, str):
        return ""
    s = re.sub(r"\$.*?\$", "", name).upper().replace("\u2212", "-").replace("\u2013", "-")
    s = re.sub(r"^CL\*\s*", "", s)
    s = re.sub(r"\s+", "", s).replace("--", "-").replace("AZV", "AV").replace("HDE", "HD")
    s = re.sub(r"^(AV|VFTS|BI|N11|LH|HD|BD|CPD|ALS|HR|SK-\d+)0+(\d)", r"\1\2", s)
    return s

def split_aliases(s):
    """'AV14,SK9=W3' -> ['AV14','SK9','W3'] (keys). LaTeX and '=' handled."""
    s = re.sub(r"\$.*?\$", "", str(s)).replace("=", ",")
    return [star_key(a) for a in s.split(",") if a.strip()]

def alias_index(df, col):
    """Map every alias key in df[col] -> row index (first occurrence wins)."""
    idx = {}
    for i, a in df[col].items():
        for k in split_aliases(a):
            idx.setdefault(k, i)
    return idx

def spec_lum(logTeff, logg):
    """log10(L_spec/L_spec_sun) with L_spec = Teff^4/g (cgs), sun = 5777^4/(274*100)."""
    return np.log10((10.0 ** np.asarray(logTeff)) ** 4 / (10.0 ** np.asarray(logg)) / ELL_SUN)

def logg_from_loggF(loggF, Teff_K):
    """Invert flux-weighted gravity: log g = log g_F + 4 log10(Teff/1e4 K)."""
    return np.asarray(loggF) + 4 * np.log10(np.asarray(Teff_K) / 1e4)

def astropy_config_fix(base=".cache"):
    """Redirect astropy/XDG config+cache dirs into the workspace (sandbox blocks ~/.astropy)."""
    for k, v in [("XDG_CONFIG_HOME", "xdgcfg"), ("XDG_CACHE_HOME", "xdgcache"),
                 ("ASTROPY_CONFIGDIR", "astropy_cfg"), ("ASTROPY_CACHEDIR", "astropy_cache")]:
        p = os.path.abspath(os.path.join(base, v))
        os.makedirs(p, exist_ok=True)
        os.environ[k] = p

def xmatch_coords(df_a, ra_a, de_a, df_b, ra_b, de_b, radius_arcsec=2.0):
    """Nearest-neighbour sky match; returns (idx_into_b, sep_arcsec) aligned to df_a (NaN if none)."""
    astropy_config_fix()
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    ra1 = pd.to_numeric(df_a[ra_a], errors="coerce").values
    de1 = pd.to_numeric(df_a[de_a], errors="coerce").values
    ra2 = pd.to_numeric(df_b[ra_b], errors="coerce").values
    de2 = pd.to_numeric(df_b[de_b], errors="coerce").values
    gb = np.isfinite(ra2) & np.isfinite(de2)
    ga = np.isfinite(ra1) & np.isfinite(de1)
    idx = np.full(len(df_a), -1, dtype=int)
    sep = np.full(len(df_a), np.nan)
    if ga.sum() == 0 or gb.sum() == 0:
        return idx, sep
    cb = SkyCoord(ra2[gb] * u.deg, de2[gb] * u.deg)
    ca = SkyCoord(ra1[ga] * u.deg, de1[ga] * u.deg)
    i, d2d, _ = ca.match_to_catalog_sky(cb)
    ok = d2d.arcsec < radius_arcsec
    bidx = np.where(gb)[0][i]
    aidx = np.where(ga)[0]
    idx[aidx[ok]] = bidx[ok]
    sep[aidx[ok]] = d2d.arcsec[ok]
    return idx, sep
```

---

## Appendix B — Theory predictions summary (from `theory_predictions.md`)

| observable | sub-surface FeCZ prediction | core-IGW prediction | observed |
|---|---|---|---|
| α₀ vs ℒ | increases steeply (∝ F_c/F toward Eddington limit) | increases slowly with mass/L | +1.56 in log ℒ ✓FeCZ |
| ν_char vs T_eff at fixed ℒ | increases (shorter turnover in hotter, thinner zone) | decreases / no dependence | +1.53 ✓FeCZ |
| ν_char vs ℒ | weak | decreases (thicker envelope, more damping) | −0.23, ~flat |
| onset | yes, above FeCZ formation luminosity | none (core convection in all stars) | at log ℒ ≈ 3.16 ✓FeCZ |
| evolution at fixed M | amplitude grows (FeCZ deepens/strengthens) | amplitude decreases (damping grows) | α₀ ×50 ✓FeCZ |
| α₀–v_macro | correlated (shared driver) | no necessary link | ρ +0.52, partial +0.19 ✓FeCZ |
| α₀–rotation | possible weak modifier | none (Anders+2023) | weak + on all three |
| metallicity | weakens at low Z, absent below threshold | independent of Z | LMC ≈ SMC; MW confounded |
| absolute amplitude | ~observed (Cantiello+2021 within ×3 in ν) | 0.06–0.21 µmag (Anders+2023) | 64–3600 µmag ✗IGW |
