# STATUS: where the project stands (restart point for a new session)

*Last updated 2026-09-25 17:30 EDT. Read this first, then `models/GRID_LOG.md` (grids) and the `analysis_mesa/REPORT_*.md`
files named below. Standing rules: org CLAUDE.md LAWs (Slurm: read-only commands only, no polling; scoped searches); log all
grid work in `models/GRID_LOG.md`; never rebuild or modify `models/template_new`; commit/push only when the user asks
(message ends with the Co-Authored-By line); never add a reference without a resolved DOI.*

## 1. Running right now (check, don't poll)
| what | where | state at 16:45 on 2026-09-25 |
|---|---|---|
| v2 grids (Slurm, user-submitted) | `models/grids_v2/{MW,LMC,SMC}/w*` | MW and LMC: 24–29 of 34 masses done per sub-grid, high masses still running. **SMC: w0.0 and w0.2 running, w0.4 and w0.6 pending** (one `squeue -u $USER` at ~15:00). Status file `models/grids_v2/grid_status.txt` (last written 10:11; regenerate with the grid scripts). |
| composition bracket (local, not Slurm) | `/home/mcantiello/rednoise_tests/zbracket_Z0{10,18}/MW/w0.0/M{10..30}` | 12 runs via `run_zbracket.sh` since 10:59; at models 1150–1270; 24 h timeout. MS only (stop X_c = 1e-3). See GRID_LOG entry "composition bracket". |

## 2. Done on 2026-09-25 (after commit f555c5f; all uncommitted)
| topic | report | code → outputs | one-line result |
|---|---|---|---|
| Magnetic stars | `analysis_mesa/REPORT_magnetic.md` | `magnetic_test.py` → `data_obs/magnetic_{test,stars}.csv` | 30 Shen 2023 magnetic OB stars: α₀ +0.22 dex *above* the non-magnetic plane; no B_d trend; only 3–6 stars above model B_shutoff. No suppression. Two literature refs flagged "to verify". |
| GP ↔ Lorentzian ν_char | `REPORT_nuchar_Z.md` | `digitize_bdw22.py` → `data_obs/bdw22_nuchar.csv`; `nuchar_Z_test.py` | BDW22 figures digitised (validated to 0.003 dex): GP − Lor = −0.02 dex median, 0.19 scatter. Bowman 2024 LMC/SMC at matched position −0.03 to −0.12 dex. Paper §4.6 "method-blocked" claim is obsolete. |
| Van Daele 2026 replication | `REPORT_vandaele.md` | `vandaele_fit.py` (env-configurable), `digitize_vandaele.py`, `vandaele_compare.py` | Not reproducible star by star (fit-metric choice alone: 0.17 dex); ensemble slopes only partly reproduced (table corrected after the placement bug). |
| Homogeneous MW vs SMC ν_char | `REPORT_mw_smc.md` | `tess_mw_download.py` (lightkurve venv), `mw_smc_nuchar.py`, `smc_instrumental_check.py` | 171 Galactic stars / 532 SPOC sectors; our pipeline matches Shen 2024 per sector (−0.02 dex). Raw SMC deficit −0.15 to −0.33 dex (30 stars). PSF floor is not separable from SLF at the periodogram level. |
| Noise-floor model | `REPORT_floor.md` | `smc_floor_model.py` (FLOOR_STEP=calib/fits/inject), `floor_analysis.py` (+ `expected_R()`, `expected_detection()`) → `data_obs/floor_*.csv` | Whittle fits with an empirical floor from the 38 non-variable stars; 3328 injections. Only 24% of SLF-class sectors are above the floor (17% false positives). Galactic-like SLF is disfavoured (expected naive log ν −0.35 vs −0.69 observed): SMC SLF is ≲ 0.1× Galactic amplitude or ≳ 0.4 dex slower. Model-dependent. |
| Van Daele correspondence | `Vandaele_correspondence/README.md` | — | Email thread logged (timeline, v1 → published changes, action items, paper implications). **Private: keep out of public repos** (not yet in `.gitignore`; user to decide). |
| Email to Van Daele | `Vandaele_correspondence/email_draft_2026-09-25.md` | `make_fig_vandaele_email.py` → `fig{1,2,3}_*.png` | **Draft, not sent.** User to review tone and whether to include the forward-model result; could be placed into a Gmail draft on request. |

**Bug fixed today (affected earlier numbers, now corrected everywhere):** `kernel.xmatch_coords` returns −1 for unmatched
rows. `vandaele_compare.placed()` had treated −1 as a match, which gave 37 SMC stars the last Bestenlehner row. Only 54 of 91
BLOeM PSF stars have Bestenlehner+2025 parameters; the rest (mostly cool supergiants, Patrick+2025) are unplaced.

Earlier on 2026-09-25 (committed in f555c5f and before): LMC v_macro shift (`REPORT_lmc_vmac.md`), PB25 calibration and
α₀-onset robustness (`REPORT_pb25.md`), v_mic (`REPORT_vmic.md`), rotation and LMC regime shift (`REPORT_rotation_Z.md`,
`REPORT_regime.md`), data survey (`handoff/DATA_SURVEY_2026-09-25.md`), MESA upgrade handoff (`handoff/MESA_UPGRADE_HANDOFF.md`).

## 2b. Reusable skills created from this project (2026-09-25)
User-level, auto-loaded in every project from `~/.claude/skills/`:
- `tess-slf-fitting` (SKILL.md + `slf_fit.py`: Whittle fits, prewhitening, floor model, injection, lightkurve notes);
- `figure-digitize` (SKILL.md + `digitize.py`: tick calibration, colour markers, marker lines, validation);
- `citation-verify` (SKILL.md + `verify_refs.py`: CrossRef/arXiv audit of a .bib).
Update draft for the claude.ai-synced `vizier-catalog-harmonize` (−1 sentinel fix, `take_matched`, arXiv-source tables):
`~/.claude/skill_drafts/vizier-catalog-harmonize/` (see README_UPLOAD.md; the user uploads it). The project's own
`handoff/project_handoff/skill/kernel.py` docstring is fixed.

## 3. Waiting on the user
1. Review and send the Van Daele email (and optionally cc Dominic Bowman).
2. Whether to add `Vandaele_correspondence/` to `.gitignore` (recommended: private emails and a referee-stage draft).
3. Commit the day's work (nothing since f555c5f is committed; `.gitignore` already gained `.venv_tess/`).
4. Earlier open decision: paper revisions below are not yet made. The user decided on 2026-09-25 to migrate §5.6 to v2 once MW
   v2 is complete.

## 4. Next steps (in order)
1. **When the SMC v2 grid finishes:** rerun the extraction and regime chain
   (`RN_GRID=../models/grids_v2 RN_DATA=data_v2 RN_FIG=figures_v2 python3 extract_mesa.py SMC/w0.0`, then
   `regime_crossings.py`). Then:
   - (a) SMC FeCZ-present map with the Van Daele (SLF / no SLF), full BLOeM and Bowman 2024 stars overlaid, reporting the
     fraction in the no-FeCZ window (correspondence action item 1; a paper figure);
   - (b) model Δlog ν_c (SMC − MW) at the SMC stars, to compare with the floor-model constraint;
   - (c) FeCZ vs He II/He I CZ at each Van Daele star, to test whether their blue/green F-test split matches a model zone switch.
2. **Bowman & Van Daele 2024 window stars** (correspondence action item 2): which B24 SMC/LMC stars sit in the no-FeCZ window;
   distance from the edge, α₀, detection strength (the user recalls most are marginal, e.g. TIC 181887485).
3. **Floor model, decisive version:** measure the floor independently. Options: LEMONS/PSF light curves of non-target or
   background sources in the same FFIs (needs FFI cutouts, e.g. TESScut via astroquery in `.venv_tess`), or FFI photometry of
   fainter Galactic OB stars made the same way. Also try γ ≠ 2 injections, and add the 37 unplaced stars with Patrick+2025
   parameters.
4. **Composition bracket:** when the 12 local runs finish, run `regime_crossings.py` and `scenarios.py` with
   `RN_GRID=/home/mcantiello/rednoise_tests/zbracket_Z010` (and Z018); log the result in GRID_LOG.
5. **Paper revisions (`paper_apj/main.tex`), none done yet:**
   - abstract, §4.2 and §5.1: only v_macro has a robust onset; the α₀ onset is not robust (REPORT_pb25);
   - §5.3: the 30 µmag floor is unsupported;
   - §4.6 and limitation (1): BDW22 calibration and the matched-position metallicity results (REPORT_nuchar_Z, REPORT_mw_smc,
     REPORT_floor);
   - new text discussing Van Daele+2026 and Bowman & Van Daele 2024: the selection-bias argument, the PSF floor, and a
     collegial tone (see correspondence README);
   - a magnetic-star paragraph (REPORT_magnetic);
   - v_mic (REPORT_vmic) and the LMC v_macro shift (REPORT_lmc_vmac);
   - §5.6 and the appendix migrated to v2 once MW v2 is complete;
   - **bibliography errors found by `citation-verify` (refs.bib, not yet fixed):**
     - `bestenlehner2025`: the title is paper XII, but DOI 10.1051/0004-6361/202452491 (A&A 695, A198) is paper **XI**,
       "Pipeline-determined physical properties of Magellanic Cloud OB stars". Fix the title, and check which paper the text
       means.
     - `deburgos2024`: the title says "IACOB XII. New grid-based automatic tool…", but the DOI is "IACOB **X**. Large-scale
       quantitative spectroscopic analysis of Galactic luminous…". Decide which paper is cited and fix the entry.
     - Audit CSV: rerun `python3 ~/.claude/skills/citation-verify/verify_refs.py paper_apj/refs.bib --out paper_apj/refs_check.csv`.

## 5. Environment and gotchas
- Python with pandas: `bash -lc 'module load python; python3 …'`; system python3 lacks pandas.
- lightkurve/astroquery: `/mnt/home/mcantiello/work/rednoise/.venv_tess/bin/python` (venv on module python with
  system site-packages; only lightkurve was pip-installed). Needed only for `tess_mw_download.py`.
- Data on Ceph: `/mnt/ceph/users/mcantiello/rednoise/vandaele2026/` (Zenodo 20540863 light curves, `gaia_pos.tsv`, logs) and
  `/mnt/ceph/users/mcantiello/rednoise/tess_mw/` (`SPOC/*.txt` = time, PDCSAP, SAP; `download_log.csv`; lightkurve cache).
- `vandaele_fit.py` settings via env: `VD_LC`, `VD_GLOB`, `VD_COL` (1 PDCSAP, 2 SAP), `VD_TAG` (output suffix),
  `VD_MATCH_A0CW` (noise-degrade), `VD_FMIN`/`VD_FMAX` (default 0.1/40), `VD_FLUX` (1 = e⁻/s). Outputs:
  `data_obs/vandaele_fits{TAG}.csv` (SMC without a tag; `_mw`, `_mwdeg`, `_mwsap` Galactic;
  `_fullrange_mmag` the first SMC run).
- Single-sector ν_char depends on fit choices at the 0.1–0.2 dex level. Always compare galaxies with one pipeline and state
  the variant. The Whittle fits (`smc_floor_model.py`) avoid the linear-vs-log ambiguity.
- LaTeX: `module load texlive`. Paper backups: `main.tex.pre_mesa`, `.pre_lit`.
- MESA r26.04.1 at `/mnt/home/mcantiello/mesa-26.04.1`, SDK 26.6.1; see GRID_LOG for environment lines.
