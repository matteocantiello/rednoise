# Manuscript bundle — "Surface turbulence and stochastic low-frequency variability in massive stars share a sub-surface driver"

## Contents
- `main.tex` — AASTeX 6.3.1 two-column manuscript
- `aastex631.cls`, `aasjournal.bst` — class and bibliography style (from journals.aas.org)
- `refs.bib` — 37 references; all DOIs verified against CrossRef
- `tables/tab_*.tex` — deluxetables (input by main.tex)
- `tables/table_rednoise_shrd.{mrt,csv}` — MRT: 342 red-noise stars on the sHRD
- `tables/table_macroturbulence.{mrt,csv}` — MRT: 832 IACOB stars
- `figures/*.pdf` — 7 main + 3 appendix figures
- `make_figures.py` — regenerates every figure from the project artifacts
- `outline.md` — narrative outline / figure arc

## Build
```
pdflatex main && bibtex main && pdflatex main && pdflatex main
```
or `latexmk -pdf main`.

## Placeholders to fill (marked `\todo{}` in red)
1. Author list and affiliations (`\author`, `\affiliation`, `\shortauthors`)
2. Section 5.6 — MESA model comparison (currently the three target numbers only)
3. Acknowledgments
4. Check the `\facilities` list against the spectrographs actually used by the IACOB sources

## Provenance
Every number in the text was checked against the CSV tables saved as project artifacts
(`multivariate_fits_v2.csv`, `evolution_test_full.csv`, `rotation_test.csv`,
`fecz_onset_test.csv`, `hrd_gradients.csv`, `rednoise_sHRD_extended_evol.csv`,
`macroturbulence_evol.csv`) by a scripted cross-check (79 assertions).
