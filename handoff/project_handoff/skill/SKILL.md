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
