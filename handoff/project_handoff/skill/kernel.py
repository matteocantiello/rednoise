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
    """Nearest-neighbour sky match; returns (idx_into_b, sep_arcsec) aligned to df_a.
    Unmatched rows get idx = -1 (NOT NaN); guard with idx >= 0 before any .iloc (iloc[-1] = last row)."""
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
