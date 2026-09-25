#!/usr/bin/env python3
"""
Download SPOC 2-min TESS light curves of the Galactic red-noise stars (Bowman+2020 by TIC; Shen+2024 by coordinates), for a
homogeneous refit with the pipeline used on the Van Daele SMC light curves (vandaele_fit.py).

Per star at most NMAX sectors: for Shen 2024 their fitted sectors first (J/ApJS/275/2 table2 'Sec'), for Bowman 2020
sectors <= 13 first (their Cycle-1 data), then the earliest others. Each sector is written as a text file
'<name>_sector<NN>_SPOC.txt' with columns time [BTJD], PDCSAP flux, SAP flux [e-/s] (quality-flagged cadences removed).
Run with the lightkurve venv:  ../.venv_tess/bin/python tess_mw_download.py
"""
import os
import sys
import re
import warnings
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

warnings.filterwarnings('ignore')
import lightkurve as lk                 # noqa: E402

sys.path.insert(0, '/mnt/home/mcantiello/work/rednoise/handoff/project_handoff/skill')
from kernel import vz_read_tsv          # noqa: E402

OUT = '/mnt/ceph/users/mcantiello/rednoise/tess_mw/SPOC'
RAW = '/mnt/home/mcantiello/work/rednoise/handoff/data_survey_raw/vizier'
B20 = '/mnt/home/mcantiello/work/rednoise/handoff/project_handoff/data/bowman2020_A36_tablea2.tsv'
NMAX = 4
os.makedirs(OUT, exist_ok=True)


def targets():
    b = pd.read_csv(B20, sep='\t', comment='#')
    b = b[pd.to_numeric(b.TIC, errors='coerce').notna()]
    rows = [dict(name=n.strip(), query=f'TIC {int(t)}', pref=list(range(1, 14)), sample='Bowman2020')
            for n, t in zip(b.Name, b.TIC.astype(int))]
    t1 = vz_read_tsv(f'{RAW}/shen24_t1.tsv'); t2 = vz_read_tsv(f'{RAW}/shen24_t2.tsv')
    sec = t2.groupby(t2.Name.str.strip()).Sec.apply(lambda s: sorted(pd.to_numeric(s, errors='coerce').dropna().astype(int)))
    for _, r in t1.iterrows():
        n = r.Name.strip()
        rows.append(dict(name=n, query=f'{float(r._RA):.5f} {float(r._DE):.5f}', pref=sec.get(n, []), sample='Shen2024'))
    d = pd.DataFrame(rows).drop_duplicates('name')
    return d


def plain(a):
    """astropy Masked / numpy masked -> float array with NaN where masked."""
    if hasattr(a, 'unmasked'):
        return np.where(np.asarray(a.mask), np.nan, np.asarray(a.unmasked, float))
    return np.ma.filled(np.ma.asarray(a, dtype=float), np.nan)


def fname(name, sector):
    return f"{OUT}/{re.sub(r'[^A-Za-z0-9+-]', '_', name)}_sector{sector:02d}_SPOC.txt"


def one(r):
    try:
        return _one(r)
    except Exception as e:
        return dict(name=r['name'], sample=r['sample'], status=f'error: {e}'[:120], n=0)


def _one(r):
    try:
        sr = lk.search_lightcurve(r['query'], author='SPOC', exptime=120, radius=5 if ' ' in r['query'] and 'TIC' not in r['query'] else None)
    except Exception as e:
        return dict(name=r['name'], status=f'search failed: {e}'[:120], n=0)
    if len(sr) == 0:
        return dict(name=r['name'], status='no SPOC 2-min data', n=0)
    tab = sr.table.to_pandas()
    tics = tab.target_name.astype(str).unique()
    if len(tics) > 1:                        # coordinate search hit several targets: keep the nearest
        tab = tab[tab.target_name.astype(str) == str(tab.sort_values('distance').target_name.iloc[0])]
    secs = [int(re.search(r'(\d+)', str(m)).group(1)) for m in tab.mission]
    order = sorted(range(len(secs)), key=lambda i: (secs[i] not in r['pref'], secs[i]))[:NMAX]
    got = []
    for i in order:
        f = fname(r['name'], secs[i])
        if os.path.exists(f):
            got.append(secs[i]); continue
        try:
            lc = sr[int(tab.index[i])].download(quality_bitmask='default', flux_column='pdcsap_flux',
                                                 download_dir='/mnt/ceph/users/mcantiello/rednoise/tess_mw/cache')
            t = np.asarray(lc.time.value, float); p = plain(lc.flux.value); s = plain(lc.sap_flux.value)
            ok = np.isfinite(t) & np.isfinite(p) & np.isfinite(s)
            np.savetxt(f, np.c_[t[ok], p[ok], s[ok]], fmt='%.7f %.4f %.4f')
            got.append(secs[i])
        except Exception as e:
            print(r['name'], secs[i], 'download failed', str(e)[:100], flush=True)
    return dict(name=r['name'], sample=r['sample'], tic=tab.target_name.iloc[0], status='ok', n=len(got),
                sectors=' '.join(map(str, got)), pref_hit=sum(s in r['pref'] for s in got))


def main():
    d = targets()
    print(f'{len(d)} targets', flush=True)
    with ProcessPoolExecutor(4) as ex:
        res = list(ex.map(one, d.to_dict('records')))
    res = pd.DataFrame(res)
    res.to_csv('/mnt/ceph/users/mcantiello/rednoise/tess_mw/download_log.csv', index=False)
    print(res.status.value_counts().to_string()); print('light curves:', res.n.sum())


if __name__ == '__main__':
    main()
