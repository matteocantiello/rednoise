#!/usr/bin/env python3
"""
Compare our refits (vandaele_fit.py) with the values digitised from Van Daele+2026 (digitize_vandaele.py), star-sector by
star-sector. Stars are placed with Bestenlehner+2025 (J/MNRAS/540/3523, table A1; matched to the Gaia DR3 IDs by position,
2 arcsec), the parameters used in their figure. Writes data_obs/vandaele_compare.csv.
"""
import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, '/mnt/home/mcantiello/work/rednoise/handoff/project_handoff/skill')
from kernel import vz_read_tsv, astropy_config_fix, xmatch_coords   # noqa: E402
RAW = '/mnt/home/mcantiello/work/rednoise/handoff/data_survey_raw/vizier'
VD = '/mnt/ceph/users/mcantiello/rednoise/vandaele2026'


def placed():
    astropy_config_fix()
    g = pd.read_csv(f'{VD}/gaia_pos.tsv', sep='\t', header=None, names=['gaia', 'ra', 'de', 'G'])
    b = vz_read_tsv(f'{RAW}/bestenlehner25_a1.tsv')
    for c in ('_RAJ2000', '_DEJ2000', 'Teff', 'logL', 'logg'):
        b[c] = pd.to_numeric(b[c], errors='coerce')
    idx, sep = xmatch_coords(g, 'ra', 'de', b, '_RAJ2000', '_DEJ2000', radius_arcsec=2.0)
    g['bidx'] = idx; g['sep'] = sep
    ok = g.bidx >= 0                     # xmatch_coords marks unmatched rows with -1
    bb = b.iloc[g.bidx[ok].astype(int)].reset_index(drop=True)
    g.loc[ok, 'BLOeM'] = bb.BLOeM.values
    g.loc[ok, 'lT'] = np.log10(bb.Teff.values * 1e3)
    g.loc[ok, 'lL'] = bb.logL.values
    g.loc[ok, 'logg'] = bb.logg.values
    g.loc[ok, 'SpType'] = bb.SpType.values
    print(f'Gaia stars {len(g)}, placed with Bestenlehner {ok.sum()}')
    return g


def main():
    g = placed()
    f = pd.read_csv(f'{HERE}/data_obs/vandaele_fits{os.environ.get("VD_TAG", "")}.csv').merge(g[['gaia', 'BLOeM', 'lT', 'lL', 'logg', 'SpType']], on='gaia')
    raw = pd.read_csv(f'{HERE}/data_obs/vandaele_digitised_raw.csv')
    rows = []
    for (metric, pw), ff in f.groupby(['metric', 'prewhitened']):
        for q, col in (('nu', 'nuchar'), ('a', 'a0cw')):
            ff = ff.assign(a0cw=ff.alpha0 / ff.Cw)
            D = raw[raw.qty == q]
            for _, r in ff[np.isfinite(ff.lT)].iterrows():
                # a value of theirs belongs to this star-sector if it sits at the star's log Teff in the Teff panel
                # and at its log L in the L panel, with the same ordinate in both
                T = D[(D.sector == r.sector) & (D.side == 'T') & (abs(D.x - r.lT) < 0.004)]
                L = D[(D.sector == r.sector) & (D.side == 'L') & (abs(D.x - r.lL) < 0.01)]
                cand = [(abs(a - b), 0.5 * (a + b)) for a in T.y for b in L.y if abs(a - b) < 0.008]
                if len(cand):
                    theirs = min(cand)[1]
                    rows.append(dict(metric=metric, prewhitened=pw, qty=q, gaia=r.gaia, sector=r.sector, lT=r.lT, lL=r.lL,
                                     ours=np.log10(r[col]), theirs=theirs, n_cand=len(cand), n_pw=r.n_pw,
                                     nuchar=r.nuchar))
    c = pd.DataFrame(rows)
    c['diff'] = c.ours - c.theirs
    c.to_csv(f'{HERE}/data_obs/vandaele_compare{os.environ.get("VD_TAG", "")}.csv', index=False, float_format='%.4f')
    s = c.groupby(['qty', 'metric', 'prewhitened']).apply(lambda x: pd.Series(dict(
        n=len(x), median_diff=x['diff'].median(), mad=np.median(abs(x['diff'] - x['diff'].median())),
        frac_within_0p05=(abs(x['diff']) < 0.05).mean(), frac_within_0p1=(abs(x['diff']) < 0.1).mean(),
        spearman=spearmanr(x.ours, x.theirs).correlation)))
    pd.set_option('display.width', 200)
    print(s.to_string(float_format=lambda v: f'{v:.3f}'))


if __name__ == '__main__':
    main()
