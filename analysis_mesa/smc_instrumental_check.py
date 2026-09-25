#!/usr/bin/env python3
"""
Is the low SMC nu_char (mw_smc_nuchar.py) instrumental? Checks on the Van Daele PSF light curves with our pipeline:
(a) stars Van Daele class as having no significant variability (blank class in their Table A1): what do our fits return?
    If they pass the cuts with nu_char like the SLF stars, the PSF noise floor itself looks like 'SLF'.
(b) SLF-star residuals from the Galactic plane (PDCSAP and SAP) against TESS magnitude: instrumental red noise grows
    for fainter stars.
(c) alpha0/C_w of our SMC fits by class: the red-noise excess over white noise of 'constant' stars.
Writes data_obs/smc_instrumental_check.csv.
"""
import os
import re
import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mw_smc_nuchar as ms          # noqa: E402
TEX = '/mnt/home/mcantiello/work/rednoise/handoff/data_survey_raw/vandaele2026/src/paper_bloemphot.tex'


def classes():
    rows = []
    for l in open(TEX):
        m = re.match(r'^(\d-\d{3}) & (\d+) & ([\d.]+) & ([^&]*)&(.*)', l)
        if m:
            c = re.sub(r'%.*', '', m.group(5)).replace('\\\\', '').strip()
            cls = 'SLF' if c.startswith('SLF') else ('none' if c == '' else ('pulsator' if 'pulsator' in c else 'binary/rot'))
            rows.append(dict(BLOeM=m.group(1), gaia=int(m.group(2)), Tmag=float(m.group(3)), SpT=m.group(4).strip(), cls=cls))
    return pd.DataFrame(rows)


def main():
    c = classes()
    print(c.cls.value_counts().to_dict())
    f = pd.read_csv(f'{HERE}/data_obs/vandaele_fits.csv').merge(c, on='gaia')
    rows = []
    for metric in ('lin', 'log'):
        g = f[(f.metric == metric) & f.prewhitened]
        ok = (g.nuchar >= 0.15) & (g.dBIC_vs_white > 0)
        for cl, s in g.groupby('cls'):
            k = ok[s.index]
            rows.append(dict(check='(a,c) our fits by Van Daele class', metric=metric, cls=cl, n_sectors=len(s),
                             frac_pass=k.mean(), median_lognu_pass=np.log10(s.nuchar[k]).median(),
                             median_log_a0cw=np.log10(s.alpha0 / s.Cw)[k].median(), median_Tmag=s.Tmag.median()))
    # (b) residuals vs Tmag, using the plane from each MW variant
    E = pd.read_csv(f'{HERE}/data_obs/vandaele_fits.csv')
    for tag in ('_mw', '_mwsap'):
        fmw = pd.read_csv(f'{HERE}/data_obs/vandaele_fits{tag}.csv')
        for metric in ('lin', 'log'):
            mw = ms.place_mw(ms.per_star(fmw, metric)); mw = mw[np.isfinite(mw.lT) & np.isfinite(mw.lL)]
            A = np.c_[np.ones(len(mw)), mw.lL, mw.lT]; cc = np.linalg.lstsq(A, mw.lnu, rcond=None)[0]
            smc = ms.place_smc(ms.per_star(E, metric))
            smc['gaia'] = smc.id.str.replace('GAIA DR3 ', '').astype(int)
            smc = smc.merge(c, on='gaia')
            smc = smc[np.isfinite(smc.lT) & np.isfinite(smc.lL)]
            smc['res'] = smc.lnu - np.c_[np.ones(len(smc)), smc.lL, smc.lT] @ cc
            for cl, s in (('SLF', smc[smc.cls == 'SLF']), ('none', smc[smc.cls == 'none']), ('all', smc)):
                if len(s) < 4:
                    continue
                br = s[s.Tmag < 12.8]
                rows.append(dict(check=f'(b) residual vs Tmag, MW{tag}', metric=metric, cls=cl, n_sectors=len(s),
                                 median_residual=s.res.median(), spearman_res_Tmag=spearmanr(s.Tmag, s.res).correlation,
                                 spearman_p=spearmanr(s.Tmag, s.res).pvalue, median_residual_Tmag_lt_12p8=br.res.median(),
                                 n_bright=len(br)))
    r = pd.DataFrame(rows)
    r.to_csv(f'{HERE}/data_obs/smc_instrumental_check.csv', index=False, float_format='%.3f')
    pd.set_option('display.width', 250)
    print(r.to_string(index=False, float_format=lambda v: f'{v:+.2f}'))


if __name__ == '__main__':
    main()
