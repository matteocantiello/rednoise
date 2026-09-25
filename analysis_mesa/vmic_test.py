#!/usr/bin/env python3
"""
Microturbulence (Markova, Cantiello & Grassitelli 2025, A&A 701, A297; J/A+A/701/A297) against the FeCZ picture.

(1) hinge fit of log v_mic vs log L_spec with a log Teff covariate (as for v_macro), hot and all-B/O subsets;
(2) cross-match with the IACOB v_macro compilation: correlation, and partial correlation after removing log L_spec and log Teff;
(3) the MW v2 grid (non-rotating, MS + post-MS): v_mic against the model FeCZ v_c,max at each star's position (beta = 1 with a
    floor, and free beta), against the empirical L_spec + Teff plane.
Writes data_obs/vmic_test.csv and figures_v2/fig_vmic.{png,pdf}.
"""
import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from lmc_vmac_test import hinge_boot, GridMap, floor_fit, bic, vz_read_tsv, star_key, RAW

HERE = os.path.dirname(os.path.abspath(__file__))
OBS = '/mnt/home/mcantiello/work/rednoise/handoff/project_handoff/data'
OUT = f'{HERE}/data_obs'


def load():
    rows = []
    for t, puls in (('tablee1', False), ('tablee2', True)):
        m = vz_read_tsv(f'{RAW}/markova25_{t}.tsv')
        m['pulsating'] = puls
        rows.append(m)
    m = pd.concat(rows, ignore_index=True)
    m['lT'] = np.log10(m.Teff.astype(float) * 1e3)
    m['lL'] = m.logLsp.astype(float)
    m['vmic'] = m.Vmic.astype(float)
    m['key'] = m.ID.astype(str).map(star_key)
    return m


def main():
    m = load()
    hotall = m[(m.lT > 4.0)]
    print(f'Markova+25: {len(m)} stars; log Teff > 4.0: {len(hotall)} ({hotall.pulsating.sum()} pulsating); '
          f'v_mic = 0: {(hotall.vmic <= 0).sum()}')
    d = hotall[hotall.vmic > 0]
    rows = []
    for lab, s in (('log Teff > 4.0, non-pulsating', d[~d.pulsating]), ('log Teff > 4.0, all', d),
                   ('log Teff > 4.3, non-pulsating', d[(d.lT > 4.3) & ~d.pulsating])):
        f = hinge_boot(s.lL.values, s.lT.values, np.log10(s.vmic.values))
        rows.append(dict(test=f'hinge: {lab}', **f))
    # (2) v_mic vs v_macro
    V = pd.read_csv(f'{OBS}/macroturbulence_evol.csv')
    V['key'] = V.key.map(star_key)
    x = d.merge(V[['key', 'vmac']], on='key', how='inner')
    x = x[np.isfinite(x.vmac) & (x.vmac > 0)]
    lv, lm = np.log10(x.vmic), np.log10(x.vmac)
    A = np.c_[np.ones(len(x)), x.lL, x.lT]
    rv = lv - A @ np.linalg.lstsq(A, lv, rcond=None)[0]
    rm = lm - A @ np.linalg.lstsq(A, lm, rcond=None)[0]
    rows.append(dict(test='v_mic vs v_macro (IACOB cross-match)', n=len(x), spearman=spearmanr(lv, lm).correlation,
                     partial_spearman=spearmanr(rv, rm).correlation, slope=np.polyfit(lm, lv, 1)[0]))
    # (3) model FeCZ velocity
    gm = GridMap('MW')
    s = d[~d.pulsating]
    X = gm(s.lT.values, s.lL.values)
    sel = np.isfinite(X)
    y = np.log10(s.vmic.values[sel])
    B = np.c_[np.ones(sel.sum()), s.lL.values[sel], s.lT.values[sel]]
    bres = y - B @ np.linalg.lstsq(B, y, rcond=None)[0]
    yhat, p = floor_fit(X[sel], y)
    pos = X[sel] > 0
    cfree = np.polyfit(np.log10(X[sel][pos]), y[pos], 1)
    rows.append(dict(test='model v_c,max (MW v2) -> v_mic, beta = 1 with floor', n=int(sel.sum()),
                     dBIC=bic(y - yhat, 2) - bic(bres, 3), floor_kms=10 ** p[0], norm=10 ** p[1],
                     rms=np.std(y - yhat), rms_plane=np.std(bres),
                     spearman=spearmanr(X[sel][pos], y[pos]).correlation, free_slope=cfree[0],
                     n_no_FeCZ=int((~pos).sum())))
    res = pd.DataFrame(rows)
    res.to_csv(f'{OUT}/vmic_test.csv', index=False, float_format='%.3f')
    pd.set_option('display.width', 220)
    print(res.to_string(index=False, float_format=lambda v: f'{v:.2f}'))

    # figure: v_mic, v_macro and model v_c,max vs L_spec for hot stars
    fig, ax = plt.subplots(1, 2, figsize=(7.25, 3.0))
    fig.subplots_adjust(wspace=0.3)
    edges = np.arange(1.0, 4.6, 0.2)
    for s_, col, lab in ((d[~d.pulsating], '#2a78d6', r'$v_{\rm mic}$ (Markova+25, non-pulsating)'),
                         (V[(V.logTeff_sp > 4.0) & (V.vmac > 0)].rename(columns={'logL_sp': 'lL', 'vmac': 'v'}), '#eb6834',
                          r'$v_{\rm macro}$ (IACOB)')):
        vv = s_['vmic'] if 'vmic' in s_ else s_['v']
        b = []
        for a, c in zip(edges[:-1], edges[1:]):
            k = s_.lL.between(a, c)
            if k.sum() >= 5:
                b.append((0.5 * (a + c), *np.percentile(vv[k], [50, 16, 84])))
        b = np.array(b)
        ax[0].errorbar(b[:, 0], b[:, 1], yerr=[b[:, 1] - b[:, 2], b[:, 3] - b[:, 1]], fmt='o-', color=col, ms=3.5,
                       lw=1, label=lab)
    Xa = gm(d.lT.values, d.lL.values)
    b = []
    for a, c in zip(edges[:-1], edges[1:]):
        k = (d.lL.values >= a) & (d.lL.values < c) & np.isfinite(Xa) & (Xa > 0)
        if k.sum() >= 5:
            b.append((0.5 * (a + c), np.median(Xa[k])))
    b = np.array(b)
    ax[0].plot(b[:, 0], b[:, 1], 'k--', lw=1.2, label=r'model FeCZ $v_{c,\max}$ at the $v_{\rm mic}$ stars')
    ax[0].axvline(2.96, color='0.5', ls=':', lw=0.8)
    ax[0].set_yscale('log'); ax[0].set_xlabel(r'$\log\,\mathcal{L}/\mathcal{L}_\odot$')
    ax[0].set_ylabel(r'velocity [km s$^{-1}$]'); ax[0].set_title(r'$\log T_{\rm eff}>4.0$ (medians)', loc='left')
    ax[0].legend(frameon=False, fontsize=6)
    ax[1].scatter(x.vmac, x.vmic, s=8, c=x.lL, cmap='viridis', lw=0)
    ax[1].set_xscale('log'); ax[1].set_yscale('log')
    ax[1].set_xlabel(r'$v_{\rm macro}$ [km s$^{-1}$]'); ax[1].set_ylabel(r'$v_{\rm mic}$ [km s$^{-1}$]')
    r = res[res.test.str.startswith('v_mic vs v_macro')].iloc[0]
    ax[1].set_title(f"cross-match N = {int(r.n)}: rho = {r.spearman:+.2f}, partial {r.partial_spearman:+.2f}",
                    loc='left', fontsize=7.5)
    for ext in ('png', 'pdf'):
        fig.savefig(f'{HERE}/figures_v2/fig_vmic.{ext}', dpi=200, bbox_inches='tight', facecolor='white')
    print('saved figures_v2/fig_vmic.png')


if __name__ == '__main__':
    main()
