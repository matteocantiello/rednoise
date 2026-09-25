#!/usr/bin/env python3
"""
LMC vs MW macroturbulence: (1) method offset between Serebriakova's and the IACOB v_macro scales, estimated from
Galactic stars at matched (log L_spec, log Teff) below the onset, where no metallicity effect is expected in either
sample; (2) the luminosity shift Delta of the LMC relative to the Galactic relation, fitted as
    log v_LMC(L) = log v_MW(L - Delta) + log k
with k a free scale (absorbs any method offset) and v_MW a smooth fit to the Galactic running median; 16-84% from
500 bootstrap resamples of both samples. The v2 grids predict Delta = +0.60 (onset proxy, v_c = 3 km/s) and +0.24
(saturation proxy, Gamma_Fe = 0.8) (REPORT_rotation_Z.md section 4).
Writes $RN_DATA/lmc_vmac_shift.csv.
"""
import os
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from lmc_vmac_test import load_mw, load_lmc, vz_read_tsv, RAW, D2

rng = np.random.default_rng(11)
TLO, THI = 4.0, 4.53          # common log Teff range of the two samples (LMC hot end ~4.53)


def serebriakova_mw():
    s = vz_read_tsv(f'{RAW}/serebriakova24.tsv')
    ra, de = pd.to_numeric(s._RAJ2000, errors='coerce'), pd.to_numeric(s._DEJ2000, errors='coerce')
    mw = ~((ra > 60) & (ra < 100) & (de < -63) & (de > -75))
    vm = pd.to_numeric(s.Vmacro, errors='coerce')
    ok = mw & np.isfinite(vm) & (vm > 0)
    return pd.DataFrame(dict(lT=np.log10(s.Teff[ok].astype(float)), lL=s.Ls[ok].astype(float), vmac=vm[ok]))


def running_median(x, y, edges):
    c, m = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        k = (x >= a) & (x < b)
        if k.sum() >= 8:
            c.append(0.5 * (a + b)); m.append(np.median(y[k]))
    return np.array(c), np.array(m)


def fit_shift(mw, lmc, lmin=3.4, fix_k=None):
    edges = np.arange(1.8, 4.7, 0.1)
    cx, cy = running_median(mw.lL.values, np.log10(mw.vmac.values), edges)
    f = lambda L: np.interp(L, cx, cy, left=np.nan, right=np.nan)
    d = lmc[lmc.lL >= lmin]
    x, y = d.lL.values, np.log10(d.vmac.values)

    def res(p):
        k = p[1] if fix_k is None else fix_k
        r = y - (f(x - p[0]) + k)
        return np.where(np.isfinite(r), r, 1.0)          # penalise shifts that leave the MW relation's range
    best = min((least_squares(res, [d0, 0.0]) for d0 in np.arange(-0.4, 1.01, 0.1)), key=lambda r: np.sum(r.fun ** 2))
    x_ = best.x.copy()
    if fix_k is not None:
        x_[1] = fix_k
    return x_, len(d), np.sqrt(np.mean(best.fun ** 2))


def main():
    mw_i, lmc = load_mw(), load_lmc()
    mw_i = mw_i[mw_i.lT.between(TLO, THI)]
    lmc_h = lmc[lmc.lT.between(TLO, THI)]
    rows = []
    # (1) method offset below the onset: Serebriakova Galactic vs IACOB, log L 2.0-2.8, same Teff range
    sm = serebriakova_mw()
    for lo, hi in ((2.0, 2.8), (2.2, 2.8)):
        a = sm[sm.lL.between(lo, hi) & sm.lT.between(TLO, THI)]
        b = mw_i[mw_i.lL.between(lo, hi)]
        if len(a) >= 3 and len(b) >= 3:
            off = np.median(np.log10(a.vmac)) - np.median(np.log10(b.vmac))
            bs = [np.median(np.log10(rng.choice(a.vmac, len(a)))) - np.median(np.log10(rng.choice(b.vmac, len(b))))
                  for _ in range(2000)]
            rows.append(dict(test=f'method offset (Serebriakova MW - IACOB), log L {lo}-{hi}', value=off,
                             lo=np.percentile(bs, 16), hi=np.percentile(bs, 84), n=f'{len(a)} vs {len(b)}'))
    # (2) shift of the LMC relation
    for lmin in (3.4, 3.6):
        (dl, lk), n, rms = fit_shift(mw_i, lmc_h, lmin)
        bs = []
        for _ in range(500):
            m1 = mw_i.iloc[rng.integers(0, len(mw_i), len(mw_i))]
            l1 = lmc_h.iloc[rng.integers(0, len(lmc_h), len(lmc_h))]
            bs.append(fit_shift(m1, l1, lmin)[0])
        bs = np.array(bs)
        rows.append(dict(test=f'LMC shift Delta log L (LMC stars log L >= {lmin})', value=dl,
                         lo=np.percentile(bs[:, 0], 16), hi=np.percentile(bs[:, 0], 84), n=n, rms=rms))
        rows.append(dict(test=f'LMC scale log k (LMC stars log L >= {lmin})', value=lk,
                         lo=np.percentile(bs[:, 1], 16), hi=np.percentile(bs[:, 1], 84), n=n))
    # (2b) shift with the method offset fixed at 0 (supported by (1)); full bootstrap distribution of Delta
    (dl0, _), n0, rms0 = fit_shift(mw_i, lmc_h, 3.4, fix_k=0.0)
    bs0 = np.array([fit_shift(mw_i.iloc[rng.integers(0, len(mw_i), len(mw_i))],
                              lmc_h.iloc[rng.integers(0, len(lmc_h), len(lmc_h))], 3.4, fix_k=0.0)[0][0] for _ in range(500)])
    rows.append(dict(test='LMC shift Delta log L, method offset fixed at 0 (log L >= 3.4)', value=dl0,
                     lo=np.percentile(bs0, 16), hi=np.percentile(bs0, 84), n=n0, rms=rms0))
    for lab, arr in (('free k', bs[:, 0]), ('k = 0', bs0)):
        print(f'bootstrap Delta ({lab}): P(<0.2) {np.mean(arr < 0.2):.2f}, P(0.2-0.7) {np.mean((arr >= 0.2) & (arr <= 0.7)):.2f},'
              f' P(>0.7) {np.mean(arr > 0.7):.2f}; median {np.median(arr):+.2f}')
    # (3) plain ratio at matched L, for reference (no shift)
    for lo, hi in ((3.6, 4.0), (4.0, 4.4)):
        a, b = lmc_h[lmc_h.lL.between(lo, hi)], mw_i[mw_i.lL.between(lo, hi)]
        rows.append(dict(test=f'median log v_LMC - log v_MW at log L {lo}-{hi} (no shift)',
                         value=np.median(np.log10(a.vmac)) - np.median(np.log10(b.vmac)), n=f'{len(a)} vs {len(b)}'))
    res = pd.DataFrame(rows)
    res.to_csv(f'{D2}/lmc_vmac_shift.csv', index=False, float_format='%.3f')
    pd.set_option('display.width', 200)
    print(res.to_string(index=False, float_format=lambda v: f'{v:+.2f}'))


if __name__ == '__main__':
    main()


def figure():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 8, "legend.fontsize": 6.5})
    mw_i, lmc = load_mw(), load_lmc()
    mw_i, lmc_h = mw_i[mw_i.lT.between(TLO, THI)], lmc[lmc.lT.between(TLO, THI)]
    r = pd.read_csv(f'{D2}/lmc_vmac_shift.csv').set_index('test')
    d0 = r.loc['LMC shift Delta log L, method offset fixed at 0 (log L >= 3.4)']
    edges = np.arange(1.8, 4.7, 0.1)
    cx, cy = running_median(mw_i.lL.values, np.log10(mw_i.vmac.values), edges)
    fig, ax = plt.subplots(figsize=(3.6, 3.0))
    ax.scatter(mw_i.lL, mw_i.vmac, s=4, color='#2a78d6', alpha=0.2, lw=0)
    ax.scatter(lmc_h.lL, lmc_h.vmac, s=9, color='#eb6834', alpha=0.8, lw=0, label=f'LMC, Serebriakova ({len(lmc_h)})')
    ax.plot(cx, 10 ** cy, color='#2a78d6', lw=1.8, label=f'MW running median (IACOB + Holgado, {len(mw_i)})')
    ax.fill_betweenx([1, 300], 0, 0)  # keep limits
    for dl, ls, lab in ((0.24, ':', 'MW shifted +0.24 (model saturation shift)'),
                        (0.60, '-.', 'MW shifted +0.60 (model onset shift)'),
                        (d0.value, '--', f'MW shifted {d0.value:+.2f} (fit, no method offset)')):
        ax.plot(cx + dl, 10 ** cy, color='0.25', ls=ls, lw=1.1, label=lab)
    ax.set_yscale('log'); ax.set_ylim(3, 250); ax.set_xlim(1.8, 4.8)
    ax.set_xlabel(r'$\log\,\mathcal{L}/\mathcal{L}_\odot$'); ax.set_ylabel(r'$v_{\rm macro}$ [km s$^{-1}$]')
    ax.set_title(r'Macroturbulence, $4.0<\log T_{\rm eff}<4.53$', loc='left', fontsize=8)
    ax.legend(loc='lower right', frameon=False)
    fig.savefig(f'{os.path.dirname(D2)}/figures_v2/fig_lmc_vmac_shift.png', dpi=200, bbox_inches='tight', facecolor='white')
    fig.savefig(f'{os.path.dirname(D2)}/figures_v2/fig_lmc_vmac_shift.pdf', bbox_inches='tight', facecolor='white')
    print('saved figures_v2/fig_lmc_vmac_shift.png')


if __name__ == '__main__' and os.environ.get('RN_FIG_ONLY'):
    figure()
