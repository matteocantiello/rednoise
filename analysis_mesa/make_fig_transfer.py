#!/usr/bin/env python3
"""Summary figure of the FeCZ -> surface transfer scenarios (run scenarios.py first)."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from test_transfer import OBSSETS, V, A, N, predict, D
from scenarios import SCEN, evaluate, NOFECZ
import transfer_models as tm
from test_transfer import sub

plt.rcParams.update({"font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8, "legend.fontsize": 6.3,
                     "xtick.labelsize": 7, "ytick.labelsize": 7})
FIG = os.environ.get('RN_FIG', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figures'))
SHOW = {'S3t M_t F_c waves, v saturated': ('#2a78d6', '-'), 'S0 MLT at the FeCZ': ('#eb6834', '--'),
        'S4 core IGW (undamped)': ('#1baf7a', ':')}
LBL = {'S3t M_t F_c waves, v saturated': r'waves $F_w=\mathcal{M}_t F_c$, $v$ saturated at $c_s$',
       'S0 MLT at the FeCZ': r'MLT: $v_{c,\max}$, $F_c/F$, $t_c$',
       'S4 core IGW (undamped)': 'core IGW (undamped)'}
scen = pd.read_csv(f'{D}/scenarios_MW_w0.0.csv').set_index('scenario')


def binned(x, y, edges, nmin=4):
    out = []
    for a, b in zip(edges[:-1], edges[1:]):
        k = (x >= a) & (x < b) & np.isfinite(y)
        if k.sum() >= nmin:
            out.append((0.5 * (a + b), np.median(y[k]), *np.percentile(y[k], [16, 84])))
    return np.array(out)


fig, axes = plt.subplots(2, 2, figsize=(7.25, 5.8))
fig.subplots_adjust(left=0.08, right=0.98, bottom=0.08, top=0.95, wspace=0.3, hspace=0.38)
eL = np.arange(2.7, 4.46, 0.15)
for ax, obs, idx, ylab, title in ((axes[0, 0], 'vmac_hot', 1, r'$\log v_{\rm macro}$ [km s$^{-1}$]', r'Macroturbulence, $\log T_{\rm eff}>4.3$'),
                                  (axes[0, 1], 'alpha0', 2, r'$\log\alpha_0$ [$\mu$mag]', r'Red-noise amplitude, $\log T_{\rm eff}>4.3$')):
    o = OBSSETS[obs]
    hot = o['x'] > 4.3
    ob = binned(o['y'][hot], o['v'][hot], eL)
    ax.errorbar(ob[:, 0], ob[:, 1], yerr=[ob[:, 1] - ob[:, 2], ob[:, 3] - ob[:, 1]], fmt='o', color='k', ms=3.5,
                lw=0.8, label='observed (MS stars, median)', zorder=5)
    for sname, (c, ls) in SHOW.items():
        key = SCEN[sname][idx - 1 if idx == 1 else 1]
        cands = V if obs.startswith('vmac') else A
        ev, sel, yhat = evaluate(obs, cands, key)
        full = np.full(len(o['v']), np.nan); full[sel] = yhat
        b = binned(o['y'][hot], full[hot], eL)
        ax.plot(b[:, 0], b[:, 1], color=c, ls=ls, lw=1.6, label=f"{LBL[sname]} ($\\Delta$BIC {ev['dBIC']:+.0f})")
    if obs == 'vmac_hot':
        cs = predict(tm.inputs(sub)['c_s'] / 1e5, False, o['x'], o['y'])
        b = binned(o['y'][hot], np.log10(cs[hot]), eL)
        ax.plot(b[:, 0], b[:, 1], color='0.5', lw=0.9, ls='-.', label='photospheric $c_s$ (model)')
    for x in (3.16, 3.70):
        ax.axvline(x, color='#e34948', ls=':', lw=0.8)
    ax.set_xlabel(r'$\log\,\mathcal{L}/\mathcal{L}_\odot$'); ax.set_ylabel(ylab); ax.set_title(title, loc='left')
    ax.legend(loc='upper left', frameon=False)
axes[0, 0].set_ylim(0.6, 2.6); axes[0, 1].set_ylim(1.0, 4.6)

# frequencies vs Teff
ax = axes[1, 0]
o = OBSSETS['nuchar']
eT = np.arange(4.2, 4.66, 0.05)
ob = binned(o['x'], o['v'], eT)
ax.errorbar(ob[:, 0], ob[:, 1], yerr=[ob[:, 1] - ob[:, 2], ob[:, 3] - ob[:, 1]], fmt='o', color='k', ms=3.5, lw=0.8,
            label=r'observed $\nu_{\rm char}$ (median)', zorder=5)
for key, c, ls, lab in (('1/(2pi t_c)', '#2a78d6', '-', r'FeCZ turnover $\times10^{a}$'),
                        ('core turnover', '#1baf7a', ':', r'core turnover $\times10^{a}$'),
                        ('nu_max scaling', '0.5', '-.', r'$\nu_{\max}\propto g\,T_{\rm eff}^{-1/2}$ $\times10^{a}$')):
    ev, sel, yhat = evaluate('nuchar', N, key)
    full = np.full(len(o['v']), np.nan); full[sel] = yhat
    b = binned(o['x'], full, eT)
    ax.plot(b[:, 0], b[:, 1], color=c, ls=ls, lw=1.6, label=f"{lab}, $a$={ev['a']:+.2f} ($\\Delta$BIC {ev['dBIC']:+.0f})")
ax.set_xlim(4.66, 4.2); ax.set_xlabel(r'$\log T_{\rm eff}$ [K]'); ax.set_ylabel(r'$\log\nu$ [d$^{-1}$]')
ax.set_title('Characteristic frequency, MS stars', loc='left'); ax.legend(loc='lower left', frameon=False)

# ratios
ax = axes[1, 1]
obsr = [(5.31, 4.91, 5.77), (6.37, 4.88, 8.62), (2.19, 1.76, 2.66)]
xs = np.arange(3)
ax.errorbar(xs - 0.24, [r[0] for r in obsr], yerr=[[r[0] - r[1] for r in obsr], [r[2] - r[0] for r in obsr]],
            fmt='o', color='k', ms=5, lw=1, label='observed')
for j, (sname, (c, ls)) in enumerate(SHOW.items()):
    r = scen.loc[sname]
    v = [r[f'ratio_{a}'] for a in ('logL', 'tau', 'logM')]
    lo = [r.get(f'ratio_{a}_lo', np.nan) for a in ('logL', 'tau', 'logM')]
    hi = [r.get(f'ratio_{a}_hi', np.nan) for a in ('logL', 'tau', 'logM')]
    ax.errorbar(xs - 0.08 + 0.14 * j, v, yerr=[np.subtract(v, lo), np.subtract(hi, v)], fmt='D', color=c, ms=4, lw=1,
                label=LBL[sname])
ax.set_xticks(xs); ax.set_xticklabels([r'along $\log\mathcal{L}$', r'along $\tau$', r'along $\log M$'])
ax.set_ylim(0, 10); ax.set_ylabel(r'$\partial\log\alpha_0/\partial\log v_{\rm macro}$')
ax.set_title('Sensitivity ratios (predictions at the observed stars)', loc='left')
ax.legend(loc='upper right', frameon=False)
for a_, l in zip(axes.flat, 'abcd'):
    a_.text(-0.13, 1.03, l, transform=a_.transAxes, fontweight='bold', fontsize=9)
for ext in ('png', 'pdf'):
    fig.savefig(f'{FIG}/fig_transfer.{ext}', dpi=200, bbox_inches='tight', facecolor='white')
print('saved', f'{FIG}/fig_transfer.png')
