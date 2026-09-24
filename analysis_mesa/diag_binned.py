"""Binned (hot, MS) observed v_macro / alpha0 / nu_char vs log L_spec with candidate predictors."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
import transfer_models as tm
from test_transfer import OBSSETS, predict, V, A, N
plt.rcParams.update({"font.size": 8})
edges = np.arange(2.7, 4.45, 0.15)
def binned(x, y):
    out = []
    for a, b in zip(edges[:-1], edges[1:]):
        k = (x >= a) & (x < b) & np.isfinite(y)
        if k.sum() >= 4: out.append((0.5*(a+b), np.median(y[k]), *np.percentile(y[k], [16, 84])))
    return np.array(out)
fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
sets = [('vmac_hot', V, ['v_c,max', 'Mach x c_s', 'core IGW, v_g=c_s', 'C09 microturb', 'wave GK, v_g=c_s', 'KE density', 'wave M_x F_c, v_g=v'], 'km/s', False),
        ('alpha0', A, ['F_c/F', 'M_t F_c/F', 'F_c/F x cells', 'F_c/F x depth/R', 'F_w(GK)/F*', 'core IGW flux/F*', 'M_s [C09 microturb]'], 'norm. at 3.7', True),
        ('nuchar', N, ['1/(2pi t_c)', 'v_c/(2pi dR)', 'dynamical', 'nu_max scaling', 'core turnover'], '1/d', False)]
for ax, (obs, cands, keys, unit, norm) in zip(axes, sets):
    o = OBSSETS[obs]; hot = o['x'] > 4.3
    ob = binned(o['y'][hot], o['v'][hot])
    ref = np.interp(3.7, ob[:, 0], ob[:, 1]) if norm else 0
    ax.errorbar(ob[:, 0], ob[:, 1]-ref, yerr=[ob[:, 1]-ob[:, 2], ob[:, 3]-ob[:, 1]], fmt='ko', ms=4, label='observed', zorder=5)
    for k in keys:
        X = predict(cands[k][0], not k.startswith('core') and k not in ('dynamical', 'nu_max scaling'), o['x'][hot], o['y'][hot])
        lx = np.log10(np.where(X > 0, X, np.nan))
        b = binned(o['y'][hot], lx)
        if len(b) == 0: continue
        r = np.interp(3.7, b[:, 0], b[:, 1]) if norm else 0
        ax.plot(b[:, 0], b[:, 1]-r, '-', lw=1.3, label=k)
    if obs == 'vmac_hot':
        cs = predict(np.asarray(tm.inputs(__import__('test_transfer').sub)['c_s'])/1e5, False, o['x'][hot], o['y'][hot])
        b = binned(o['y'][hot], np.log10(cs)); ax.plot(b[:, 0], b[:, 1], 'k--', lw=1, label='surface c_s')
    ax.set_title(f'{obs} (log Teff>4.3, MS stars): log [{unit}]'); ax.set_xlabel('log L_spec')
    for x in (3.16, 3.70): ax.axvline(x, color='r', ls=':', lw=0.8)
    ax.legend(fontsize=6)
fig.tight_layout(); fig.savefig('figures/diag_binned.png', dpi=140)
print('saved')
