#!/usr/bin/env python3
"""Proposed Fig. 8: MESA FeCZ models vs observations (run compare_obs.py first)."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

HERE = os.path.dirname(os.path.abspath(__file__))
D, FIG = f'{HERE}/data', f'{HERE}/figures'
plt.rcParams.update({"font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8, "legend.fontsize": 6.5,
                     "xtick.labelsize": 7, "ytick.labelsize": 7})
C0, C2, C4 = '#2a78d6', '#eb6834', '#1baf7a'   # model w0.0, w0.2, w0.4
ONSET, ONSET_LO, ONSET_HI, SAT = 3.16, 2.94, 3.44, 3.70

ms = pd.read_csv(f'{D}/mesa_ms.csv')
onset = pd.read_csv(f'{D}/mesa_onset.csv')
ratios = pd.read_csv(f'{D}/mesa_ratios.csv')
nuc = pd.read_csv(f'{D}/mesa_nuc_offset.csv')
I = np.load(f'{D}/fig8_inputs.npz')
ob_v, ob_a = I['ob_v'], I['ob_a']


def model_binned(w, col, edges=np.arange(2.5, 4.46, 0.1), fill0=False):
    sub = ms[(ms.Z == 'MW') & (ms.w == w) & (ms.logTeff > 4.3)]
    y = sub[col].values
    if fill0:
        y = np.nan_to_num(y, nan=0.0)
    out = []
    for a, b in zip(edges[:-1], edges[1:]):
        k = sub.logLspec.between(a, b).values & np.isfinite(y)
        if k.sum() >= 20:
            out.append((0.5 * (a + b), *np.percentile(y[k], [50, 16, 84])))
    return np.array(out)


def label(ax, s):
    ax.text(-0.14, 1.03, s, transform=ax.transAxes, fontweight='bold', fontsize=9)


def hinges(ax):
    ax.axvspan(ONSET_LO, ONSET_HI, color='#e34948', alpha=0.08, lw=0)
    for x in (ONSET, SAT):
        ax.axvline(x, color='#e34948', ls=':', lw=0.9)
    x0 = onset.query("Z=='MW' and w==0.0").L_first_FeCZ.iloc[0]
    ax.axvline(x0, color='0.45', ls='--', lw=0.8)


fig, axes = plt.subplots(2, 2, figsize=(7.25, 5.6))
fig.subplots_adjust(left=0.08, right=0.98, bottom=0.08, top=0.95, wspace=0.42, hspace=0.36)

# (a) velocities
ax = axes[0, 0]
hinges(ax)
for w, c, ls in ((0.0, C0, '-'), (0.4, C4, '--')):
    mb = model_binned(w, 'FeCZ_vmax_kms', fill0=True)
    mb = mb[mb[:, 1] > 0]
    if w == 0.0:
        ax.fill_between(mb[:, 0], np.clip(mb[:, 2], 1e-2, None), mb[:, 3], color=c, alpha=0.18, lw=0)
    ax.plot(mb[:, 0], mb[:, 1], color=c, ls=ls, lw=1.5,
            label=f'model $v_{{c,\\max}}$, $\\omega/\\omega_c={w:g}$')
ax.errorbar(ob_v[:, 0], 10**ob_v[:, 1], yerr=[10**ob_v[:, 1] - 10**ob_v[:, 2], 10**ob_v[:, 3] - 10**ob_v[:, 1]],
            fmt='o', color='k', ms=3.5, lw=0.8, capsize=0, label=r'observed $v_{\rm macro}$ (median)')
ax.set_yscale('log'); ax.set_ylim(0.1, 200); ax.set_xlim(2.5, 4.45)
ax.set_xlabel(r'$\log\,\mathcal{L}/\mathcal{L}_\odot$'); ax.set_ylabel(r'velocity [km s$^{-1}$]')
ax.set_title(r'Velocities, $\log T_{\rm eff}>4.3$, MS', loc='left')
ax.text(ONSET + 0.03, 170, 'observed\nonset', color='#c03030', fontsize=5.8, va='top')
ax.text(SAT + 0.03, 170, 'observed\nsaturation', color='#c03030', fontsize=5.8, va='top')
ax.text(onset.query("Z=='MW' and w==0.0").L_first_FeCZ.iloc[0] - 0.03, 60, 'model FeCZ\nappears', color='0.35',
        fontsize=5.8, ha='right', va='top')
ax.legend(loc='lower right', frameon=False)
label(ax, 'a')

# (b) amplitude shape, normalised at log L = 3.7
ax = axes[0, 1]
hinges(ax)
ref = SAT
oa = ob_a.copy()
oa_ref = np.interp(ref, oa[:, 0], oa[:, 1])
mb = model_binned(0.0, 'FeCZ_FcF', fill0=True)
mb = mb[mb[:, 1] > 0]
lm = np.log10(mb[:, 1:])
m_ref = np.interp(ref, mb[:, 0], lm[:, 0])
ax.fill_between(mb[:, 0], lm[:, 1] - m_ref, lm[:, 2] - m_ref, color=C0, alpha=0.18, lw=0)
ax.plot(mb[:, 0], lm[:, 0] - m_ref, color=C0, lw=1.5, label=r'model $F_c/F$, $\omega/\omega_c=0$')
mv = model_binned(0.0, 'FeCZ_vmax_kms', fill0=True); mv = mv[mv[:, 1] > 0]
ax.plot(mv[:, 0], np.log10(mv[:, 1]) - np.interp(ref, mv[:, 0], np.log10(mv[:, 1])), color=C0, lw=1.0, ls=':',
        label=r'model $v_{c,\max}$')
ax.errorbar(oa[:, 0], oa[:, 1] - oa_ref, yerr=[oa[:, 1] - oa[:, 2], oa[:, 3] - oa[:, 1]], fmt='s', color='k',
            ms=3.5, lw=0.8, capsize=0, label=r'observed $\alpha_0$ (median)')
ov_ref = np.interp(ref, ob_v[:, 0], ob_v[:, 1])
ax.plot(ob_v[:, 0], ob_v[:, 1] - ov_ref, 'o', mfc='white', mec='k', ms=3.2, mew=0.8, label=r'observed $v_{\rm macro}$')
ax.set_xlim(2.5, 4.45); ax.set_ylim(-4.5, 1.5)
ax.set_xlabel(r'$\log\,\mathcal{L}/\mathcal{L}_\odot$')
ax.set_ylabel(r'$\log Q - \log Q(\log\mathcal{L}=3.7)$')
ax.set_title('Shapes, normalised at the saturation point', loc='left')
ax.legend(loc='lower right', frameon=False)
label(ax, 'b')

# (c) nu offset map
ax = axes[1, 0]
GT, GL, gp, ok, mm = I['GT'], I['GL'], I['gp_nu'], I['gp_ok'], I['model_nu_w0']
dlt = np.where(ok & np.isfinite(mm), gp - mm, np.nan)
med = nuc.query("w==0.0 and region=='OB (logT>4.3)'").median_offset_dex.iloc[0]
cf = ax.contourf(GT, GL, dlt, levels=np.linspace(med - 0.8, med + 0.8, 17), cmap='RdBu_r',
                 norm=TwoSlopeNorm(med, med - 0.8, med + 0.8), extend='both')
cs = ax.contour(GT, GL, np.where(np.isfinite(mm), mm, np.nan), levels=[-1.5, -1.0, -0.5, 0.0],
                colors='k', linewidths=0.6, linestyles='solid')
ax.clabel(cs, fmt=lambda v: f'{10**v:.2g}', fontsize=5.5)
for m in (10, 15, 20, 30, 50):
    t = ms[(ms.Z == 'MW') & (ms.w == 0.0) & (ms.Minit == m)]
    ax.plot(t.logTeff, t.logLspec, color='0.5', lw=0.5)
    ax.text(t.logTeff.iloc[0], t.logLspec.iloc[0] - 0.04, f'{m}', fontsize=5, color='0.4', ha='center', va='top')
cb = fig.colorbar(cf, ax=ax, pad=0.02, fraction=0.05)
cb.set_label(r'$\log\nu_{\rm char}^{\rm GP}-\log\nu_c^{\rm model}$ [dex]', fontsize=6.5)
cb.ax.tick_params(labelsize=6)
cb.set_ticks(np.round(np.arange(med - 0.8, med + 0.81, 0.4), 1))
ax.set_xlim(4.72, 3.9); ax.set_ylim(2.4, 4.6)
ax.set_xlabel(r'$\log T_{\rm eff}$ [K]'); ax.set_ylabel(r'$\log\,\mathcal{L}/\mathcal{L}_\odot$')
ax.set_title(r'$\nu_{\rm char}$ vs FeCZ $\nu_c$ (contours: $\nu_c$ in d$^{-1}$)', loc='left')
ax.text(0.97, 0.04, f'median offset {med:+.2f} dex (×{10**med:.1f}), $\\log T_{{\\rm eff}}>4.3$',
        transform=ax.transAxes, ha='right', fontsize=5.8, bbox=dict(fc='white', ec='none', alpha=0.8, pad=1))
label(ax, 'c')

# (d) sensitivity ratios
ax = axes[1, 1]
xs = np.arange(len(ratios))
ax.axhline(3, color='0.5', ls='--', lw=0.8)
ax.text(-0.45, 3.1, r'$F_c/F\propto\rho v_c^3$ at fixed $\rho$', fontsize=5.8, color='0.4', ha='left', va='bottom')
ax.errorbar(xs - 0.2, ratios.obs_ratio, yerr=[ratios.obs_ratio - ratios.obs_ratio_lo, ratios.obs_ratio_hi - ratios.obs_ratio],
            fmt='o', color='k', ms=5, lw=1, capsize=0, label=r'observed $\partial\log\alpha_0/\partial\log v_{\rm macro}$')
for dx, w, c in ((0.0, 0.0, C0), (0.12, 0.2, C2), (0.24, 0.4, C4)):
    r = ratios[f'w{w}_ratio']
    ax.errorbar(xs + dx, r, yerr=[r - ratios[f'w{w}_ratio_lo'], ratios[f'w{w}_ratio_hi'] - r], fmt='D', color=c,
                ms=4, lw=1, capsize=0, label=f'model $\\partial\\log(F_c/F)/\\partial\\log v_c$, $\\omega/\\omega_c={w:g}$')
ax.set_xticks(xs); ax.set_xticklabels([r'along $\log\mathcal{L}$', r'along $\tau$', r'along $\log M$'])
ax.set_xlim(-0.5, 2.6); ax.set_ylim(0, 9.5)
ax.set_ylabel('amplitude-to-velocity sensitivity ratio')
ax.set_title('Sensitivity ratios', loc='left')
ax.legend(loc='upper right', frameon=False, fontsize=5.8)
label(ax, 'd')

for ext in ('pdf', 'png'):
    fig.savefig(f'{FIG}/fig8_mesa.{ext}', dpi=200, bbox_inches='tight', facecolor='white')
print('saved', f'{FIG}/fig8_mesa.pdf')
