"""Regime parameters of the FeCZ vs spectroscopic luminosity (hot MS models), with the observed hinges.

Panels: (a) Eddington factor at the Fe peak Gamma_Fe; (b) max MLT velocity in the FeCZ and its Mach number;
(c) tau_Fe / tau_crit; (d) pseudo-Mach number Y_Fe. Median and 16-84% band over hot (log Teff > 4.3) MS
models in 0.1-dex bins of log L_spec; the 3D Athena++ models (Schultz+2023, Table 1) are overplotted in (a)
and (c). Observed onset 3.16 [2.94, 3.44] (shaded) and v_macro saturation 3.70 (dotted).

Usage (from analysis_mesa/): python3 make_fig_regime.py [data_v2/regime_MWw0.0.csv] [figures_v2]
"""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

f = sys.argv[1] if len(sys.argv) > 1 else 'data_v2/regime_MWw0.0.csv'
out = sys.argv[2] if len(sys.argv) > 2 else 'figures_v2'
C = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100']
INK, INK2, MUTED, GRID, AXIS = '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#c3c2b7'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': AXIS, 'axes.labelcolor': INK2, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6,
                     'axes.spines.top': False, 'axes.spines.right': False, 'lines.linewidth': 1.6,
                     'legend.frameon': False, 'legend.fontsize': 7.5})

d = pd.read_csv(f)
ms = d[(d.Xc > 1e-3) & (d.Xc < 0.7) & (d.logTeff > 4.3)].copy()
ms['bin'] = np.floor(ms.logLspec / 0.1) * 0.1 + 0.05
grp = ms.groupby('bin')
keep = grp.size() >= 5
x = np.array(keep.index[keep])

# 3D Athena++ models (Schultz+2023 Table 1): log L_spec from Teff (flux-defined) and the core mass/radius
# are not tabulated, so they are placed by Gamma_Fe and tau_Fe/tau_crit only, in panels (a) and (c), at the
# log L_spec of the matching MESA position (hot MS models only: M13TAMS, M35ZAMS, M35MMS).
three_d = {'M13 TAMS': dict(Gamma=0.55, tr=0.24), 'M35 ZAMS': dict(Gamma=0.82, tr=0.02),
           'M35 MMS': dict(Gamma=0.97, tr=0.25)}


def band(ax, col, color, label, ls='-'):
    q = grp[col].quantile([0.16, 0.5, 0.84]).unstack().loc[x]
    ax.fill_between(x, q[0.16], q[0.84], color=color, alpha=0.15, lw=0)
    ax.plot(x, q[0.5], color=color, ls=ls, label=label)
    return q


fig, axs = plt.subplots(2, 2, figsize=(7.4, 5.6), sharex=True)
for ax in axs.flat:
    ax.axvspan(2.94, 3.44, color='#ebebe7', lw=0, zorder=0)
    ax.axvline(3.16, color=MUTED, lw=0.8)
    ax.axvline(3.70, color=INK2, lw=0.9, ls=(0, (1, 2)))
ax = axs[0, 0]
q = band(ax, 'Gamma_Fe', C[0], r'$\Gamma_{\rm Fe}=\kappa_{\rm Fe}L/(4\pi GMc)$')
ax.axhline(1, color=MUTED, lw=0.8)
ax.set_ylabel(r'$\Gamma$ at the Fe peak')
ax.text(3.72, 0.12, r'$v_{\rm macro}$ saturation', rotation=90, fontsize=7.5, color=INK2, va='bottom')
ax.text(3.18, 0.62, 'observed onset', rotation=90, fontsize=7.5, color=INK2, va='bottom')
for i, (nm, v) in enumerate(three_d.items()):
    ax.axhline(v['Gamma'], xmin=0.93, xmax=1.0, color=C[1], lw=2)
    ax.text(4.32, v['Gamma'], nm, fontsize=6.5, color=INK2, va='center', ha='right')
ax.set_ylim(0, 1.4)
ax.legend(loc='upper left')
ax = axs[0, 1]
band(ax, 'v_c', C[0], r'$v_{c,\max}$ [km s$^{-1}$]')
band(ax, 'mach_max', C[1], r'Mach $v_{c}/c_s$', ls='--')
ax.set_yscale('log')
ax.set_ylabel('MLT velocity, Mach number')
ax.legend(loc='upper left')
ax = axs[1, 0]
band(ax, 'tau_ratio', C[0], r'$\tau_{\rm Fe}/\tau_{\rm crit}$ (MESA)')
for nm, v in three_d.items():
    ax.axhline(v['tr'], xmin=0.93, xmax=1.0, color=C[1], lw=2)
ax.text(4.32, 0.25 * 1.15, 'M13 TAMS, M35 MMS', fontsize=6.5, color=INK2, va='bottom', ha='right')
ax.text(4.32, 0.02 * 1.15, 'M35 ZAMS', fontsize=6.5, color=INK2, va='bottom', ha='right')
ax.set_yscale('log')
ax.set_ylim(5e-3, 3)
ax.axhline(1, color=MUTED, lw=0.8)
ax.set_ylabel(r'$\tau_{\rm Fe}/\tau_{\rm crit}$')
ax.set_xlabel(r'$\log\,\mathcal{L}/\mathcal{L}_\odot$')
ax.legend(loc='upper left')
ax = axs[1, 1]
band(ax, 'Y_Fe', C[0], r'pseudo-Mach $\mathcal{Y}_{\rm Fe}$')
ax.set_yscale('log')
ax.axhline(1, color=MUTED, lw=0.8)
ax.set_ylabel(r'$\mathcal{Y}$ at the Fe peak')
ax.set_xlabel(r'$\log\,\mathcal{L}/\mathcal{L}_\odot$')
ax.legend(loc='upper left')
for ax in axs.flat:
    ax.set_xlim(2.5, 4.35)
fig.suptitle('FeCZ regime parameters, MESA v2 (MW, non-rotating, hot MS); 3D models of Schultz+2023 at right',
             fontsize=8.5, color=INK)
fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(f'{out}/fig_regime.{ext}', dpi=170, bbox_inches='tight')
print(f'wrote {out}/fig_regime.png/.pdf')
