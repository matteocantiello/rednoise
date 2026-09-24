#!/usr/bin/env python3
"""HRDs of the rednoise MESA grid, ZAMS onward: one row per metallicity, one column per
omega/omega_crit, plus the MW MLT++ differential sub-grid overlaid on MW w0.0.

Track colour = initial mass; the MS part is drawn thicker, and the end point is marked
by status (complete / running / failed).
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import grid_io as g

W = ['w0.0', 'w0.2', 'w0.4', 'w0.6']
Z = ['MW', 'LMC', 'SMC']
norm = LogNorm(5, 120)
cmap = plt.get_cmap('plasma_r')
END = {'done': ('o', '#0ca30c'), 'running': ('>', '#2a78d6')}  # anything else = failed


def draw(ax, z, w, tracks, lw_ms=1.0, alpha=0.9, color=None, ls='-', mark_end=True):
    for m, d in sorted(tracks.items()):
        i0 = g.find_zams(d['center_h1'])
        it = g.find_tams(d['center_h1'])
        lt, ll = d['log_Teff'][i0:], d['log_L'][i0:]
        c = color or cmap(norm(m))
        n_ms = (it - i0) if it is not None else len(lt)
        ax.plot(lt[:n_ms + 1], ll[:n_ms + 1], color=c, lw=lw_ms, alpha=alpha, ls=ls)
        ax.plot(lt[n_ms:], ll[n_ms:], color=c, lw=0.45, alpha=alpha * 0.8, ls=ls)
        if mark_end:
            st = g.model_status(z, w, m)['status']
            mk, mc = END.get(st, ('x', '#d03b3b'))
            ax.plot(lt[-1], ll[-1], mk, color=mc, ms=3.5, mew=1.0, zorder=5)


fig, axes = plt.subplots(3, 4, figsize=(17, 12), sharex=True, sharey=True)
for i, z in enumerate(Z):
    for j, w in enumerate(W):
        ax = axes[i, j]
        tr = g.load_subgrid(z, w)
        draw(ax, z, w, tr)
        ax.set_title(f'{g.ZTITLE[z]}, $\\omega/\\omega_c={w[1:]}$  ({len(tr)} tracks)', fontsize=10)
        ax.grid(alpha=0.2, ls='--')
        if j == 0:
            ax.set_ylabel(r'$\log L/L_\odot$')
        if i == 2:
            ax.set_xlabel(r'$\log T_{\rm eff}$ [K]')
axes[0, 0].set_xlim(4.85, 3.45)
axes[0, 0].set_ylim(2.4, 6.7)

# end-point legend
for lab, (mk, mc) in [('complete', END['done']), ('running', END['running']), ('failed', ('x', '#d03b3b'))]:
    axes[0, 0].plot([], [], mk, color=mc, ms=5, label=lab)
axes[0, 0].plot([], [], color='0.3', lw=1.2, label='main sequence')
axes[0, 0].plot([], [], color='0.3', lw=0.45, label='post-MS')
axes[0, 0].legend(loc='lower left', fontsize=8, framealpha=0.9)

sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
cb = fig.colorbar(sm, ax=axes, pad=0.01, fraction=0.02, aspect=40)
cb.set_label(r'Initial mass [$M_\odot$]')
cb.set_ticks([5, 7, 10, 15, 20, 30, 50, 80, 120])
cb.set_ticklabels(['5', '7', '10', '15', '20', '30', '50', '80', '120'])
fig.suptitle('Rednoise MESA grid — HRD from ZAMS (Dutch winds, no MLT++)', fontsize=14, x=0.45, y=0.93)
out = f'{g.GRID}/hrd_grids.png'
fig.savefig(out, dpi=150, bbox_inches='tight', facecolor='white')
print('Saved', out)
