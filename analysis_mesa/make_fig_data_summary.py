#!/usr/bin/env python3
"""Summary of the observational data by galaxy, with the v2 model main sequences and predicted thresholds.

(a) red-noise stars on the sHRD (colour = galaxy, marker = fitting method), MW and LMC v2 ZAMS/TAMS
(b) IACOB macroturbulence on the sHRD (Galactic only)
(c) nu_char vs Teff by galaxy and method (the Lorentzian vs GP confound)
(d) log L_spec distribution of hot (log Teff > 4.3) red-noise stars per galaxy, with the model-predicted onset
    (v_c = 3 km/s) and saturation (Gamma_Fe = 0.8) for each metallicity (regime_crossings.py)

Usage (from analysis_mesa/):  python3 make_fig_data_summary.py  -> figures_v2/fig_data_summary.{png,pdf}
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
OBS = '/mnt/home/mcantiello/work/rednoise/handoff/project_handoff/data'
D2, FIG = f'{HERE}/data_v2', f'{HERE}/figures_v2'
plt.rcParams.update({"font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8, "legend.fontsize": 6.3,
                     "xtick.labelsize": 7, "ytick.labelsize": 7})
COL = {'MW': '#2a78d6', 'LMC': '#eb6834', 'SMC': '#1baf7a'}
METHOD = {'Bowman2020_Galactic': 'Lorentzian', 'Shen2024_Galactic': 'Lorentzian', 'Bowman2019b_LMC': 'Lorentzian',
          'Bowman2024_LMC': 'GP', 'Bowman2024_SMC': 'GP', 'Ma2024_LMC_BSG': 'mod. Lorentzian',
          'DornWallenstein2020_YSG': 'damped random walk', 'DornWallenstein2020_RSG': 'damped random walk'}
GAL = {s: ('LMC' if 'LMC' in s else 'SMC' if 'SMC' in s else 'MW') for s in METHOD}
MARK = {'Lorentzian': 'o', 'GP': 's', 'mod. Lorentzian': 'D', 'damped random walk': '^'}

E = pd.read_csv(f'{OBS}/rednoise_sHRD_extended_evol.csv')
E['galaxy'] = E['sample'].map(GAL)
E['method'] = E['sample'].map(METHOD)
V = pd.read_csv(f'{OBS}/macroturbulence_evol.csv')
cross = {z: pd.read_csv(f'{D2}/regime_{z}w0.0_crossings.csv') for z in ('MW', 'LMC', 'SMC')
         if os.path.exists(f'{D2}/regime_{z}w0.0_crossings.csv')}


def ms_band(z):
    """ZAMS and TAMS lines (log Teff, log L_spec) of a v2 sub-grid extract, if present."""
    fn = f'{D2}/mesa_ms_{z}w0.0.csv'
    if not os.path.exists(fn):
        return None
    m = pd.read_csv(fn)
    m = m[m.ms_complete]
    g = m.sort_values('tau').groupby('Minit')
    zams, tams = g.first().sort_index(), g.last().sort_index()
    return zams, tams


fig, axes = plt.subplots(2, 2, figsize=(7.25, 6.4))
fig.subplots_adjust(left=0.08, right=0.98, bottom=0.08, top=0.95, wspace=0.28, hspace=0.34)


def hrd_axes(ax):
    ax.set_xlim(4.8, 3.5); ax.set_ylim(2.0, 4.75)
    ax.set_xlabel(r'$\log T_{\rm eff}$ [K]'); ax.set_ylabel(r'$\log\,\mathcal{L}/\mathcal{L}_\odot$')


def draw_ms(ax, label=True):
    for z, ls in (('MW', '-'), ('LMC', '--'), ('SMC', ':')):
        b = ms_band(z)
        if b is None:
            continue
        zams, tams = b
        for t, lw in ((zams, 0.9), (tams, 0.6)):
            ax.plot(t.logTeff, t.logLspec, color=COL[z], ls=ls, lw=lw, alpha=0.7, zorder=1)
        if label:
            ax.plot([], [], color=COL[z], ls=ls, lw=0.9, label=f'{z} model ZAMS / TAMS (v2)')


# (a) red-noise stars on the sHRD
ax = axes[0, 0]
draw_ms(ax)
for (gal, meth), d in E.groupby(['galaxy', 'method']):
    ax.scatter(d.lT, d.lL, s=12, marker=MARK[meth], color=COL[gal], ec='white', lw=0.3, zorder=3, alpha=0.9)
hrd_axes(ax)
n = E.groupby('galaxy').size()
ax.set_title(f"Red-noise stars on the sHRD: MW {n.get('MW', 0)}, LMC {n.get('LMC', 0)}, SMC {n.get('SMC', 0)}",
             loc='left')
h = [Line2D([], [], ls='', marker='o', color=COL[g], label=g) for g in ('MW', 'LMC', 'SMC')]
h += [Line2D([], [], ls='', marker=m, color='0.4', mfc='none', label=k) for k, m in MARK.items()]
leg1 = ax.legend(handles=h, loc='lower right', frameon=False, ncol=2, bbox_to_anchor=(1.0, 0.0))
ax.add_artist(leg1)
ax.legend(loc='upper right', frameon=False)

# (b) macroturbulence (Galactic only)
ax = axes[0, 1]
draw_ms(ax, label=False)
Vs = V[np.isfinite(V.vmac) & np.isfinite(V.logL_sp)]
sc = ax.scatter(Vs.logTeff_sp, Vs.logL_sp, c=np.log10(Vs.vmac), s=6, cmap='viridis', vmin=0.8, vmax=2.1,
                lw=0, zorder=3)
cb = fig.colorbar(sc, ax=ax, pad=0.02, fraction=0.05)
cb.set_label(r'$\log v_{\rm macro}$ [km s$^{-1}$]', fontsize=7); cb.ax.tick_params(labelsize=6)
hrd_axes(ax); ax.set_xlim(4.8, 3.9)
ax.set_title(f'IACOB macroturbulence: Galactic only (N = {len(Vs)})', loc='left')
ax.text(0.03, 0.04, 'no homogeneous LMC/SMC $v_{\\rm macro}$ in the compilation', transform=ax.transAxes,
        fontsize=6, color='0.35')

# (c) nu_char vs Teff
ax = axes[1, 0]
for (gal, meth), d in E.groupby(['galaxy', 'method']):
    ax.scatter(d.lT, np.log10(d.nuchar), s=12, marker=MARK[meth], color=COL[gal], ec='white', lw=0.3, alpha=0.9,
               label=f'{gal}, {meth} ({len(d)})')
ax.set_xlim(4.8, 3.5); ax.set_xlabel(r'$\log T_{\rm eff}$ [K]'); ax.set_ylabel(r'$\log\nu_{\rm char}$ [d$^{-1}$]')
ax.set_title(r'$\nu_{\rm char}$ by galaxy and fitting method', loc='left')
ax.legend(loc='lower left', frameon=False, ncol=2, fontsize=5.6)

# (d) log L distribution of hot stars, with predicted thresholds
ax = axes[1, 1]
bins = np.arange(2.2, 4.8, 0.1)
hot = E[E.lT > 4.3]
for i, z in enumerate(('MW', 'LMC', 'SMC')):
    d = hot[hot.galaxy == z]
    ax.hist(d.lL, bins=bins, histtype='step', lw=1.4, color=COL[z], label=f'{z} red-noise stars, log T_eff > 4.3 ({len(d)})')
    c = cross.get(z)
    if c is not None:
        on = c.query("parameter == 'v_c' and threshold == 3.0").logLspec.iloc[0]
        sat = c.query("parameter == 'Gamma_Fe' and threshold == 0.8").logLspec.iloc[0]
        yy = 0.93 - 0.07 * i
        ax.plot([on, sat], [yy, yy], color=COL[z], lw=1.0, transform=ax.get_xaxis_transform())
        for x, mk in ((on, 'v'), (sat, 'D')):
            ax.plot(x, yy, marker=mk, color=COL[z], ms=5, transform=ax.get_xaxis_transform())
        ax.text(sat + 0.04, yy, f'{z} model', color=COL[z], fontsize=5.8, va='center',
                transform=ax.get_xaxis_transform())
    elif z == 'SMC':
        ax.text(4.75, 0.93 - 0.07 * i, 'SMC model: grid still running', color=COL[z], fontsize=5.8, va='center',
                ha='right', transform=ax.get_xaxis_transform())
Vh = V[(V.logTeff_sp > 4.3) & np.isfinite(V.logL_sp)]
ax.hist(Vh.logL_sp, bins=bins, histtype='stepfilled', color='0.85', zorder=0,
        label=f'MW IACOB $v_{{\\rm macro}}$, log T_eff > 4.3 ({len(Vh)}) /4', weights=np.full(len(Vh), 0.25))
for x, lab in ((3.16, 'obs. onset (MW)'), (3.70, 'obs. saturation (MW)')):
    ax.axvline(x, color='0.3', ls=':', lw=0.8)
    ax.text(x - 0.02, 0.03, lab, transform=ax.get_xaxis_transform(), fontsize=5.6, va='bottom', ha='right',
            color='0.3', rotation=90, bbox=dict(fc='white', ec='none', alpha=0.8, pad=0.5))
ax.set_xlabel(r'$\log\,\mathcal{L}/\mathcal{L}_\odot$'); ax.set_ylabel('stars per 0.1 dex')
ax.set_title('Luminosity coverage vs predicted onset / saturation', loc='left')
ax.set_ylim(0, 48)
h, l = ax.get_legend_handles_labels()
h += [Line2D([], [], ls='', marker='v', color='0.4', label=r'model onset ($v_c$ = 3 km/s)'),
      Line2D([], [], ls='', marker='D', color='0.4', label=r'model saturation ($\Gamma_{\rm Fe}$ = 0.8)')]
ax.legend(handles=h, loc='center left', bbox_to_anchor=(0.0, 0.5), frameon=False, fontsize=5.6)
for a_, l_ in zip(axes.flat, 'abcd'):
    a_.text(-0.12, 1.04, l_, transform=a_.transAxes, fontweight='bold', fontsize=9)
os.makedirs(FIG, exist_ok=True)
for ext in ('png', 'pdf'):
    fig.savefig(f'{FIG}/fig_data_summary.{ext}', dpi=200, bbox_inches='tight', facecolor='white')
print('saved', f'{FIG}/fig_data_summary.png', '| model thresholds for:', list(cross))
