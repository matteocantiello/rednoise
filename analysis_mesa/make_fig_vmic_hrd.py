#!/usr/bin/env python3
"""
Microturbulence (Markova, Cantiello & Grassitelli 2025, J/A+A/701/A297; tables E1 non-pulsating + E2 pulsating) on the sHRD,
and against the subsurface convection zones of the v2 MW grid (non-rotating, MS + post-MS).

fig_vmic_shrd    v_mic on the sHRD (like panel b of fig_data_summary for v_macro), all 1823 stars; MW v2 ZAMS/TAMS and
                 the model boundaries of the FeCZ and the H-recombination zone (HI CZ).
fig_vmic_zones   row 1: model v_c,max of each subsurface zone (FeCZ, He II, He I, H I, and the largest of them) interpolated
                 between tracks, masked where the zone is absent; v_mic stars on the same colour scale.
                 row 2: v_mic against each model velocity at the stars' positions (stars on the model grid), with Spearman
                 rho, the median log ratio, and the fraction of stars where the zone is absent.
Writes figures_v2/fig_vmic_{shrd,zones}.{png,pdf} and data_obs/vmic_zones.csv.
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import vmic_test                      # noqa: E402

D2, FIG = f'{HERE}/data_v2', f'{HERE}/figures_v2'
ZONES = [('FeCZ', 'Fe bump (FeCZ)'), ('HeII', 'He II'), ('HeI', 'He I'), ('HI', 'H I'), ('MAX', 'largest of all zones')]
NORM = LogNorm(0.3, 100)
CMAP = 'viridis'
TB = np.arange(3.55, 4.76, 0.01)
LB = np.arange(2.0, 4.66, 0.025)
SX, SY = 0.03, 0.08
plt.rcParams.update({'font.size': 7.5, 'axes.titlesize': 8, 'axes.labelsize': 7.5, 'legend.fontsize': 6,
                     'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5})


def models():
    d = pd.concat([pd.read_csv(f'{D2}/mesa_ms_MWw0.0.csv').assign(phase='MS'),
                   pd.read_csv(f'{D2}/mesa_post_MWw0.0.csv').assign(phase='post')], ignore_index=True)
    d = d[np.isfinite(d.logTeff) & np.isfinite(d.logLspec)].reset_index(drop=True)
    v = np.zeros(len(d))
    for z in ('FeCZ', 'HeII', 'HeI', 'HI'):
        vz = np.where(d[f'{z}_present'].astype(bool), d[f'{z}_vmax_kms'].fillna(0).values, 0.0)
        v = np.maximum(v, vz)
    d['MAX_vmax_kms'] = v
    d['MAX_present'] = v > 0
    return d


class Field:
    """Interpolate one zone's log v_c,max between tracks; presence interpolated separately; masked off the grid."""
    def __init__(self, d, zone):
        pts = np.c_[d.logTeff / SX, d.logLspec / SY]
        self.tree = cKDTree(pts)
        pres = d[f'{zone}_present'].astype(bool).values
        v = d[f'{zone}_vmax_kms'].values.astype(float)
        ok = pres & np.isfinite(v) & (v > 0)
        self.fv = LinearNDInterpolator(pts[ok], np.log10(v[ok])) if ok.sum() > 20 else None
        self.fp = LinearNDInterpolator(pts, pres.astype(float))

    def __call__(self, lT, lL):
        P = np.c_[np.asarray(lT) / SX, np.asarray(lL) / SY]
        near = self.tree.query(P)[0] < 1.2
        p = self.fp(P)
        v = 10 ** self.fv(P) if self.fv is not None else np.full(len(P), np.nan)
        v = np.where(np.nan_to_num(p) >= 0.5, v, 0.0)          # 0 = zone absent
        return np.where(near, v, np.nan), np.where(near, p, np.nan)


def ms_lines(d):
    m = d[(d.phase == 'MS') & d.ms_complete.astype(bool)].sort_values('tau').groupby('Minit')
    return m.first().sort_index(), m.last().sort_index()


def shrd(ax, xlim=(4.72, 3.55), ylim=(-1.0, 4.5)):
    ax.set_xlim(*xlim); ax.set_ylim(*ylim)
    ax.set_xlabel(r'$\log T_{\rm eff}$ [K]'); ax.set_ylabel(r'$\log\,\mathcal{L}/\mathcal{L}_\odot$')


def fig_shrd(m, d, F):
    fig, ax = plt.subplots(1, 1, figsize=(4.8, 4.4))
    s = m[m.vmic > 0]
    z = m[m.vmic <= 0]
    for pul, mk, lab in ((False, 'o', 'non-pulsating (table E1)'), (True, 'D', 'pulsating (table E2)')):
        q = s[s.pulsating == pul]
        sc = ax.scatter(q.lT, q.lL, c=q.vmic, cmap=CMAP, norm=NORM, s=7 if not pul else 9, marker=mk, lw=0.2,
                        edgecolors='k' if pul else 'none', zorder=3, label=f'{lab}: {len(q)}')
    ax.scatter(z.lT, z.lL, s=7, facecolors='none', edgecolors='0.4', lw=0.4, zorder=2, label=f'v$_{{\\rm mic}}$ = 0: {len(z)}')
    za, ta = ms_lines(d)
    ax.plot(za.logTeff, za.logLspec, color='#2a78d6', lw=0.9, label='MW v2 ZAMS / TAMS (≥ 5 M$_\\odot$)')
    ax.plot(ta.logTeff, ta.logLspec, color='#2a78d6', lw=0.7, ls='--')
    X, Y = np.meshgrid(0.5 * (TB[1:] + TB[:-1]), 0.5 * (LB[1:] + LB[:-1]), indexing='ij')
    for zone, col, lab in (('FeCZ', 'k', 'FeCZ present (models)'), ('HI', '#b04a00', 'H I CZ present (models)')):
        _, p = F[zone](X.ravel(), Y.ravel())
        ax.contour(X, Y, np.nan_to_num(p.reshape(X.shape), nan=0), levels=[0.5], colors=col, linewidths=0.9)
        ax.plot([], [], color=col, lw=0.9, label=lab)
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=NORM, cmap=CMAP), ax=ax, pad=0.02)
    cb.set_label(r'$v_{\rm mic}$ [km s$^{-1}$]')
    shrd(ax)
    ax.set_title(f'Microturbulence, Markova+2025 (N = {len(m)})', loc='left')
    ax.legend(loc='lower left', frameon=True, fontsize=5.3)
    for ext in ('png', 'pdf'):
        fig.savefig(f'{FIG}/fig_vmic_shrd.{ext}', dpi=200, bbox_inches='tight', facecolor='white')
    plt.close(fig)


def fig_zones(m, d, F):
    fig, ax = plt.subplots(2, 5, figsize=(15, 6.4))
    fig.subplots_adjust(wspace=0.28, hspace=0.35)
    X, Y = np.meshgrid(0.5 * (TB[1:] + TB[:-1]), 0.5 * (LB[1:] + LB[:-1]), indexing='ij')
    s = m[(m.vmic > 0)].copy()
    za, ta = ms_lines(d)
    rows = []
    for j, (zone, lab) in enumerate(ZONES):
        v, p = F[zone](X.ravel(), Y.ravel())
        v = v.reshape(X.shape)
        a = ax[0, j]
        a.pcolormesh(TB, LB, np.where(v == 0, 1.0, np.nan).T, cmap='Greys', vmin=0, vmax=3, shading='flat', rasterized=True)
        a.pcolormesh(TB, LB, np.where(v > 0, v, np.nan).T, cmap=CMAP, norm=NORM, shading='flat', rasterized=True)
        a.plot(za.logTeff, za.logLspec, color='w', lw=0.8); a.plot(ta.logTeff, ta.logLspec, color='w', lw=0.5, ls='--')
        on = s[(s.lL > 2.0)]
        a.scatter(on.lT, on.lL, c=on.vmic, cmap=CMAP, norm=NORM, s=8, edgecolors='k', linewidths=0.2, zorder=4)
        shrd(a, ylim=(2.0, 4.5)); a.set_title(f'model v$_{{c,\\max}}$: {lab}', loc='left')
        if j:
            a.set_ylabel('')
        # star-by-star
        vs, ps = F[zone](s.lT.values, s.lL.values)
        g = np.isfinite(vs)
        a = ax[1, j]
        pos = g & (vs > 0)
        a.scatter(np.clip(vs[pos], 1e-4, None), s.vmic.values[pos], c=s.lT.values[pos], cmap='coolwarm', vmin=3.8, vmax=4.6,
                  s=7, lw=0)
        a.plot([1e-4, 300], [1e-4, 300], color='0.5', lw=0.7, ls='--')
        a.set_xscale('log'); a.set_yscale('log'); a.set_xlim(1e-4 if zone == 'HeI' else 0.01, 300); a.set_ylim(0.3, 50)
        a.set_xlabel(f'model v$_{{c,\\max}}$ ({lab}) [km s$^{{-1}}$]'); a.set_ylabel(r'$v_{\rm mic}$ [km s$^{-1}$]' if j == 0 else '')
        rho = spearmanr(vs[pos], s.vmic.values[pos]).correlation if pos.sum() > 5 else np.nan
        lr = np.median(np.log10(s.vmic.values[pos] / vs[pos])) if pos.any() else np.nan
        a.text(0.03, 0.97, f'N = {pos.sum()} with this zone\n{int((g & (vs == 0)).sum())} on grid without it\n'
               f'ρ = {rho:+.2f}\nmedian log(v$_{{\\rm mic}}$/v$_c$) = {lr:+.2f}', transform=a.transAxes, va='top', fontsize=6,
               bbox=dict(facecolor='white', alpha=0.8, lw=0))
        a.set_title('v$_{\\rm mic}$ vs model, star by star', loc='left', fontsize=7)
        rows.append(dict(zone=zone, n_on_grid=int(g.sum()), n_zone_present=int(pos.sum()), n_zone_absent=int((g & (vs == 0)).sum()),
                         spearman=rho, median_log_ratio=lr,
                         spearman_hot=spearmanr(vs[pos & (s.lT.values > 4.3)], s.vmic.values[pos & (s.lT.values > 4.3)]).correlation
                         if (pos & (s.lT.values > 4.3)).sum() > 5 else np.nan,
                         spearman_cool=spearmanr(vs[pos & (s.lT.values < 4.0)], s.vmic.values[pos & (s.lT.values < 4.0)]).correlation
                         if (pos & (s.lT.values < 4.0)).sum() > 5 else np.nan,
                         n_cool=int((pos & (s.lT.values < 4.0)).sum())))
    cax = fig.add_axes([0.915, 0.56, 0.008, 0.32])
    fig.colorbar(plt.cm.ScalarMappable(norm=NORM, cmap=CMAP), cax=cax, label='velocity [km s$^{-1}$]')
    cax2 = fig.add_axes([0.915, 0.11, 0.008, 0.32])
    fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(3.8, 4.6), cmap='coolwarm'), cax=cax2, label=r'star $\log T_{\rm eff}$')
    for ext in ('png', 'pdf'):
        fig.savefig(f'{FIG}/fig_vmic_zones.{ext}', dpi=180, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    return pd.DataFrame(rows)


def main():
    m = vmic_test.load()
    d = models()
    F = {z: Field(d, z) for z, _ in ZONES}
    fig_shrd(m, d, F)
    r = fig_zones(m, d, F)
    r.to_csv(f'{HERE}/data_obs/vmic_zones.csv', index=False, float_format='%.3f')
    pd.set_option('display.width', 200)
    print(r.to_string(index=False, float_format=lambda v: f'{v:+.2f}'))


if __name__ == '__main__':
    main()
