#!/usr/bin/env python3
"""
Sub-surface convection properties of the rednoise MESA grid mapped on the spectroscopic
HRD (and the classical HRD), in the style of Cantiello et al. 2021.

Outputs (in figures/):
  sHRD_FeCZ_<Z>_<w>.png, HRD_FeCZ_<Z>_<w>.png   8-panel figure per sub-grid
  sHRD_overview_<prop>.png                       one property, Z (rows) x omega (columns),
                                                 common colour scale
Usage:
  python3 plot_sHRD.py                 # everything
  python3 plot_sHRD.py MW w0.0         # one sub-grid only
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.spatial import cKDTree
import grid_io as g

FIGDIR = f'{g.GRID}/figures'
os.makedirs(FIGDIR, exist_ok=True)

PROPS = [  # column, label, cmap
    ('v_FeCZ_max',          r'$v_{\rm conv}^{\rm max}$ (FeCZ) [cm/s]',  'viridis'),
    ('mach_FeCZ_max',       r'Mach$_{\rm max}$ (FeCZ)',                 'inferno'),
    ('FeCZ_B_shutoff_conv', r'$B_{\rm shutoff}$ (FeCZ) [G]',            'plasma'),
    ('FeCZ_Fc_max',         r'$F_{c,\rm max}$ (FeCZ) [erg/cm$^2$/s]',   'cividis'),
    ('b_FeCZ_max',          r'$B_{\rm eq}^{\rm max}$ (FeCZ) [G]',       'magma'),
    ('turnover_FeCZ',       r'$\tau_{\rm conv}$ (FeCZ) [s]',            'YlOrRd'),
    ('v_HeII_max',          r'$v_{\rm conv}^{\rm max}$ (HeII) [cm/s]',  'viridis'),
    ('mach_HeII_max',       r'Mach$_{\rm max}$ (HeII)',                 'inferno'),
]
LABEL = {p: l for p, l, _ in PROPS}
CMAP = {p: c for p, _, c in PROPS}
OVERVIEW = ['v_FeCZ_max', 'mach_FeCZ_max', 'turnover_FeCZ', 'FeCZ_B_shutoff_conv']
LABEL_MASSES = [5, 7, 10, 15, 20, 30, 50, 80, 120]

AXES = {  # coordinate system: (y column function, y range, y label)
    'sHRD': ((1.5, 5.0), r'$\log\,\mathcal{L}/\mathcal{L}_\odot$'),
    'HRD':  ((2.0, 6.5), r'$\log\,L/L_\odot$'),
}
XR = (3.5, 4.75)
NGRID = 150
MASK_DX, MASK_DY = 0.06, 0.15   # max distance (log Teff, log L) from a model point


def ycoord(d, kind):
    return g.spec_ell(d['log_Teff'], d['log_g']) if kind == 'sHRD' else d['log_L']


def build_cloud(tracks, kind, npts=300):
    """ZAMS-onward point cloud, npts points per track spaced uniformly in arc length on the
    diagram (index-uniform sampling leaves the fast Hertzsprung-gap crossings nearly empty)."""
    x, y, props = [], [], {p: [] for p, _, _ in PROPS}
    for m, d in sorted(tracks.items()):
        i0 = g.find_zams(d['center_h1'])
        lt, yy = d['log_Teff'][i0:], ycoord(d, kind)[i0:]
        step = np.hypot(np.diff(lt) / MASK_DX, np.diff(yy) / MASK_DY)
        cum = np.r_[0, np.cumsum(np.nan_to_num(step))]
        sel = i0 + np.unique(np.minimum(np.searchsorted(cum, np.linspace(0, cum[-1], npts)), len(lt) - 1))
        x.append(d['log_Teff'][sel])
        y.append(ycoord(d, kind)[sel])
        for p in props:
            props[p].append(d[p][sel] if p in d else np.full(len(sel), np.nan))
    return np.concatenate(x), np.concatenate(y), {p: np.concatenate(v) for p, v in props.items()}


def interp(x, y, v):
    ok = np.isfinite(v) & (v > 0)   # MESA flags "no zone" with -1e99 / 0
    if ok.sum() < 50:
        return None, None
    z = np.log10(v[ok])
    GT, GL = np.meshgrid(np.linspace(*XR, NGRID), np.linspace(-1, 7, 4 * NGRID))
    zi = griddata((x[ok], y[ok]), z, (GT, GL), method='linear')
    # blank cells far from any model (linear interpolation otherwise fills the whole
    # convex hull, e.g. across the gap left by failed tracks)
    tree = cKDTree(np.c_[x[ok] / MASK_DX, y[ok] / MASK_DY])
    dist, _ = tree.query(np.c_[GT.ravel() / MASK_DX, GL.ravel() / MASK_DY])
    zi[dist.reshape(GT.shape) > 1] = np.nan
    return (GT, GL, zi), z


def overlay_tracks(ax, tracks, kind):
    for m, d in sorted(tracks.items()):
        i0 = g.find_zams(d['center_h1'])
        yy = ycoord(d, kind)
        ax.plot(d['log_Teff'][i0:], yy[i0:], color='0.3', lw=0.35, alpha=0.5, zorder=2)
        if m in LABEL_MASSES:
            ax.text(d['log_Teff'][i0], yy[i0] - 0.08, f'{int(m)}', fontsize=5, color='0.35',
                    ha='center', va='top', zorder=5)


def filled(ax, grid, levels, cmap):
    GT, GL, zi = grid
    return ax.contourf(GT, GL, zi, levels=levels, cmap=cmap, extend='both')


def log_cbar(fig, cf, ax, label, fs=7):
    cb = fig.colorbar(cf, ax=ax, pad=0.02, fraction=0.05, aspect=20)
    ticks = np.unique(np.round(cf.levels[::3], 1))
    cb.set_ticks(ticks)
    cb.set_ticklabels([f'$10^{{{t:.1f}}}$' for t in ticks])
    cb.ax.tick_params(labelsize=6)
    cb.set_label(label, fontsize=fs, labelpad=2)
    return cb


def eight_panel(z, w, tracks, kind):
    x, y, props = build_cloud(tracks, kind)
    yr, ylab = AXES[kind]
    fig, axes = plt.subplots(2, 4, figsize=(16, 8.5))
    rot = 'non-rotating' if w == 'w0.0' else f'$\\omega/\\omega_c = {w[1:]}$'
    fig.suptitle(f"Sub-surface convection on the {'spectroscopic ' if kind == 'sHRD' else ''}HRD — "
                 f'{g.ZTITLE[z]}, {rot} ({len(tracks)} tracks)' + g.GRID_TAG, fontsize=13, fontweight='bold', y=0.97)
    for j, (p, lab, cmap) in enumerate(PROPS):
        ax = axes[j // 4, j % 4]
        grid, zz = interp(x, y, props[p])
        if grid is None:
            ax.text(0.5, 0.5, 'insufficient data', transform=ax.transAxes, ha='center')
        else:
            levels = np.linspace(np.percentile(zz, 2), np.percentile(zz, 98), 15)
            log_cbar(fig, filled(ax, grid, levels, cmap), ax, lab)
        overlay_tracks(ax, tracks, kind)
        ax.set_xlim(XR[::-1]); ax.set_ylim(yr)
        ax.tick_params(labelsize=7)
        if j % 4 == 0: ax.set_ylabel(ylab, fontsize=9)
        if j >= 4: ax.set_xlabel(r'$\log T_{\rm eff}$ [K]', fontsize=9)
    plt.subplots_adjust(left=0.05, right=0.98, top=0.91, bottom=0.08, wspace=0.35, hspace=0.25)
    out = f'{FIGDIR}/{kind}_FeCZ_{z}_{w}.png'
    fig.savefig(out, dpi=170, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print('  saved', out)


def overview(all_tracks, prop, kind='sHRD'):
    """prop on the sHRD for every Z x omega, one colour scale (2-98% over all sub-grids)."""
    Z, W = ['MW', 'LMC', 'SMC'], ['w0.0', 'w0.2', 'w0.4', 'w0.6']
    clouds = {k: build_cloud(t, kind) for k, t in all_tracks.items() if k[0] in Z}
    allz = np.concatenate([np.log10(c[2][prop][np.isfinite(c[2][prop]) & (c[2][prop] > 0)])
                           for c in clouds.values()])
    levels = np.linspace(np.percentile(allz, 2), np.percentile(allz, 98), 17)
    yr, ylab = AXES[kind]
    fig, axes = plt.subplots(3, 4, figsize=(16, 11.5), sharex=True, sharey=True)
    for i, z in enumerate(Z):
        for j, w in enumerate(W):
            ax = axes[i, j]
            if (z, w) in clouds:
                x, y, props = clouds[(z, w)]
                grid, _ = interp(x, y, props[prop])
                if grid is not None:
                    cf = filled(ax, grid, levels, CMAP[prop])
                overlay_tracks(ax, all_tracks[(z, w)], kind)
            ax.set_title(f'{g.ZTITLE[z]}, $\\omega/\\omega_c={w[1:]}$', fontsize=10)
            ax.set_xlim(XR[::-1]); ax.set_ylim(yr)
            if j == 0: ax.set_ylabel(ylab)
            if i == 2: ax.set_xlabel(r'$\log T_{\rm eff}$ [K]')
    cb = fig.colorbar(cf, ax=axes, pad=0.01, fraction=0.02, aspect=40)
    ticks = np.unique(np.round(levels[::2], 1))
    cb.set_ticks(ticks); cb.set_ticklabels([f'$10^{{{t:.1f}}}$' for t in ticks])
    cb.set_label(LABEL[prop])
    fig.suptitle(f'{LABEL[prop]} across metallicity and rotation (ZAMS onward)' + g.GRID_TAG, fontsize=14, x=0.45, y=0.94)
    out = f'{FIGDIR}/{kind}_overview_{prop}.png'
    fig.savefig(out, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print('  saved', out)


if __name__ == '__main__':
    subgrids = [tuple(sys.argv[1:3])] if len(sys.argv) >= 3 else g.SUBGRIDS
    all_tracks = {}
    for z, w in subgrids:
        tr = g.load_subgrid(z, w)
        if len(tr) < 5:
            continue
        all_tracks[(z, w)] = tr
        for kind in ('sHRD', 'HRD'):
            eight_panel(z, w, tr, kind)
    if len(subgrids) > 1:
        for p in OVERVIEW:
            overview(all_tracks, p)
